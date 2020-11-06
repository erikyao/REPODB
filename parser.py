import os
import pandas as pd
from MyChemInfoClient import batch_query_drugbank_names

"""
This script parses entries of the repoDB csv file (named `full.csv`) and 
outputs documents as SmartAPI requires.

Each document represents a unique drug using drugbank id as primary key (`_id`).
One drug may have multiple indications, and these indications will be grouped into one document.

See https://github.com/SmartAPI/smartAPI/issues/85 for more details.

Example: Suppose we have the following entries in the repoDB csv.

| drug_name     | drugbank_id | ind_name                                           | ind_id   | NCT | status   | phase | DetailedStatus |
|---------------|-------------|----------------------------------------------------|----------|-----|----------|-------|----------------|
| enalapril     | DB00584     | Asymptomatic left ventricular systolic dysfunction | C3698411 | NA  | Approved | NA    | NA             |
| enalapril     | DB00584     | Hypertensive disease                               | C0020538 | NA  | Approved | NA    | NA             |

The expected output is:

```
{
    "_id": "DB00584",
    "repodb": {
        "drugbank": "DB00584",
        "name": "enalapril",
        "indications": [
            {
                "name": "Asymptomatic left ventricular systolic dysfunction",
                "umls": "C3698411",
                "NCT": "NA",
                "status": "Approved",
                "phase": "NA"
                "detailed_status": "NA"
            },
            {
                "name": "Hypertensive disease",
                "umls": "C0020538",
                "NCT": "NA",
                "status": "Approved",
                "phase": "NA"
                "detailed_status": "NA"
            }
        ]
    }
}
```
"""

class IndicationEntry:
    def __init__(self, series):
        """Init an indication entry from a pandas Series object."""
        self.name = series.ind_name
        self.umls = series.ind_id
        self.NCT = series.NCT
        self.status = series.status
        self.phase = series.phase

        # Some entries will have a line break followed by 4 spaces inside its `DetailedStatus` field
        # Here we replace such a sequence with a single space
        self.detailed_status = series.DetailedStatus.replace("\n    ", " ")

class DrugEntry:
    def __init__(self, drugbank_id, drug_name):
        """Init a drug entry with 2 fields, `drugbank_id` and `drug_name`."""
        self.id = drugbank_id
        self.name = drug_name

class RepoDBDoc:
    def __init__(self, drug_entry, indication_entries):
        """
        A RepoDB doc is composed of a drug entry and a list of indication entries.

        Args:
            drug_entry (DrugEntry): a DrugEntry object
            indication_entries (list): a list of IndicationEntry objects
        """
        self.drug = drug_entry
        self.indications = indication_entries

    def to_dict(self):
        """Represent the RepoDB document as a dictionary."""
        ret_dict = {
            "_id": self.drug.id,
            "repodb": {
                "drugbank": self.drug.id,
                "name": self.drug.name,
                "indications": [vars(ind) for ind in self.indications]  # `vars(obj)` is equivalent to `obj.__dict__`
            }
        }

        return ret_dict

def revise_drugbank_name(repoDB_df):
    """
    Revise the drugbank name of the original repoDB csv.

    Args:
        repoDB_df (pandas.DataFrame): the original dataframe from the repoDB csv

    Returns:
        repoDB_df (pandas.DataFrame): the revised dataframe
    """

    """
    Addendum (2020-11-05): we found a possible glitch in the original repoDB csv file that one drugbank ID
    may correspond to multiple drug names. E.g. drugbank ID "DB00002" is found on the following entries:

    | drug_name               | drugbank_id | ind_name                           | ind_id   | NCT         | status     | phase           | DetailedStatus                   |
    |-------------------------|-------------|------------------------------------|----------|-------------|------------|-----------------|----------------------------------|
    | cetuximab               | DB00002     | Squamous cell carcinoma of mouth   | C0585362 | NA          | Approved   | NA              | NA                               |
    | cetuximab               | DB00002     | Squamous cell carcinoma of nose    | C3163899 | NA          | Approved   | NA              | NA                               |
    | cetuximab               | DB00002     | Squamous cell carcinoma of pharynx | C1319317 | NA          | Approved   | NA              | NA                               |
    | cetuximab               | DB00002     | Laryngeal Squamous Cell Carcinoma  | C0280324 | NA          | Approved   | NA              | NA                               |
    | cetuximab               | DB00002     | Malignant tumor of colon           | C0007102 | NA          | Approved   | NA              | NA                               |
    | NA                      | DB00002     | Non-Small Cell Lung Carcinoma      | C0007131 | NCT00203931 | Terminated | Phase 2         | Slow accrual and evidence...     |
    | dexamethasone phosphate | DB00002     | Multiple Myeloma                   | C0026764 | NCT00368121 | Terminated | Phase 2         | Lack of recruitable patients     |
    | cetuximab               | DB00002     | Non-Small Cell Lung Carcinoma      | C0007131 | NCT00694603 | Terminated | Phase 2         | Slow accrual                     |
    | cetuximab               | DB00002     | Squamous cell carcinoma            | C0007137 | NCT01794845 | Terminated | Phase 2         | Early termination due to lack... |
    | NA                      | DB00002     | Squamous cell carcinoma            | C0007137 | NCT02298595 | Withdrawn  | Phase 1/Phase 2 | NA                               |

    It's not clear how NA was produced in the 1st column but for the 7th entry (dexamethasone phosphate) 
    Kevin found it's actually a combo drug with cetuximab. 

    To maintain the one-to-one relationship between `drug_name` and `drugbank_id`, here we apply the following policy:

    - Discard the original `drug_name` column in the repoDB csv file
    - Use mychem.info API to find a correct drug name for each drugbank ID
    """

    """
    Addendum (2020-11-06): we found that some drug bank IDs are not found in mychem.info (e.g. "DB12430")

    To maintain the one-to-one relationship between `drug_name` and `drugbank_id`, here we apply the following policy:

    - Temporarily discard the original `drug_name` column in the repoDB csv file
    - Use mychem.info API to find a drug name for each drugbank ID
        - If not found, use the original `drug_name`
            - If the original `drug_name` is NA or not unique, manipulate manually
    """

    drugbank_ids = repoDB_df.drugbank_id.unique()
    id_name_map = batch_query_drugbank_names(drugbank_ids)

    for _, row in repoDB_df.iterrows():
        if row.drugbank_id in id_name_map:
            row.drug_name = id_name_map[row.drugbank_id]

    assert is_one_to_one(repoDB_df, "drug_name", "drugbank_id"), "drugname and drugbank_id are not 1-to-1 after manipulation"

    return repoDB_df

def is_one_to_one(df, col1, col2):
    """
    Check if two columns are one-to-one in a dataframe.

    A one-to-one relation between two columns is like:

    | col1 | col2      |
    |------|-----------|
    | a    | apple     |
    | b    | banana    |
    | c    | cranberry |

    Same with the bijection in functions.

    Args:
        df (pandas.DataFrame): a pandas dataframe
        col1 (str): name of the 1st column
        col2 (str): name of the 2nd column
    Returns:
        result (boolean): True if one-to-one otherwise False
    """
    df = df.drop_duplicates(subset=[col1, col2])
    is_injective = (df.groupby(col1)[col2].count().max() == 1)
    is_surjective = (df.groupby(col2)[col1].count().max() == 1)

    return is_injective and is_surjective

def load_data(data_folder):
    repoDB_file = os.path.join(data_folder, "full.csv")

    # "NA" strings in the csv will be preserved instead of being converted to `np.nan`
    repoDB_df = pd.read_csv(repoDB_file, na_filter=False)

    # Revise the drugbank name of the original repoDB csv.
    repoDB_df = revise_drugbank_name(repoDB_df)

    for drug_tuple, indication_dataframe in repoDB_df.groupby(["drugbank_id", "drug_name"], as_index=False):
        # each group key is transformed into a drug entry
        drug_entry = DrugEntry(*drug_tuple)
        # each indication inside a certain group is transformed into an Indication Entry
        indication_entries = [IndicationEntry(series) for _, series in indication_dataframe.iterrows()]

        repoDB_doc = RepoDBDoc(drug_entry, indication_entries)
        yield repoDB_doc.to_dict()
