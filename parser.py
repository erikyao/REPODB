import os
import pandas as pd

"""
This script parses entries of the repoDB csv file (named `full.csv`) and 
outputs documents as SmartAPI requires.

Each document represent a unique drug using drugbank id as primary key (`_id`).
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
                "name": "Hypertensive disease ",
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

        # Some entries will have a line break following 4 spaces inside its `DetailedStatus` field; 
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
        A RepoDB doc is composed from a drug entry and a list of indication entries.

        Args:
            drug_entry (DrugEntry): a DrugEntry object
            indication_entries (list): a list of IndicationEntry objects 
        """
        self.drug = drug_entry
        self.indications = indication_entries

    def to_dict(self):
        """Represent the RepoDB document as a dictionary."""
        repodb_dict = {
            "drugbank": self.drug.id,
            "name": self.drug.name,
            "indications": [vars(ind) for ind in self.indications]  # `vars(obj)` is equivalent to `obj.__dict__`
        }

        ret_dict = {
            "_id": self.drug.id,
            "repodb": repodb_dict
        }

        return ret_dict

def load_data(data_folder):
    repoDB_file = os.path.join(data_folder, "full.csv")
    repoDB_df = pd.read_csv(repoDB_file, na_filter=False)

    for drug_tuple, indication_dataframe in repoDB_df.groupby(["drugbank_id", "drug_name"], as_index=False):
        # each group key is transformed into a drug entry
        drug_entry = DrugEntry(*drug_tuple)
        # each indication inside a certain group is transformed into an Indication Entry
        indication_entries = [IndicationEntry(series) for _, series in indication_dataframe.iterrows()]

        repoDB_doc = RepoDBDoc(drug_entry, indication_entries)
        yield repoDB_doc.to_dict()
