import os
from collections import ChainMap
import pandas as pd
import requests


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

############################################
# Part 1 - Client to query mychem.info API #
############################################

def query_drugbank_name(drugbank_id):
    """
    Find the drugbank name of the input `drugbank_id`
    See https://docs.mychem.info/en/latest/doc/chem_annotation_service.html#get-request for the specification of the API

    Args:
        drugbank_id (str): a drugbank ID like "DB00002"

    Returns:
        drugbank_name (str): a drugbank name like 'Cetuximab' or None if not found
    """
    url = 'http://mychem.info/v1/chem/{}'.format(drugbank_id)
    params = {"fields": "drugbank.name"}
    response = requests.get(url, params=params)  # GET method required

    # raise an HTTPError if the HTTP request returned an unsuccessful status code
    response.raise_for_status()

    json_response = response.json()

    """
    E.g. the json response of http://mychem.info/v1/chem/DB00002?fields=drugbank.name is

    ```
    {
        "_id": "DB00002",
        "_version": 1,
        "drugbank": {"_license": "http://bit.ly/2PSfZTD", "name": "Cetuximab"}
    }
    ```

    The json response of http://mychem.info/v1/chem/DB12430?fields=drugbank.name is

    ```
    {"code": 404, "success": false, "error": "ID 'DB12430' not found"}
    ```
    """
    if "drugbank" not in json_response:
        return None

    drugbank_name = json_response["drugbank"]["name"]
    return drugbank_name

def _query_drugbank_names(drugbank_ids):
    """
    Find the drugbank names of the input `drugbank_ids`
    See https://docs.mychem.info/en/latest/doc/chem_annotation_service.html#batch-queries-via-post for the specification of the API
    This is a private function that ignores the payload limit (1000 IDs at most in one query)

    Args:
        drugbank_ids (list): a list of drugbank IDs like ["DB00002", "DB01234"]

    Returns:
        id_name_map (dict): a mapping of drugbank ids and names like {'DB00002': 'Cetuximab', 'DB01234': 'Dexamethasone'}
    """
    url = 'http://mychem.info/v1/chem'
    data = {"ids": ",".join(drugbank_ids), "fields": "drugbank.name"}
    response = requests.post(url, data=data)  # POST method required

    # raise an HTTPError if the HTTP request returned an unsuccessful status code
    response.raise_for_status()

    json_response = response.json()

    """
    E.g. the json response when querying `["DB00002", "DB12430"]` is

    ```
    [
        {
            'query': 'DB00002',
            '_id': 'DB00002',
            '_version': 1,
            'drugbank': {'_license': 'http://bit.ly/2PSfZTD', 'name': 'Cetuximab'}
        },
        {
            'query': 'DB12430',
            'notfound': True
        }
    ]
    ```
    """
    id_name_map = {entry["query"]: (entry["drugbank"]["name"] if "drugbank" in entry else None) 
                   for entry in json_response}

    return id_name_map

def batch_query_drugbank_names(drugbank_ids, batch_size=1000):
    """
    Find the drugbank names of the input `drugbank_ids` by batches.
    See https://docs.mychem.info/en/latest/doc/chem_annotation_service.html#batch-queries-via-post for the specification of the API
    The API has a payload limit that at most 1000 IDs can be posted in one query.
    Therefore if `len(drugbank_ids)` is greater than 1000, this function will partition the ID list into batches, 
    call the API multiple times, and merge the results.

    Args:
        drugbank_ids (list): a list of drugbank IDs like ["DB00002", "DB01234"]

    Returns:
        id_name_map (dict): a mapping of drugbank ids and names like {'DB00002': 'Cetuximab', 'DB01234': 'Dexamethasone'}
    """
    if batch_size is None or len(drugbank_ids) <= batch_size:
        return _query_drugbank_names(drugbank_ids)

    drugbank_id_batches = [drugbank_ids[i:i+batch_size] for i in range(0, len(drugbank_ids), batch_size)]
    id_name_map_batches = [_query_drugbank_names(id_batch) for id_batch in drugbank_id_batches]

    id_name_map = dict(ChainMap(*id_name_map_batches))
    return id_name_map

##########################################
# Part 2 - Util to revise the repoDB csv #
##########################################

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
    may correspond to multiple drug names. E.g. drugbank ID "DB00002" is found in the following entries:

    | drug_name               | drugbank_id | ind_name                           | ind_id   | NCT         | status     | phase           | DetailedStatus      |
    |-------------------------|-------------|------------------------------------|----------|-------------|------------|-----------------|---------------------|
    | cetuximab               | DB00002     | Squamous cell carcinoma of mouth   | C0585362 | NA          | Approved   | NA              | NA                  |
    | cetuximab               | DB00002     | Squamous cell carcinoma of nose    | C3163899 | NA          | Approved   | NA              | NA                  |
    | cetuximab               | DB00002     | Squamous cell carcinoma of pharynx | C1319317 | NA          | Approved   | NA              | NA                  |
    | cetuximab               | DB00002     | Laryngeal Squamous Cell Carcinoma  | C0280324 | NA          | Approved   | NA              | NA                  |
    | cetuximab               | DB00002     | Malignant tumor of colon           | C0007102 | NA          | Approved   | NA              | NA                  |
    | NA                      | DB00002     | Non-Small Cell Lung Carcinoma      | C0007131 | NCT00203931 | Terminated | Phase 2         | Slow accrual and ...|
    | dexamethasone phosphate | DB00002     | Multiple Myeloma                   | C0026764 | NCT00368121 | Terminated | Phase 2         | Lack of ...         |
    | cetuximab               | DB00002     | Non-Small Cell Lung Carcinoma      | C0007131 | NCT00694603 | Terminated | Phase 2         | Slow accrual ...    |
    | cetuximab               | DB00002     | Squamous cell carcinoma            | C0007137 | NCT01794845 | Terminated | Phase 2         | Early termination...|
    | NA                      | DB00002     | Squamous cell carcinoma            | C0007137 | NCT02298595 | Withdrawn  | Phase 1/Phase 2 | NA                  |

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

    assert is_one_to_one(repoDB_df, "drug_name", "drugbank_id"), "drug_name and drugbank_id are not 1-to-1 after manipulation"

    return repoDB_df

def is_one_to_one(df, col1, col2):
    """
    Check if two columns are one-to-one in a dataframe.

    A one-to-one relation between two columns is like (same with the bijection in functions):

    | col1 | col2      |
    |------|-----------|
    | a    | apple     |
    | b    | banana    |
    | c    | cranberry |

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

#####################################
# Part 3 - Parser of the repoDB csv #
#####################################

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
