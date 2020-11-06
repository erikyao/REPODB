import requests
from collections import ChainMap

def query_drugbank_name(drugbank_id):
    """
    Find the drugbank name of the input `drugbank_id`
    See https://docs.mychem.info/en/latest/doc/chem_annotation_service.html#get-request for the specification of the API

    Args:
        drugbank_id (str): a drugbank ID like "DB00002"

    Returns:
        drugbank_name (str): a drugbank name like 'Cetuximab'
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
