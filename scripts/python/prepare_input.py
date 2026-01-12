import requests
import os
import argparse

# Function to download the top AF2 database entry's .pdb and .json files for a given UniProt ID
def download_af2_model(uniprot_id, out_dir):
    """
    Downloads an AlphaFold2 (AF2) model by querying the AlphaFold DB API using the UniProt ID,
    and saves the corresponding PDB file and PAE .json file directly to out_dir using the filenames in the urls.

    Args:
        uniprot_id (str): The UniProt ID (e.g., 'Q8NB16').
        out_dir (str): The directory to save the downloaded PDB and JSON file.

    Returns:
        str or None: The URL to the downloaded PDB file or None if not found or download failed.
    """

    # Query the AlphaFold API to get pdbUrl and paeDocUrl
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        api_response = requests.get(api_url)
        if api_response.status_code != 200:
            print(f"Failed to query AlphaFold API for {uniprot_id} (HTTP {api_response.status_code})")
            return None
        api_json = api_response.json()
        if not api_json or 'pdbUrl' not in api_json[0] or 'paeDocUrl' not in api_json[0]:
            print(f"No pdbUrl or paeDocUrl found in AlphaFold API response for {uniprot_id}")
            return None
        pdb_url = api_json[0]['pdbUrl']
        paeDocUrl = api_json[0]['paeDocUrl']
    except Exception as e:
        print(f"Error fetching/parsing AlphaFold API response: {e}")
        return None

    try:
        os.makedirs(out_dir, exist_ok=True)

        # Download the PDB file
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code != 200:
            print(f"Failed to download PDB file for {uniprot_id} at {pdb_url} (HTTP {pdb_response.status_code})")
            return None
        pdb_filename = f"{uniprot_id}.pdb"
        save_path_pdb = os.path.join(out_dir, pdb_filename)
        with open(save_path_pdb, "wb") as f:
            f.write(pdb_response.content)

        # Download the PAE JSON file
        pae_filename = f"{uniprot_id}.json"
        pae_path = os.path.join(out_dir, pae_filename)
        pae_response = requests.get(paeDocUrl)
        if pae_response.status_code != 200:
            print(f"Failed to download PAE JSON file for {uniprot_id} at {paeDocUrl} (HTTP {pae_response.status_code})")
            return None
        with open(pae_path, "wb") as f:
            f.write(pae_response.content)

        return pdb_filename
    except Exception as e:
        print(f"Failed to download or save PDB or PAE file for {uniprot_id}: {e}")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare DPAM-AI input files.")
    parser.add_argument("--list_path", help="Path to the _struc.list file containing the UniProt IDs of the models to download", required=True)

    args = parser.parse_args()

    # Path to the _struc.list file containing the UniProt IDs of the models to download
    list_path = args.list_path
    list_path = os.path.abspath(list_path)

    # Path to save the downloaded PDB and PAE JSON files
    out_dir = os.path.dirname(list_path)
    list_name = str.replace(list_path, "_struc.list", "")

    download_dir = os.path.join(out_dir, list_name)
    os.makedirs(download_dir, exist_ok=True)

    # Download the PDB and PAE JSON files for each UniProt ID in the _struc.list file
    with open(list_path, "r") as f:
        for line in f:
            uniprot_id = line.strip()
            pdb_url = download_af2_model(uniprot_id, download_dir)

            # If pdb_url is None, remove the uniprot_id from the list
            if pdb_url is None:
                print(f"Failed to download PDB file for {uniprot_id}, removing from list")
                with open(list_path, "r") as f:
                    lines = f.readlines()
                with open(list_path, "w") as f:
                    for line in lines:
                        if line.strip() != uniprot_id:
                            f.write(line)