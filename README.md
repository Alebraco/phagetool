# phagetool
**An evolution-proof tool that identifies phage cocktails to treat bacterial infections**

The code is divided into two main parts: Data Collection and Bacteria.

### 1. Data Collection
This part of the code involves defining functions to fetch data from the NCBI database. It includes:

- **Protein ID Retrieval** ([protein_ids.py](https://github.com/Alebraco/phagetool/blob/main/data_collection/protein_ids.py)):
    - Function: `retrieve_ids`
    - Description: Fetches IDs from an NCBI database based on a given query.
- **Protein Names** ([protein_names.py](https://github.com/Alebraco/phagetool/blob/main/data_collection/protein_names.py)):
    - Functions: `retrieve_titles`, `fix_unnamed`
    - Description: Retrieves protein names and accession numbers for given IDs. Also renames unnamed proteins (empty strings) for consistency.
- **Sequence Retrieval** ([sequences.py](https://github.com/Alebraco/phagetool/blob/main/data_collection/sequences.py)):
    - Function: `fetch_sequences`
    - Description: Fetches amino acid sequences for given accession numbers.
- **Receptors Function** ([receptors_function.py](https://github.com/Alebraco/phagetool/blob/main/data_collection/receptors_function.py)):
    - Function: `receptors`
    - Description: Combines previously defined functions to retrieve receptor proteins names and sequences of a specific pathogenic species.
- **Pathogen Processing** ([pathogen_processing.py](https://github.com/Alebraco/phagetool/blob/main/data_collection/pathogen_processing.py)):
    - Function: `read_pathogens`, `query_pathogens`
    - Description: Reads a list of pathogenic species and obtains information on their receptor proteins.
- **Title Processing** ([title_processing.py](https://github.com/Alebraco/phagetool/blob/main/data_collection/title_processing.py)):
    - Functions: `unique_titles`, `protdict`
    - Description: Organize protein data for analysis, based on unique titles or amino acid sequences.

### 2. Pathogenic Bacteria
This section is focused on processing information related to pathogenic bacteria. It includes:

- **Downloading Pathogenic Bacteria List**: Fetching a list of pathogenic bacteria from Barlett et al. (2022).
  
- **Clean and Prepare Bacteria List**: Removing unwanted characters from the bacteria list and combining genus and species names.

- **Phage Hosts Data Collection**: `phageid_to_host` fetches information about phages, including their IDs and host.

- **Pathogenic Hosts Selection**: `select_hosts` compares phage hosts to the pathogenic bacteria list and keeps the information of phages with pathogenic hosts.

### Utils 
This directory contains utility functions that mainly support the Data Collection part.

- **User Interaction** ([data_management.py](https://github.com/Alebraco/phagetool/blob/main/utils/data_management.py)):
    - Description: Manages user interactions, such as prompting for the user email for NCBI Entrez.
- **Data Management** ([user_interaction.py](https://github.com/Alebraco/phagetool/blob/main/utils/user_interaction.py)):
    - Description: Functions for reading and storing data in JSON format.
  

This script is a tool for extracting, processing, and organizing biological data from the NCBI database, focusing on pathogenic bacteria and their associated protein receptors. It utilizes the Entrez programming utilities from BioPython and employs various Python libraries for data handling and processing.
