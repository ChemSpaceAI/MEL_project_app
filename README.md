# MEL Enumeration Pipeline

This project provides a multi-step pipeline for generating, processing, and enumerating MEL data from custom synthons and reactions

## 🧩 Project Structure
```
MEL_Project_app/
├── mel_package/                  
│   ├── __init__.py
│   ├── 1_generate.py              
│   ├── 2_process_generation.py              
│   ├── 3_enumerate.py             
│   └── functions/              
│       └── ...                  
│
├── configs/                   
│   └── Config_For_*.json
|   └── Filter_Config.json
│
├── data/                       
│   ├── FINAL_data/
│   ├── Generated_MEL_on_FINAL_data/
│   ├── Test_Input_1/
│   └── Test_Output_1/
│
├── scripts/                  
│   └── run_generate.sh
│   └── run_process_generation.sh
│   └── run_enumeration.sh
│
├── notebooks/                 
│   └── generate_test_sdf.ipynb
│
├── requirements.txt
├── README.md

```
---

## 🧪 1. Generate MELs

Takes a reaction table and synthon libraries to generate MELs.

### Comments 


### 🧾 JSON Config Example
```json
{
    "synthons_folder_path":"data/SYNTHONS_by_reaction_id_modified_for_MEL",
    "reaction_file_path": "data/REACTION_file_for_MEL.tsv",
    "external_minimal_synthons_path":  "data/Minimal_Caps_by_reaction_id_for_generation_MEL",

    "output_folder_path": "data/Test_Output/Generation_Test",

    "n_processes":30,
    "custom_reaction_id": false
}
```
### Run
```bash
bash scripts/run_enumeration.sh 
```
---
## 🧹 2. Process MELs

Processes generated MEL files into manageable chunks by MEL type.


### Comments 

### 🧾 JSON Config Example (Config_For_Processing.json)

```json
{
  "input_dir": "data/Generated_MEL_on_FINAL_data/enumerated_MEL",
  "output_dir": "data/Generated_MEL_on_FINAL_data/processed_enumerated_MEL",
  "chunk_size": 40000,
  "mel_types": ["3_component", "2_component", "bridge"]
}
```
### Run
```bash
bash scripts/run_process_generation.sh configs/Config_For_Processing.json
```
---
## 🔄 3. Enumerate MEL

Fully enumerates final molecules from MEL structures.


### Comments 
mel_type can be : "2_component", "3_component", "bridge"
iteration_level: 1, 2
sdf to proces sshoud have mel_synthon_id column

if not filter ocnfig type false

### 🧾 JSON Config Example (Config_For_Enumeration.json)
```json
{    
    "reactions_file_path": "data/REACTION_file_for_MEL.tsv",
    "syntons_folder_path": "data/SYNTHONS_by_reaction_id_modified_for_MEL",
    "minimal_caps_folder_path": "data/Minimal_Caps_by_reaction_id_for_enumeration_MEL",

    "filter_config_file_path": "configs/Filter_Config.json",

    "output_folder_path": "data/Test_Output/2_comp",
    "path_to_sdf_to_process": "data/Test_Input/2_comp.sdf",

    "mel_type": "2_component",
    "iteration_level": 1,

    "N_cores":50,
    "N_compounds_to_Enumerate": 5000000
}
```

### Run
```bash
bash scripts/run_enumeration.sh
```