import os
import json
import logging
import argparse
from datetime import datetime
import time
import re
import pandas as pd
import numpy as np
import os
import sys
import platform
import json 

def setup_logger(file_name: str, log_folder: str) -> logging.Logger:

    # Ensure the log folder exists
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    
    # Create a logger instance with a unique name to prevent conflicts
    logger = logging.getLogger(f'logger_{file_name}')
    logger.setLevel(logging.INFO)
    
    # Define the full path for the log file
    log_file_path = os.path.join(log_folder, f'{file_name}.log')

    # Create a file handler to write logs to the specified file in append mode
    handler = logging.FileHandler(log_file_path, mode='w')  
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    # Set the log format
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    
    return logger

def read_json_config(json_file_path: str) -> dict:

    with open(json_file_path, 'r') as file:
        config = json.load(file)

    return config

def get_reaction_ids_to_process_old(reaction_id_list: str, 
                                synthons_folder: str) -> list[str]:
  
    
    all_reaction_ids = [
        os.path.basename(x).replace('.txt','') for x in os.listdir(synthons_folder)
    ]
    
    if reaction_id_list == -1: # processing all reactions
        return all_reaction_ids
    elif isinstance(reaction_id_list, list) and len(reaction_id_list) > 0:
        if set(reaction_id_list).issubset(set(all_reaction_ids)):
            return reaction_id_list
        else:
            raise ValueError("Reaction ids are not present in the input")
    else:
        raise ValueError("The reaction_id has to be a list of ids or -1.")

def get_reaction_ids_to_process(reaction_id_list_path: str, synthons_folder: str) -> list[str]:
   
    if reaction_id_list_path: 

        # Read reaction IDs from CSV
        reaction_id_df = pd.read_csv(reaction_id_list_path)
        
        # Check if the expected column exists
        if 'reaction_id' not in reaction_id_df.columns:
            raise ValueError("CSV file must contain a 'reaction_id' column.")
        
        reaction_id_list = reaction_id_df['reaction_id'].astype(str).tolist()

        # Get all reaction IDs from the synthons folder
        all_reaction_ids = {
            os.path.basename(x).replace('.txt', '') for x in os.listdir(synthons_folder)
        }

        # Return only reaction IDs that are present in the synthons folder
        return [rid for rid in reaction_id_list if rid in all_reaction_ids]

    else:
    # Get all reaction IDs from the synthons folder
            all_reaction_ids = {
                os.path.basename(x).replace('.txt', '') for x in os.listdir(synthons_folder)
            }
            return list(all_reaction_ids)
    
def log_script_info(logger, config):
    """
    Logs system and Python information, experiment configuration, and launch time.
    """
    launch_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  
    logger.info('-'*150)
    logger.info('Input INFO')
    logger.info(f"Script name: {sys.argv[0]}")
    logger.info(f"Script started at {launch_time}")
    logger.info(f"Command-line arguments: {sys.argv}")  

    logger.info(f"Configuration: {json.dumps(config, indent=4)}")
    if config.get('filter_config_path', False):
        logger.info(f"Filters Configuration: {json.dumps(read_json_config(config['filter_config_path']), indent=4)}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"System version: {platform.platform()}")
    logger.info(f"Machine architecture: {platform.machine()}")
    logger.info(f"Python executable: {sys.executable}")
    logger.info(f"Current working directory: {os.getcwd()}")
    logger.info(f"Number of CPUs: {config['N_cores']}/{os.cpu_count()}")
    logger.info('-'*150)
    #logger.info(f"Environment variables: {json.dumps(dict(os.environ), indent=4)}")

def calculate_full_space_size(synthons_df: pd.DataFrame, column_name='synton#') -> int:

    space_size = (
        synthons_df.groupby('reaction_id')[column_name]
        .apply(lambda x: x.value_counts().prod()).sort_values(ascending=False)
    )
    return space_size

def setup_report_logger(output_folder='DEFAULt', file_name='Report',  mode: str = 'w') -> logging.Logger:
    """
    Sets up a logger to log script execution details with time format only.
    """
    # Ensure the log folder exists
    if not os.path.exists(output_folder):
         os.makedirs(output_folder, exist_ok=True)
    
    logger = logging.getLogger(f'{file_name}')
    logger.setLevel(logging.INFO)
    
    log_file_path = os.path.join(output_folder, f'{file_name}.log')
    handler = logging.FileHandler(log_file_path, mode=mode)  
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    formatter = logging.Formatter('- %(message)s')
    handler.setFormatter(formatter)
    
    return logger

def mel_id_generator(full_synthons_IDs_list):
    for full_synthon_ID in full_synthons_IDs_list:
        mel_IDs = full_synthon_ID.split('|')
        for mel_ID in mel_IDs:
            yield mel_ID


def get_reaction_category(filename: str, categories: list) -> dict:
    """
    Determine the reaction category based on the filename.
    """
    filename = filename.lower()
    for category in categories:
        if category in filename:
            return category
    return None