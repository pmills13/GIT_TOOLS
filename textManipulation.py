# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 06:55:25 2024

@author: mfund
"""

import numpy as np

def read_fluent_out_to_dict(filename):
    """
    Reads a text file and returns a dictionary.
    
    The file should have field names on the third line in the format:
    ("Iteration" "quantity 1" "quantity 2" "quantity 3")
    
    The values are space-delimited starting from the fourth line.
    
    Args:
        filename (str): The path to the input text file.
    
    Returns:
        dict: A dictionary where keys are field names and values are numpy arrays.
    """
    
    with open(filename, 'r') as file:
        # Skip the first two lines
        next(file)
        next(file)
        
        # Read the third line and extract the field names
        field_names_line = next(file).strip()
        
        # Remove parentheses and split by spaces to get field names
        field_names = field_names_line.replace('(', '').replace(')', '').replace('"', '').split()
        
        # Read the rest of the file (data values)
        data = np.loadtxt(file)
    
    # Create a dictionary with field names as keys and columns of data as numpy arrays
    data_dict = {field_names[i]: data[:, i] for i in range(len(field_names))}
    
    return data_dict

import csv
import os

def read_csvs_to_dict(var1, var2, folder):
    nested_dict = {}

    # Iterate over the files from 1 to 6
    for i in range(1, 7):
        filename = f"{var1}_{var2}_{i}.csv"
        file_path = os.path.join(folder, filename)

        # Ensure the file exists before proceeding
        if not os.path.exists(file_path):
            print(f"File {file_path} does not exist. Skipping...")
            continue
        
        with open(file_path, mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            
            # Initialize lists for var1 and var2 for the current file
            var1_list = []
            var2_list = []
            
            # For each row in the CSV, add values to var1 and var2 lists
            for row in csv_reader:
                var1_list.append(float(row['x']))
                var2_list.append(float(row[' y']))
        
        # Add data to the dictionary with the key being the integer i
        nested_dict[i] = {var1: var1_list, var2: var2_list}
    
    return nested_dict

def read_fluent_xy_to_dict(filepath, var1, var2):
    nested_dict = {}
    current_block = None  # To keep track of the current data block (e.g., "x_d_1")

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()

            # Skip empty lines or titles/labels
            if not line or line.startswith('(title') or line.startswith('(labels'):
                continue

            # Detect new data block (header) and extract the key
            if line.startswith('((xy/key/label'):
                # Extract block key (e.g., "x_d_1" from '((xy/key/label "x_d_1")')
                current_block = line.split('"')[1]
                nested_dict[current_block] = {var1: [], var2: []}
                continue

            # Process valid data lines within a block
            if current_block:
                values = line.split()  # Split values by tab or spaces
                if len(values) == 2:  # Ensure exactly two values are present
                    try:
                        # Attempt to convert both values to float, skipping non-numeric lines
                        value1 = float(values[0])
                        value2 = float(values[1])
                        nested_dict[current_block][var1].append(value1)  # First column as var1
                        nested_dict[current_block][var2].append(value2)  # Second column as var2
                    except ValueError:
                        # Skip any lines that have non-numeric values
                        continue

    return nested_dict

import re

def extract_log_residuals(log_file_path, case_files, residual_keys):
    """
    Extracts the last iteration residuals for each case in a Fluent log file.

    Parameters:
        log_file_path (str): Path to the Fluent output log file.
        case_files (list): List of case file names (without extensions) to extract residuals for.
        residual_keys (list): List of residual names corresponding to values in the log.

    Returns:
        dict: A dictionary with case names as keys and dictionaries of residuals as values.
    """
    # Initialize a dictionary to store residuals for each case
    case_residuals = {case: None for case in case_files}

    # Regular expressions
    case_pattern = r'file/write-case-data\s+"(.+)"'  # Matches the case name when writing the case
    residual_pattern = (
        r"\[\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z\]:\s+(\d+)" + 
        r"\s+([\deE.+-]+)" * len(residual_keys)  # Adjust to match the number of residuals
    )

    # Read the file into memory for processing
    with open(log_file_path, "r") as file:
        lines = file.readlines()[::-1]  # Reverse for backward search

    # Process each case individually
    for case_name in case_files:
        found_residuals = None
        found_case = False

        for i_line in range(len(lines)):
            # Match residuals first

            # Match case name from "file/write-case-data"
            case_match = re.search(case_pattern, lines[i_line])
            if case_match:
                found_case_name = case_match.group(1).strip()
                
                
                if found_case_name == case_name:
                
                    for i_lines_below in range(0,100):
                        residual_match = re.search(residual_pattern, lines[i_line+i_lines_below])
                        
                        if residual_match:
                            residual_values = residual_match.groups()
                            residuals = dict(zip(["iteration"] + residual_keys, residual_values))
                            case_residuals[case_name] = residuals  # Store residuals
                            found_residuals=True
                            break  # Stop searching for this case
    
                    if found_residuals:
                        break

    return case_residuals


def extract_log_residuals_test(log_file_path, case_files, residual_keys):
    """
    Extracts the last iteration residuals for each case in a Fluent log file.

    Parameters:
        log_file_path (str): Path to the Fluent output log file.
        case_files (list): List of case file names (without extensions) to extract residuals for.
        residual_keys (list): List of residual names corresponding to values in the log.

    Returns:
        dict: A dictionary with case names as keys and dictionaries of residuals as values.
    """
    # Initialize a dictionary to store residuals for each case
    case_residuals = {case: None for case in case_files}

    # Regular expressions
    case_pattern = r'file/write-case-data\s+"(.+)"'  # Matches the case name when writing the case
    
    # Two possible patterns, with or without timestamp
    residual_pattern_with_timestamp = (
        r"\[\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z\]:\s+(\d+)" +
        r"\s+([\deE.+\-]+)" * len(residual_keys)  # Adjust to match the number of residuals
    )
    residual_pattern_without_timestamp = (
        r"(\d+)" + r"\s+([\deE.+\-]+)" * len(residual_keys)
    )

    # Read the file into memory for processing
    with open(log_file_path, "r") as file:
        lines = file.readlines()[::-1]  # Reverse for backward search

    # Process each case individually
    for case_name in case_files:
        found_residuals = None

        for i_line in range(len(lines)):
            # Match case name from "file/write-case-data"
            case_match = re.search(case_pattern, lines[i_line])
            if case_match:
                found_case_name = case_match.group(1).strip()
                
                if found_case_name == case_name:
                    for i_lines_below in range(0, 100):
                        line_to_check = lines[i_line + i_lines_below]
                        
                        # Try both patterns
                        residual_match = re.search(residual_pattern_with_timestamp, line_to_check)
                        if not residual_match:
                            residual_match = re.search(residual_pattern_without_timestamp, line_to_check)
                        
                        if residual_match:
                            residual_values = residual_match.groups()
                            residuals = dict(zip(["iteration"] + residual_keys, residual_values))
                            case_residuals[case_name] = residuals  # Store residuals
                            found_residuals = True
                            break  # Stop searching for this case
                    
                    if found_residuals:
                        break

    return case_residuals
