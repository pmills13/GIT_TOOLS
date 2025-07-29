# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 05:20:07 2024

@author: mfund
"""

import pandas as pd
import os

class DataFrameMerger:
    def __init__(self, file_path: str):
        """
        Initializes the DataFrameMerger by reading the original dataframe from an Excel file.
        :param file_path: Path to the Excel file containing the original dataframe.
        """
        self.file_path = file_path
        self.df = pd.read_excel(file_path)

    def append_new_runs(self, new_run_df: pd.DataFrame):
        """
        Appends new runs to the original dataframe, ensuring no duplicate Run entries are added.
        :param new_run_df: DataFrame containing the new runs to append.
        """
        # Identify the Runs in both dataframes
        original_run_ids = set(self.df['Run'])
        new_run_ids = set(new_run_df['Run'])

        # Identify the non-duplicate Runs in the new dataframe
        non_duplicate_runs = new_run_ids - original_run_ids

        if non_duplicate_runs:
            # Filter out only the non-duplicate runs from the new_run_df
            filtered_new_run_df = new_run_df[new_run_df['Run'].isin(non_duplicate_runs)]

            # Append the new runs to the original dataframe
            self.df = pd.concat([self.df, filtered_new_run_df], ignore_index=True)
            print(f"\vAdded {len(filtered_new_run_df)} new runs.\n")
        else:
            print("\nNo new runs were added. All 'Run' entries already exist.\n")

    def save_to_excel(self):
        """
        Saves the updated dataframe back to the original Excel file.
        """
        with pd.ExcelWriter(self.file_path, engine='openpyxl', mode='w') as writer:
            self.df.to_excel(writer, index=False)
        print(f"Updated dataframe saved to {self.file_path}")

    def get_dataframe(self):
        """
        Returns the updated dataframe after appending new runs.
        :return: DataFrame with the new runs appended.
        """
        return self.df

# Function to handle master.xlsx
def update_master(new_df, file_path='master.xlsx'):
    # Check if the master file exists
    if os.path.exists(file_path):
        # Load the existing master file
        master_df = pd.read_excel(file_path)

        # Identify rows in new_df that have "Run" numbers not present in master_df
        #non_duplicate_df = new_df[~new_df['Run'].isin(master_df['Run'])]

        ##PCM: temporary fix is to toggle this line, to force merge results
        non_duplicate_df = new_df[new_df['Run'].isin(master_df['Run'])]

        # Append the non-duplicate rows to the master DataFrame
        #updated_master_df = pd.concat([master_df, non_duplicate_df])
        ##PCM: This should combine the results in order
        updated_master_df = pd.concat([master_df, non_duplicate_df], ignore_index = True)

        # Save the updated master DataFrame
        updated_master_df.to_excel(file_path, index=False)
        print('\nUpdated ' + str(file_path.split('\\')[-1]) +' with non-duplicate data.\n')
    else:
        # If master.xlsx doesn't exist, write new_df to it
        new_df.to_excel(file_path, index=False)
        print('\n' + str(file_path.split('\\')[-1]) + ' created with the new data.\n')