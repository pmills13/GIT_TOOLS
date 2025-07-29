import pandas as pd
import matplotlib.pyplot as plt
import itertools
import numpy as np
import os

class DataFramePlotter:
    def __init__(self, df: pd.DataFrame, saveDir = None):
        self.df = df

        ## if saveDir is passed, save it and set the save flag
        self.saveDir = saveDir
        self.isSave = False
        if not saveDir == None:
          self.isSave = True

        # Define at least 5 colors, markers, and line styles
        self.colors = ['black','b', 'g', 'r', 'c', 'm']  # At least 5 colors
        self.markers = ['o', 's', 'D', '^', 'v']  # At least 5 markers
        self.line_styles = ['solid','-', '--', '-.', ':']  # At least 5 line styles

        # Create a product of all possible combinations of colors, markers, and line styles
        self.extended_styles = list(itertools.product(self.line_styles,self.markers,self.colors))

        # Pre-assign styles to each run (Run_ID) in a structured order (non-random)
        self.run_styles = {}

    def assign_styles(self):
        """
        Assigns a unique color, marker, and line style to each Run_ID in a structured manner.
        Loops over color first, then marker, and then line style.
        """
        unique_runs = self.df['Run'].unique()
        for i, run_id in enumerate(unique_runs):
            # Cycle through the extended styles if there are more runs than styles
            self.run_styles[run_id] = self.extended_styles[i % len(self.extended_styles)]

    def assign_styles_components(self, dep_vars):
        """
        Assigns a unique color, marker, and line style to each Run_ID in a structured manner.
        Loops over color first, then marker, and then line style.
        """
        for i, var_id in enumerate(dep_vars):
            # Cycle through the extended styles if there are more runs than styles
            self.run_styles[var_id] = self.extended_styles[i % len(self.extended_styles)]

    def filter_data(self, filter_conditions: dict):
        """
        Filters the DataFrame based on provided filter conditions.
        :param filter_conditions: A dictionary where keys are column names and values are conditions.
        """
        for column, condition in filter_conditions.items():
            if isinstance(condition, list):
                self.df = self.df[self.df[column].isin(condition)]
            else:
                self.df = self.df[self.df[column] == condition]

    def construct_legend_label(self, run_data: pd.DataFrame, metadata_columns: list):
        """
        Constructs a custom legend label using the provided metadata columns.
        :param run_data: A DataFrame slice for a specific run.
        :param metadata_columns: List of DataFrame columns to be used in the legend.
        :return: A formatted legend string.
        """
        if not metadata_columns:
            # Default to Run_ID and Mach if no metadata columns are provided
            return f'Run {run_data["Run_ID"].iloc[0]} (Mach {run_data["Mach"].iloc[0]})'

        # Combine the metadata into a formatted string
        metadata_values = [f'{col}: {run_data[col].iloc[0]}' for col in metadata_columns]
        return ', '.join(metadata_values)

    def plot_sweep(self, dep_vars: list, independent_vars=None, metadata_columns=None, title: str = None, grid=True, dep_vars_rename = None, independent_vars_rename = None):
        """
        Plots dependent variables against independent variables.
        If independent_vars is a single variable, all dep_vars are plotted against it.
        If independent_vars is a list, each dep_var is plotted against the corresponding independent variable.

        :param dep_vars: List of dependent variable column names (e.g., CL, CD, CM).
        :param independent_vars: Single independent variable or list of independent variables.
        :param metadata_columns: List of column names to be used in the legend.
        :param title: Title of the entire plot (optional).
        :param grid: Option to include grid on plots (default is True).
        """
        if independent_vars is None:
            raise ValueError("You must specify at least one independent variable.")

        # If a single independent variable is provided, convert it to a list to make handling easier
        if not isinstance(independent_vars, list):
            independent_vars = [independent_vars] * len(dep_vars)

            if independent_vars_rename != None:
                independent_vars_rename = [independent_vars_rename] * len(dep_vars)

        if len(independent_vars) != len(dep_vars):
            raise ValueError("The number of independent variables must match the number of dependent variables or be a single variable.")

        unique_runs = self.df['Run'].unique()

        num_vars = len(dep_vars)
        num_rows = (num_vars + 1) // 2

        fig, axes = plt.subplots(num_rows, 2, figsize=(14, 5 * num_rows),dpi=150)
        axes = axes.flatten() if num_vars > 1 else [axes]

        # Ensure styles are assigned before plotting
        self.assign_styles()

        for i, dep_var in enumerate(dep_vars):
            independent_var = independent_vars[i]
            ax = axes[i]
            for run_id in unique_runs:
                run_data = self.df[self.df['Run'] == run_id]

                # Get the style combination for the current run (from pre-assigned styles)
                style, marker, color = self.run_styles[run_id]

                # Construct the legend label using the provided metadata columns
                label = self.construct_legend_label(run_data, metadata_columns)

                ax.plot(run_data[independent_var], run_data[dep_var],
                        linestyle=style, marker=marker,
                        label=label,
                        color=color)  # Ensure color is applied consistently

            #print(independent_vars_rename)
            if independent_vars_rename == None:
                ax.set_xlabel(independent_var)
            else:
                ax.set_xlabel(independent_vars_rename[i])

            if dep_vars_rename == None:
                ax.set_ylabel(dep_var)
            else:
                ax.set_ylabel(dep_vars_rename[i])

            ax.legend()#loc='upper left')
            ax.grid(grid)

        # Remove unused axes if the number of dependent variables is odd
        if num_vars % 2 != 0:
            fig.delaxes(axes[-1])

        # Set the overall plot title
        title = ''
        for i_dep_var in dep_vars:
          ## when saving a file, '/' is an illegal character. replace it with 'o'
          if '/' in i_dep_var:
            i_dep_var = i_dep_var.replace('/','o')
          title += i_dep_var + '_'
        title += 'VS_' + independent_vars[0]
        if title:
            fig.suptitle(title, fontsize=16)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        if self.isSave == True:
          fig.savefig(os.path.join(self.saveDir, title + ".png"))
          plt.close(fig)
        else:
          plt.show()


    def plot_sweep_components(self,run_id,dep_vars_list, independent_var,metadata_columns = None, grid=True, body_length_m = np.nan, x_CG_m=np.nan):
        """
        Plots dependent variables (grouped into subplots) against a single independent variable.

        :param dep_vars_list: List of lists, where each sublist contains dependent variable names to be plotted together.
        :param independent_var: The independent variable name.
        :param data: DataFrame containing the variables.
        :param title: Title of the entire plot (optional).
        :param grid: Option to include grid on plots (default is True).
        """
        # Number of subplots needed
        num_subplots = len(dep_vars_list)
        num_rows = (num_subplots + 1) // 2

        # Create subplots

        fig, axes = plt.subplots(num_rows, 2, figsize=(14, 5 * num_rows),dpi=150)
        axes = axes.flatten() if num_subplots > 1 else [axes]

        # Handle case with a single subplot
        if num_subplots == 1:
            axes = [axes]

        # Loop through each group of dependent variables

        data = self.df[self.df['Run'] == run_id]


        for ax, dep_vars in zip(axes, dep_vars_list):
            self.assign_styles_components(dep_vars)

            for dep_var in dep_vars:
                # Plot each dependent variable against the independent variable
                style, marker, color = self.run_styles[dep_var]

                if dep_var.split('_')[0] in ['CXb','CYb','CZb','Clb','Cmb','Cnb','FXb','FYb','FZb','Mxb','Myb','Mzb']\
                    or dep_var.split('_')[0] in ['CXs','CYs','CZs','Cls','Cms','Cns','FXs','FYs','FZs','Mxs','Mys','Mzs','CLs','CDs']:
                        ylabel=dep_var.split('_')[0]
                        label=dep_var.split('_')[1]

                elif dep_var == 'x_NP_m':
                    label=dep_var
                    ylabel='x_m'
                else:
                    ylabel=dep_var
                    label=dep_var

                ax.plot(data[independent_var], data[dep_var], label=label,
                        linestyle=style, marker=marker,
                        color=color)

                if dep_var == 'SM_m':
                    ax.hlines(0.1*body_length_m,np.min(data[independent_var]),np.max(data[independent_var]),linestyle='--',label='stability target',color='green')
                if dep_var == 'x_NP_m':
                    ax.hlines(x_CG_m+0.1*body_length_m,np.min(data[independent_var]),np.max(data[independent_var]),linestyle='--',label='stability NP target',color='green')
                    ax.hlines(x_CG_m,np.min(data[independent_var]),np.max(data[independent_var]),linestyle='--',label='x_CG_m',color='blue')
                    ax.plot(data[independent_var],np.array(data[dep_var])-0.1*body_length_m,linestyle='--',label='stability CG target',color='orange')


            # Add labels, legend, and grid
            ax.set_xlabel(independent_var)
            ax.set_ylabel(ylabel)

            ax.legend()
            if grid:
                ax.grid(True)

        # Set the overall title
        if num_subplots % 2 != 0:
            fig.delaxes(axes[-1])

        title = self.construct_legend_label(data, metadata_columns)

        if metadata_columns:
            fig.suptitle(title, fontsize=16)

        # Adjust layout
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()

# %% Example usage:
# df = pd.read_csv('sample_data.csv')  # Load the CFD AoA sweep data
# plotter = DataFramePlotter(df)

# # Filter the data to only include specific Mach numbers or runs (optional)
# filter_conditions = {'Mach': 0.8, 'Run_ID': [1, 2]}
# plotter.filter_data(filter_conditions)

# # Example 1: Plot all dependent variables against a single independent variable ('AoA'), with custom metadata for legends
# metadata_columns = ['Run_ID', 'Mach', 'AoA']  # You can change this to include any columns you want in the legend
# plotter.plot_sweep(['CL', 'CD', 'CM', 'CL', 'CM', 'CD'], independent_vars='AoA', metadata_columns=metadata_columns, title='Aerodynamic Coefficients vs AoA')

# # Example 2: Plot each dependent variable against a different independent variable, with custom metadata
# plotter.plot_sweep(['CL', 'CD', 'CM'], independent_vars=['AoA', 'Mach', 'Run_ID'], metadata_columns=['Mach', 'Run_ID'], title='Different Independent Variables per Dependent')
