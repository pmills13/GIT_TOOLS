# -*- coding: utf-8 -*-
"""
@author: pmills
"""

import sys
import os
import subprocess
import argparse
import shutil

pathToEngineeringTools = 'G:\\My Drive\\tools'
if pathToEngineeringTools not in sys.path:
    sys.path.append(pathToEngineeringTools)

from airspeed_v1 import Airspeed
from calculateAtmosphere_v1 import *
from unitConverter import *
import TUI_scripting as tui
from textManipulation import *
from DataFramePlotter import DataFramePlotter
from DataFrameMerger import update_master
from rotationMatrices import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
from scipy.interpolate import interp1d
import json

## setup the argparser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("-case", help="Identifying string for the SpaceClaim file. If no case is entered, The code will select the most recently created mesh workflow file directory in the current directory.\nNOTE: you do not need to include \'workflow_files\' in the case name.\nUsage: -case Specific_Case_Name", type=str, required=False)
parser.add_argument("-forces", help="Optional flag that will report raw Force values in addition to Coefficients. (-forces  (no other input needed)) If not entered, defaults to Coefficients only.\nUsage: -forces", action="store_true", required=False)
parser.add_argument("-model_factor", help="Model factor for the case. \nIf it is a halfspan model, enter 0.5.\nIf it is a fullspan model, enter 1.0.\nIf It is a quarterspan model, enter 0.25.\nIf not entered, defaults to halfspan.\nUsage: -model_factor 0.5", type=float, required = False)

global x_MRC_m, y_MRC_m, z_MRC_m
global L_ref_m, S_ref_m2, model_factor, isBackPressured

###############################
## READ_INPUTS
###############################
def Read_Inputs(currDir):
  global model_factor

  args = parser.parse_args()
  model_factor = args.model_factor
  recordForces = args.forces
  case = args.case

  if model_factor == None:
    model_factor = 0.5

  if case == None:
    print('>> The case flag was not set, searching the current directory for the most recent mesh workflow folder...')
    mesh_files = []
    for file in os.listdir(currDir):
      temp = os.path.join(currDir, file)
      if temp.split('.')[-1] == 'msh':
        mesh_files.append(file)

    if len(mesh_files) == 0:
      print('>>>> ERROR... no .msh files were found in the currrent directory(\''+str(currDir)+'\'). Exiting now')
      sys.exit(0)
    else:
      case = max(mesh_files, key=os.path.getmtime)
      case_path = os.path.join(currDir, case)
      case = case.replace('.msh', '')
      print('>>>> The case file to be used is:\t'+str(case))
  else:
    caseFound = False
    for dir in os.listdir(currDir):
      temp = os.path.join(currDir, dir)
      if os.path.isdir(temp) and str(case) in temp:
        case_path = temp
        print('\n>> The case you entered was found. The code will grab data from within ' +str(case_path))
        caseFound = True

        if '_workflow_files' in case:
          case = case.replace('_workflow_files','')
        break

    if caseFound == False:
      print('!! ERROR. The case provided was not found in the name of any of the mesh workflow file directories in the current directory. Exiting now.')
      sys.exit(0)

  return(case, case_path, recordForces)

###############################
## READ_JOURNAL_OUTPUTS
###############################

def Read_Journal_Outputs(currDir):
  global x_MRC_m, y_MRC_m, z_MRC_m
  global L_ref_m, S_ref_m2, model_factor, isBackPressured

  ## check that the OUTPUTS directory was created and the 3 required files are in it
  metaFileFound = False
  caseFileFound = False
  surfFileFound = False
  for dir in os.listdir(currDir):
    tempDir1 = os.path.join(currDir, dir)
    if os.path.isdir(tempDir1):
      if '_POST' in tempDir1:
        tempDir = os.path.join(tempDir1, 'case_files')
        if os.path.exists(tempDir):
          for file in os.listdir(tempDir):
            if 'METADATA_FOR_POST' in file:
              metaFileFound = True
              postFile_Meta = os.path.join(tempDir, file)
            elif 'CASE_DATA_FOR_POST' in file:
              caseFileFound = True
              postFile_Case = os.path.join(tempDir, file)
            elif 'SURFACES_FOR_POST' in file:
              surfFileFound = True
              postFile_Surf = os.path.join(tempDir, file)
          if metaFileFound == True and caseFileFound == True and surfFileFound == True:
            outputDir = tempDir1
            break
          else:
            if metaFileFound == False:
              print('!! ERROR !! The metadata file from journal_setup.py could not be found... Exiting now.')
              sys.exit(0)
            elif caseFileFound == False:
              print('!! ERROR !! The case file from journal_setup.py could not be found... Exiting now.')
              sys.exit(0)
            elif surfFileFound == False:
              print('!! ERROR !! The surface file from journal_setup.py could not be found... Exiting now.')
              sys.exit(0)

  ## read the METADATA file
  with open(postFile_Meta, 'r') as fin:
    try:
      metadata_dict = json.load(fin)
    except Exception as error:
      print('!! ERROR !! The file ' + str(postFile_Meta) + ' could not be read... Exiting now.')
      sys.exit(0)
  fin.close()

  ## read the case file
  caseData = []
  with open(postFile_Case, 'r') as fin:
    try:
      for line in fin:
        caseData.append(line.strip().split())
    except Exception as error:
      print('!! ERROR !! The file ' + str(postFile_Case) + ' could not be read... Exiting now.')
      sys.exit(0)
  fin.close()

  case_name_list = []
  for entry in range(len(caseData)):
    header = caseData[entry][0]
    match header:
      case 'Base_case:':
        base_case = caseData[entry][1]
      case 'Journal_name:':
        journal_name = caseData[entry][1]
      case 'Case_name_list:':
        for case_name in range(1, len(caseData[entry])):
          case_name_list.append(caseData[entry][case_name])
      case 'MRC_m:':
        x_MRC_m = float(caseData[entry][1])
        y_MRC_m = float(caseData[entry][2])
        z_MRC_m = float(caseData[entry][3])
      case 'L_ref:':
        L_ref_m = float(caseData[entry][1])
      case 'S_ref:':
        S_ref_m2 = float(caseData[entry][1])
      case 'model_factor:':
        model_factor = float(caseData[entry][1])
      case 'isBackPressured:':
        isBackPressured = bool(caseData[entry][1])
      case _:
        print('>> Warning... Unexpected header (' + str(header) + ') in ' +str(postFile_Case) +'. Skipping this entry.')

  ## read the surface list file
  surfData = []
  with open(postFile_Surf, 'r') as fin:
    try:
      for line in fin:
        surfData.append(line.strip().split())
    except Exception as error:
      print('!! ERROR !! The file ' + str(postFile_Surf) + ' could not be read... Exiting now.')
      sys.exit(0)
  fin.close()

  surface_list = surfData[0][1:]
  interface_list = surfData[1][1:]

  return(base_case, journal_name, case_name_list, metadata_dict, surface_list, interface_list, outputDir)

###############################
## CREATE_WALL_OUTPUTS
###############################
def Create_Wall_Outputs(currDir, results_fileDir, base_case, case_name_list, metadata_dict, surface_list, interface_surfs, recordForces):

  global x_MRC_m, y_MRC_m, z_MRC_m
  global L_ref_m, S_ref_m2
  global model_factor
  global isBackPressured

  master_dict = {}

  its_to_average_start = -1 # all cases converged on NS residuals

  wall_groups = []
  wall_group_names = []
  for i_surface in range(len(surface_list)):
    wall_groups.append([surface_list[i_surface]])
    wall_group_names.append(surface_list[i_surface])

  wall_groups.append([wall for wall in surface_list if 'wall_aero' in wall])
  wall_groups.append([wall for wall in surface_list if 'wall_prop' in wall])
  wall_groups.append([wall for wall in surface_list if 'wall_prop' in wall or 'wall_aero' in wall])

  wall_group_names.append('wall_aero')
  wall_group_names.append('wall_prop')
  wall_group_names.append('wall_total')

  keys_to_ignore = ['FXb','FYb','FZb','MXb','MYb','MZb','MYprime','int_'] #don't write keys containing these strings to output
  ##PCM: Shouldnt we be dividing the mdot by the model factor as well ???????
  keys_to_scale = ['fx_','fy_','fz_','mx_','my_','mz_'] #all outputs containing these strings will be divided by model factor to arrive at fullspan load
  ##PCM: NEED TO FIND WAY TO SCALE !ONLY! THE APPROPRIATE STRUTS
  keys_not_to_scale = ['canard','fin', 'strake'] #outputs containing these strings won't be scaled (use for individual canards, fins, etc.)

  keys_to_correct = keys_to_scale
  quants_to_sum = keys_to_scale

  dim_keys = ['FX_','FY_','FZ_','MX_','MY_','MZ_'] #restoring correct cases
  nondim_keys = ['CX_','CY_','CZ_','Cl_','Cm_','Cn_']

  for i_case in range(len(case_name_list)):
    outputFile = os.path.join(results_fileDir, case_name_list[i_case] + '_reports.out' )
    if os.path.exists(outputFile):
      try:
        output_dict = read_fluent_out_to_dict(outputFile)
      except Exception as error:
        print('!! ERROR !! While parsing the file ' + str(outputFile) + ' the following error occured:\n')
        print(error)
        continue
    else:
      print('\n>> Warning... The file: ' +str(outputFile) + ' could not be found. Skipping this entry.\n')
      continue

    #first append the metadata you wrote to previously
    results_dict = {}
    results_dict['Base case file'] = base_case
    results_dict['Output file'] = outputFile
    ind = metadata_dict['Case file'].index(case_name_list[i_case])
    metadata_fields = list(tuple(metadata_dict.keys())) ## cast dict.keys to tuple to list to make hashable

    for i_field in range(len(metadata_fields)):
      results_dict[metadata_fields[i_field]] = metadata_dict[metadata_fields[i_field]][ind]

    results_dict['L_ref'] = L_ref_m
    results_dict['S_ref'] = S_ref_m2
    results_dict['x_MRC_m'] = x_MRC_m
    results_dict['y_MRC_m'] = y_MRC_m
    results_dict['z_MRC_m'] = z_MRC_m

    #now append the various raw quantities you want
    keys = list(output_dict.keys())
    for i_key in range(len(keys)):
      if any(key_to_scale in keys[i_key] for key_to_scale in keys_to_scale) and not 'prime' in keys[i_key]:
        if any(key_to_correct in keys[i_key] for key_to_correct in keys_to_correct):
          key_index = next((i for i, s in enumerate(keys_to_correct) if s in keys[i_key]), -1)
          key_name = keys[i_key].replace(keys_to_correct[key_index],dim_keys[key_index].replace('_','b_')) #appends body subscript
          results_dict[key_name] = np.mean(output_dict[keys[i_key]][its_to_average_start:])/model_factor
        else:
          results_dict[keys[i_key]] = np.mean(output_dict[keys[i_key]][its_to_average_start:])/model_factor
      else:
        if 'prime' in keys[i_key]:
          key_name = keys[i_key].replace('my_prime','MYprime')
          results_dict[key_name] = np.mean(output_dict[keys[i_key]][its_to_average_start:])
        else:
          results_dict[keys[i_key]] = np.mean(output_dict[keys[i_key]][its_to_average_start:])

    #now sum forces over wall groups
    for i_group in range(len(wall_groups)):
      for i_quant in range(len(quants_to_sum)):
        quant = 0
        for i_key in range(len(keys)):
          if quants_to_sum[i_quant] in keys[i_key] and any(wall in keys[i_key] for wall in wall_groups[i_group]) and 'prime' not in keys[i_key]:
            if any(key_to_scale in keys[i_key] for key_to_scale in keys_to_scale) and not any(key_not_to_scale in wall_group_names[i_group] for key_not_to_scale in keys_not_to_scale):
              quant+=np.mean(output_dict[keys[i_key]][its_to_average_start:])/model_factor #average over last x its
            else:
              quant+=np.mean(output_dict[keys[i_key]][its_to_average_start:])
        results_dict[dim_keys[i_quant].replace('_','b_') + wall_group_names[i_group]] = quant
        if 'f' in quants_to_sum[i_quant]: #nondim
          results_dict[nondim_keys[i_quant].replace('_','b_') + wall_group_names[i_group]] = quant/(results_dict['q_0_BC']*results_dict['S_ref']) #force coeff
        elif 'm' in quants_to_sum[i_quant]:
          results_dict[nondim_keys[i_quant].replace('_','b_') + wall_group_names[i_group]] = quant/(results_dict['q_0_BC']*results_dict['S_ref']*results_dict['L_ref']) #mom_coeff

    results_keys = list(results_dict.keys())
    for i_key in range(len(results_keys)): #flip cad to body and zero quarter/halfspan lat dir
      if any(key_to_flip in results_keys[i_key] for key_to_flip in ['FXb','FZb','MXb','MZb','CXb','CZb','Clb','Cnb']) and not any(key_not_to_scale in results_keys[i_key] for key_not_to_scale in keys_not_to_scale) :
        results_dict[results_keys[i_key]] = results_dict[results_keys[i_key]]*-1 #convert to body x
      if any(key_to_zero in results_keys[i_key] for key_to_zero in ['FYb','MXb','MZb','CYb','Clb','Cnb']) and model_factor < 1 and not any(key_not_to_scale in results_keys[i_key] for key_not_to_scale in keys_not_to_scale): #zero lat dir for halfspan
        results_dict[results_keys[i_key]] = 0.0
      if any(key_to_zero in results_keys[i_key] for key_to_zero in ['MYb','Cmb']) and model_factor < 1/2 and not any(key_not_to_scale in results_keys[i_key] for key_not_to_scale in keys_not_to_scale): #zero longitudinal moment for quarter span
        results_dict[results_keys[i_key]] = 0.0
      if 'MYprime' in results_keys[i_key]:
        results_dict[results_keys[i_key].replace('MYprime','CH')] = results_dict[results_keys[i_key]]/(results_dict['q_0_BC']*results_dict['S_ref']*results_dict['L_ref'])

    #now making system dependent forces/converting to stability
    for i_group in range(len(wall_group_names)):
      forces_body = np.array([results_dict['FXb_' + wall_group_names[i_group]],results_dict['FYb_' + wall_group_names[i_group]],results_dict['FZb_' + wall_group_names[i_group]]])
      moments_body = np.array([results_dict['MXb_' + wall_group_names[i_group]],results_dict['MYb_' + wall_group_names[i_group]],results_dict['MZb_' + wall_group_names[i_group]]])

      if recordForces == True:
        results_dict['A_' + wall_group_names[i_group]] = -forces_body[0]
        results_dict['S_' + wall_group_names[i_group]] = forces_body[1]
        results_dict['N_' + wall_group_names[i_group]] = -forces_body[2]

      forces_stab = np.matmul(roty(results_dict['AoA_deg']),forces_body)
      moments_stab = np.matmul(roty(results_dict['AoA_deg']),moments_body)

      if recordForces == True:
        results_dict['FXs_' + wall_group_names[i_group]] = forces_stab[0]
        results_dict['FYs_' + wall_group_names[i_group]] = forces_stab[1]
        results_dict['FZs_' + wall_group_names[i_group]] = forces_stab[2]

        results_dict['Ds_' + wall_group_names[i_group]] = -forces_stab[0]
        results_dict['Ys_' + wall_group_names[i_group]] = forces_stab[1]
        results_dict['Ls_' + wall_group_names[i_group]] = -forces_stab[2]

        results_dict['MXs_' + wall_group_names[i_group]] = moments_stab[0]
        results_dict['MYs_' + wall_group_names[i_group]] = moments_stab[1]
        results_dict['MZs_' + wall_group_names[i_group]] = moments_stab[2]

      forces_body = np.array([results_dict['CXb_' + wall_group_names[i_group]],0,results_dict['CZb_' + wall_group_names[i_group]]])
      moments_body = np.array([0,results_dict['Cmb_' + wall_group_names[i_group]],0])

      results_dict['CA_' + wall_group_names[i_group]] = -forces_body[0]
      results_dict['CS_' + wall_group_names[i_group]] = forces_body[1]
      results_dict['CN_' + wall_group_names[i_group]] = -forces_body[2]

      forces_stab = np.matmul(roty(results_dict['AoA_deg']),forces_body)
      moments_stab = np.matmul(roty(results_dict['AoA_deg']),moments_body)

      results_dict['CXs_' + wall_group_names[i_group]] = forces_stab[0]
      results_dict['CYs_' + wall_group_names[i_group]] = forces_stab[1]
      results_dict['CZs_' + wall_group_names[i_group]] = forces_stab[2]

      results_dict['CDs_' + wall_group_names[i_group]] = -forces_stab[0]
      results_dict['CYs_' + wall_group_names[i_group]] = forces_stab[1]
      results_dict['CLs_' + wall_group_names[i_group]] = -forces_stab[2]

      results_dict['Cls_' + wall_group_names[i_group]] = moments_stab[0]
      results_dict['Cms_' + wall_group_names[i_group]] = moments_stab[1]
      results_dict['Cns_' + wall_group_names[i_group]] = moments_stab[2]

    results_dict['L/D_wall_aero'] = results_dict['CLs_wall_aero']/results_dict['CDs_wall_aero']
    results_dict['L/D_wall_prop'] = results_dict['CLs_wall_prop']/results_dict['CDs_wall_prop']
    results_dict['L/D_wall_total'] = results_dict['CLs_wall_total']/results_dict['CDs_wall_total']

    #now move everything over to the master dictionary with multiple case results
    results_keys = list(results_dict.keys())
    for i_key in range(len(results_keys)):
        if not any(key_to_ignore in results_keys[i_key] for key_to_ignore in keys_to_ignore):
            try:
                master_dict[results_keys[i_key]].append(results_dict[results_keys[i_key]])
            except:
                master_dict[results_keys[i_key]] = []
                master_dict[results_keys[i_key]].append(results_dict[results_keys[i_key]])

  ## get the name of the journal file
  for file in os.listdir(results_fileDir):
    if 'SWEEP' in file:
      journal_name = os.path.join(results_fileDir, file)
      break

  if not len(master_dict) == 0:
    results_df = pd.DataFrame(master_dict)
    update_master(results_df, file_path=journal_name.replace('.jou','')+'_results.xlsx')
    master_df = pd.read_excel(journal_name.replace('.jou','')+'_results.xlsx').sort_index(axis=1)

    master_df['CLs_wall_total^2'] = master_df['CLs_wall_total']**2
  else:
      print('!!ERROR!! The master dict is empty, so no result file will be created.')



  ## CREATE CANE CURVES, ONLY IF BACK PRESSSURED
  if isBackPressured == True:
    if not len(master_dict) == 0:
      ## create the plots directory
      plotDir = os.path.join(currDir, 'PLOTS')
      if not os.path.exists(plotDir):
        os.mkdir(plotDir)

      Create_Cane_Curves(plotDir, master_df, interface_surfs)


    '''
    # %% plot some stuff
    plotter = DataFramePlotter(master_df.sort_values(by=['Run','Point']), plotDir)
    #filter_conditions = {'Converged on NS':[True]}#, 'Run_ID': [1, 2]}
    #plotter.filter_data(filter_conditions)
    title = journal_name.replace('.jou','')
    plotter.plot_sweep(['Cms_wall_total','CLs_wall_total','CDs_wall_total','L/D_wall_total'], independent_vars='AoA_deg', metadata_columns=['mach_0_BC'], title=title)
    plotter.plot_sweep(['CDs_wall_total','CDs_wall_total'], independent_vars=['CLs_wall_total','CLs_wall_total^2'], metadata_columns=['mach_0_BC'], title=title)
    plotter.plot_sweep(['CN_wall_total','CA_wall_total','CLs_wall_total','CDs_wall_total'], independent_vars='AoA_deg', metadata_columns=['mach_0_BC'], title=title)


    ## Create Cane Curves here
    # parse the master_df
    # collect all vaues of MFR and PTR across the sweep
    # organize by specific sweep
    '''

def Create_Cane_Curves(plotDir, df, interface_surfs):

  print('\n>> Generating Plots...')

  INTERFACES = []
  for i_int in range(len(interface_surfs)):
    interface = interface_surfs[i_int].split('_')[-2]
    if 'station' in interface:
      interface = interface.replace('station', '')
    elif 'sta' in interface:
      interface = interface.replace('sta', '')
    elif 'STATION' in interface:
      interface = interface.replace('STATION', '')
    elif 'STA' in interface:
      interface = interface.replace('STA','')
    elif 'st' in interface:
      interface = interface.replace('st', '')
    elif 'ST' in interface:
      interface = interface.replace('ST', '')
    INTERFACES.append(interface)

  ## grab ALL freestream values
  p_0 = list(df['p_0_BC'])
  pt_0 = list(df['pt_0'])
  mach_0 = list(df['mach_0_BC'])
  rho_0 = list(df['rho_0'])
  vmag_0 = list(df['v_mag_0'])
  alt_0 = list(df['h_ft'])

  ## get all the unique mach values
  mach_vals = set()
  for i_mach in range(len(mach_0)):
    mach_vals.add(float(mach_0[i_mach]))

  mach_vals = sorted(list(mach_vals)) ## sort the list FIRST to ensure proper order

  ## Loop over every Mach value
  for i_mach in mach_vals:

    ## create the mach directories
    currPlotDir = os.path.join(plotDir, 'Mach_'+str(i_mach))
    if not os.path.exists(currPlotDir):
      os.mkdir(currPlotDir)

    ## Loop over every interface for each mach
    for i_int in range(len(INTERFACES)):

      ## create the interface directories for each mach value
      currIntPlotDir = os.path.join(currPlotDir, 'STATION_' + str(INTERFACES[i_int]).upper())
      if not os.path.exists(currIntPlotDir):
        os.mkdir(currIntPlotDir)

      ## grab the station quantities for each mach
      mdot_ST = list(df[df['mach_0_BC'] == i_mach]['mdot_' + str(INTERFACES[i_int])])
      backPress_ST = list(df[df['mach_0_BC'] == i_mach]['BackPress_Pa'])
      pt_ST = list(df[df['mach_0_BC'] == i_mach]['pt_' + str(INTERFACES[i_int])])
      p_ST = list(df[df['mach_0_BC'] == i_mach]['p_' + str(INTERFACES[i_int])])

      ## grab the freestream quantities for each mach
      p_0 = list(df[df['mach_0_BC'] == i_mach]['p_0_BC'])
      pt_0 = list(df[df['mach_0_BC'] == i_mach]['pt_0'])
      mach_0 = list(df[df['mach_0_BC'] == i_mach]['mach_0_BC'])
      rho_0 = list(df[df['mach_0_BC'] == i_mach]['rho_0'])
      vmag_0 = list(df[df['mach_0_BC'] == i_mach]['v_mag_0'])
      alt_0 = list(df[df['mach_0_BC'] == i_mach]['h_ft'])

      ## Sort all the values in ascending order of backPressure
      sorted_vars = sorted(zip(backPress_ST, mdot_ST, pt_ST, p_ST, p_0, pt_0, mach_0, rho_0, vmag_0, alt_0))

      ## unzip back into seperate lists
      backPress_ST, mdot_ST, pt_ST, p_ST, p_0, pt_0, mach_0, rho_0, vmag_0, alt_0 = zip(*sorted_vars)

      ## convert all to lists
      backPress_ST = list(backPress_ST)
      mdot_ST = list(mdot_ST)
      pt_ST = list(pt_ST)
      p_ST = list(p_ST)
      p_0 = list(p_0)
      pt_0 = list(pt_0)
      mach_0 = list(mach_0)
      rho_0 = list(rho_0)
      vmag_0 = list(vmag_0)
      alt_0 = list(alt_0)

      currMach = mach_0[0]
      currAlt = alt_0[0]
      currInterface = INTERFACES[i_int]
      unstart_vars=Calculate_critical_value(mdot_ST, p_ST)
      maxBP = backPress_ST[unstart_vars[1]]
      caseVars = [currMach, currAlt, currInterface, maxBP]

      ## plot the different cane curves at each station
      Plot_Cane_Curves(currIntPlotDir, x_var = mdot_ST, y_var = p_ST, x_name = 'Mdot', y_name = 'static_pressure', case_vars = caseVars, unstart_vars=Calculate_critical_value(mdot_ST, p_ST))
      Plot_Cane_Curves(currIntPlotDir, x_var = mdot_ST, y_var = pt_ST, x_name = 'Mdot', y_name = 'total_pressure', case_vars = caseVars, unstart_vars=Calculate_critical_value(mdot_ST, pt_ST))

      # calculate MFR
      MFR = []
      PTR = []
      for i in range(len(rho_0)):
        MFR.append(mdot_ST[i]/ (rho_0[0] * vmag_0[0] * 0.00553))
        PTR.append(pt_ST[i] / pt_0[i])

      Plot_Cane_Curves(currIntPlotDir, x_var = MFR, y_var = PTR, x_name = 'MFR', y_name = 'PTR', case_vars = caseVars, unstart_vars=Calculate_critical_value(MFR, PTR))

def Calculate_critical_value(x_var, y_var):
  critical = x_var[0] * 0.98
  max_idx = 0
  for i in range(len(x_var)):
    if i == 0:
      if x_var[i] >= critical:
        max_idx = i
    else:
      if x_var[i] >= critical and y_var[i] >= y_var[i-1]:
        max_idx = i
      elif x_var[i] >= critical and y_var[i] <= y_var[i-1]:
        max_idx = i

  return([critical, max_idx])

def Plot_Cane_Curves(plotDir, x_var = [], y_var = [], x_name = '', y_name = '', case_vars = None, unstart_vars=None):

  mach = case_vars[0]
  alt = case_vars[1]
  interface = case_vars[2]
  maxBP = case_vars[3]

  critical = unstart_vars[0]
  max_idx = unstart_vars[1]


  plt.grid(True, which = 'major')
  plt.plot([critical, critical], [min(y_var), max(y_var)], 'k--', label = 'Unstart limit (2% deviation)')
  plt.plot(x_var, y_var, 'x-')
  plt.plot(x_var[max_idx], y_var[max_idx], 'rx', label = 'Critical Point. (Back Pressure: ' + str(maxBP) + ' Pa)')
  title = 'Mach ' +str(mach) + ' Alt ' +str(alt) + ' ft\nStation ' + str(interface) + ': ' + str(x_name) +' vs. ' +str(y_name)
  plt.title(title)
  plt.xlabel(x_name)
  plt.ylabel(y_name)


  plt.legend()
  plt.savefig(os.path.join(plotDir, 'Station_' + str(interface) + '_' +str(x_name) + '_vs_' + str(y_name) +'.png'))
  plt.close()


def Create_Interface_Outputs(outputDir, base_case, case_name_list, metadata_dict, interface_surfs):
  global x_MRC_m, y_MRC_m, z_MRC_m
  global L_ref_m, S_ref_m2
  global model_factor

  master_dict = {}

  its_to_average_start = -1 ## all cases converged on NS residuals

  for i_case in range(len(case_name_list)):
    outputFile = os.path.join(rescale_fileDir, case_name_list[i_case] + '_reports1.out' )
    if os.path.exists(outputFile):
      try:
        output_dict = read_fluent_out_to_dict(outputFile)
      except Exception as error:
        print('!! ERROR !! While parsing the file ' + str(outputFile) + ' the following error occured:\n')
        print(error)
        continue
    else:
      print('\n>> Warning... The file: ' +str(outputFile) + ' could not be found. Skipping this entry.\n')
      continue

    #first append the metadata you wrote to previously
    results_dict = {}
    results_dict['Base case file'] = base_case
    results_dict['Output file'] = outputFile
    ind = metadata_dict['Case file'].index(case_name_list[i_case])
    metadata_fields = list(tuple(metadata_dict.keys())) ## cast dict.keys to tuple to list to make hashable

    for i_field in range(len(metadata_fields)):
      results_dict[metadata_fields[i_field]] = metadata_dict[metadata_fields[i_field]][ind]

    results_dict['L_ref'] = L_ref_m
    results_dict['S_ref'] = S_ref_m2
    results_dict['x_MRC_m'] = x_MRC_m
    results_dict['y_MRC_m'] = y_MRC_m
    results_dict['z_MRC_m'] = z_MRC_m

    #now append the various raw quantities you want
    keys = list(output_dict.keys())
    for i_key in range(len(keys)):
      results_dict[keys[i_key]] = np.mean(output_dict[keys[i_key]][its_to_average_start:])

    #now move everything over to the master dictionary with multiple case results
    results_keys = list(results_dict.keys())
    for i_key in range(len(results_keys)):
      try:
        master_dict[results_keys[i_key]].append(results_dict[results_keys[i_key]])
      except:
        master_dict[results_keys[i_key]] = []
        master_dict[results_keys[i_key]].append(results_dict[results_keys[i_key]])

  ## get the name of the journal file
  for file in os.listdir(rescale_fileDir):
    _,ext = os.path.splitext(file)
    if ext.lower() == '.jou' and not '.original' in file:
      journal_name = os.path.join(rescale_fileDir, file)
      break

  if not len(master_dict) == 0:
    results_df = pd.DataFrame(master_dict)
    update_master(results_df, file_path=journal_name.replace('.jou','')+'_results_integrations.xlsx')
    master_df = pd.read_excel(journal_name.replace('.jou','')+'_results_integrations.xlsx').sort_index(axis=1)
  else:
    print('!!ERROR!! The master dict is empty, so no result file will be created.')


def move_report_files(currDir, case_list):

  ## move the report files into the REPORTS directory
  reportDir = os.path.join(currDir, 'REPORT_FILES')
  if not os.path.exists(reportDir):
    print('!! ERROR !! The REPORT_FILES directory coes not exist. 2D_journal_setup.py must be run before 2D_post. Exiting now.')
    sys.exit(0)

  report_files = []
  for file in os.listdir(currDir):
    if 'reports' in file:
      currReportFile = os.path.join(currDir, file)
      report_files.append([currReportFile, file])

  filesAlreadyMoved = False
  if len(report_files) == 0:
    report_files = []
    for file in os.listdir(reportDir):
      report_files.append(file)
    if not len(report_files) == 0:
      filesAlreadyMoved = True
    else:
      print('!! ERROR !! No report files were found in the current directory or the REPORT FILES directory. Exiting now.')
      sys.exit(0)


  ## move the files
  if filesAlreadyMoved == False:
    print('\n>> moving report files into the REPORT_FILES directory...')
    for i_file in range(len(report_files)):
      orig = report_files[i_file][0]
      new = os.path.join(reportDir, report_files[i_file][1])
      shutil.move(orig, new)

  ## check the report files against the case list
  for i_case in case_list:
    caseFound = False
    for i_file in range(len(report_files)):
      if i_case in report_files[i_file][1]:
        caseFound = True
    #if caseFound == False:
      #print('~~ WARNING ~~ The report file for case ' + i_case + ' could NOT be found. Skipping this entry...')

  ## copy the journal file into the reports dir
  journalFile = ''
  journalFound = False
  for file in os.listdir(currDir):
    if os.path.isfile(file):
      if 'SWEEP' in file:
        orig = os.path.join(currDir, file)
        new = os.path.join(reportDir, file)
        shutil.copy(orig, new)
        journalFound = True
        break

  if journalFound == False:
    print('!! ERROR !! The sweep journal was not found in ' + str(currDir) + '. Exiting now.')
    sys.exit(0)

  return(reportDir)


def CleanDir(dir):
  ## delete all files and directories from a directory
  dir = dir.replace('\\\\', '\\').replace('\"', '')
  for file in os.listdir(dir):
    temp = os.path.join(dir, file)
    if os.path.isfile(temp):
      os.remove(temp)
    elif os.path.isdir(temp):
      CleanDir(temp)
      try:
        os.rmdir(temp)
      except FileNotFoundError:
        continue

  os.rmdir(dir)


def CopyDownloadedFiles(tempDir, origDir):
  ## get all the newly downloaded files in the temp_fileDir
  tempDir = tempDir.replace('\\\\', '\\').replace('\"', '')
  origDir = origDir.replace('\\\\', '\\').replace('\"', '')

  newFiles = []
  for file in os.listdir(tempDir):
    filePath = os.path.join(tempDir, file)
    if os.path.isfile(filePath):
      newFiles.append([file, filePath])

  ## Copy ALL the newly downloaded files to the original folder
  for i_file in range(len(newFiles)):
      origFile = newFiles[i_file][1]
      copyFile = os.path.join(origDir, newFiles[i_file][0])
      shutil.copy(origFile, copyFile)

def clean_up_image_dir(currDir, case_list):

  print('\n>> Organizing images...')

  imageDir = os.path.join(currDir, 'IMAGES')
  if not os.path.exists(imageDir):
    print('!! ERROR !! The IMAGES directory does not exist. Please run 2D_journal_setup.py before 2D_post.py. Exiting now.')
    sys.exit(0)

  ## remove the case name for easier readability and shorter folder names
  caseID = case_list[0].split('_')
  for iVar in range(len(caseID)):
    if caseID[iVar] == 'Mach':
      break
  caseRemove = ''
  for i in range(iVar):
    caseRemove += caseID[i] + '_'

  ## get the contour names
  contours = set()
  origContours = []
  newImages = []
  ImageSort = False
  for i_dir in os.listdir(imageDir):
    temp = os.path.join(imageDir, i_dir)
    if os.path.isdir(temp):
      if not str(i_dir) in case_list:
        if not '_COMPARE' in temp:
          if not 'Mach_' in str(i_dir) and not 'h_' in str(i_dir):
            compareDir = os.path.join(imageDir, temp + '_COMPARE')
            if os.path.exists(compareDir): ## if the corresponding COMPARE dir already exists, move the files into the existing directory.
              currImageDir = os.path.join(imageDir, temp)
              contours.add(str(i_dir)+'_COMPARE')
              ## remove the case identifier from the image name
              for pic in os.listdir(currImageDir):
                orig = os.path.join(currImageDir, pic)
                new = os.path.join(compareDir, str(pic).replace(caseRemove,''))
                newImages.append(new)
                shutil.move(orig, new)
                ImageSort = True
            ## if the COMPARE directory does not exist, create it
            else: ## rename the directories
              orig = os.path.join(imageDir, i_dir)
              new = os.path.join(imageDir, str(i_dir) + '_COMPARE')
              os.rename(orig, new)
              contours.add(str(i_dir)+'_COMPARE')
              origContours.append(str(i_dir))
              ## remove the case identifier from the image name
              for pic in os.listdir(new):
                origPic = os.path.join(new, pic)
                newPic = os.path.join(new, str(pic).replace(caseRemove, ''))
                os.rename(origPic, newPic)
                ImageSort = True

  contours = list(contours)


  if ImageSort == True:
    ## create the individual case directories
    caseDirs = []
    for i_case in case_list:
      i_case = i_case.replace(caseRemove, '')
      caseDir = os.path.join(imageDir, str(i_case))
      if not os.path.exists(caseDir):
        os.mkdir(caseDir)
      caseDirs.append([i_case, caseDir])

    for i_contour in contours:
      contourDir = os.path.join(imageDir, i_contour)
      ## if there are no new images, it is the first time running 2D_post
      if len(newImages) == 0:
        for i_pic in os.listdir(contourDir):
          for i_case in range(len(caseDirs)):
            if str(caseDirs[i_case][0]) in str(i_pic):
              orig = os.path.join(contourDir, i_pic)
              #new = os.path.join(caseDirs[i_case][1], str(i_pic).split('_')[-1] ) ## in the individual directories, name the images just the contour type (i.e. Mach, static_pressure, etc)
              new = os.path.join(caseDirs[i_case][1], str(i_pic).replace(caseDirs[i_case][0],'').lstrip('_'))
              shutil.copy(orig, new)
      else:
        for i_pic in range(len(newImages)):
          for i_case in range(len(caseDirs)):
            i_pic_name = newImages[i_pic].split('\\')[-1]
            if str(caseDirs[i_case][0]) in i_pic_name:
              orig = newImages[i_pic]
              new = os.path.join(caseDirs[i_case][1], i_pic_name.replace(caseDirs[i_case][0],'').lstrip('_')) ## in the individual directories, name the images just the contour type (i.e. Mach, static_pressure, etc)
              shutil.copy(orig, new)


    ## arrange the COMPARE directories based on Mach and altitude combinations
    mach_alt = set() ## initialize as a set to prevent duplicates
    if len(newImages) == 0: ## If no new images, first time running and no directories will exist
      contourDir_conditions = os.path.join(imageDir, contours[0])
      for i_pic in os.listdir(contourDir_conditions):
        case = i_pic.split('_')
        mach = ''
        alt = ''
        for i_var in range(len(case)):
          if case[i_var] == 'Mach':
            mach = case[i_var+1]
          elif case[i_var] == 'h':
            alt = case[i_var+1]
            break
        if not mach == '' and not alt == '':
          mach_alt.add(tuple([mach, alt]))

    else: ## If not first time running, existing mach_alt directories will be present
      contourDir_conditions = os.path.join(imageDir, contours[0])
      ## check for existing mach_alt directories
      for dir in os.listdir(contourDir_conditions):
        tempDir = os.path.join(contourDir_conditions, dir)
        if os.path.isdir(tempDir):
          case  = dir.split('_')
          mach = ''
          alt = ''
          for i_var in range(len(case)):
            if case[i_var] == 'Mach':
              mach = case[i_var+1]
            elif case[i_var] == 'h':
              alt = case[i_var+1]
              break
          if not mach == '' and not alt == '':
            mach_alt.add(tuple([mach, alt]))

      ## check for new images
      for file in os.listdir(contourDir_conditions):
        currFile = os.path.join(contourDir_conditions, file)
        if os.path.isfile(currFile):
          case = file.split('_')
          mach = ''
          alt = ''
          for i_var in range(len(case)):
            if case[i_var] == 'Mach':
              mach = case[i_var+1]
            elif case[i_var] == 'h':
              alt = case[i_var+1]
              break
          if not mach == '' and not alt == '':
            mach_alt.add(tuple([mach, alt]))

    ## cast the set to a list to make it hashable
    mach_alt = [list(i_comb) for i_comb in mach_alt]

    ## make the mach and alt combination directories and move the pictures
    for i_contour in contours:
      contourDir = os.path.join(imageDir, i_contour)

      ## check if the mach_alt directories exist, if not create them
      MA_dirs = []
      for i_MA in range(len(mach_alt)):
        mach_alt_name = 'Mach_' + str(mach_alt[i_MA][0]) + '_h_' + str(mach_alt[i_MA][1])
        mach_alt_dir = os.path.join(contourDir, mach_alt_name)
        if not os.path.exists(mach_alt_dir):
          os.mkdir(mach_alt_dir)
        MA_dirs.append([mach_alt_name, mach_alt_dir])

      ## loop over all the images in the current directory and move them
      allPics = []
      for pic in os.listdir(contourDir):
        currPic = os.path.join(contourDir, pic)
        if os.path.isfile(currPic): ## only grab images, dont grab directories
          allPics.append(pic)

      for i_pic in allPics:
        for i_MA in range(len(MA_dirs)):
          if str(MA_dirs[i_MA][0]) in i_pic:
            orig = os.path.join(contourDir, i_pic)
            new = os.path.join(MA_dirs[i_MA][1], i_pic)
            shutil.move(orig, new)
            break

    ## remove any empty directories
    for i_dir in os.listdir(imageDir):
      temp = os.path.join(imageDir, i_dir)
      if os.path.isdir(temp):
        if not os.listdir(temp):
          os.rmdir(temp)
        if str(i_dir) in origContours:
          os.rmdir(i_dir)

if __name__ == '__main__':

  ## get the working directory
  currDir = os.getcwd()

  ## read command line arguments
  case, case_path, recordForces = Read_Inputs(currDir)

  ## read the data from the journal_setup.py files
  base_case, journal_name, case_list, metadata_dict, surface_list, interface_surfs, outputDir = Read_Journal_Outputs(currDir)

  ## check that the files exist
  results_fileDir = move_report_files(currDir, case_list)

  ## collect and output data from the wall surfaces
  Create_Wall_Outputs(currDir, results_fileDir, base_case, case_list, metadata_dict, surface_list, interface_surfs, recordForces)

  ## clean up image directory
  #clean_up_image_dir(currDir, case_list)

  ## collect and output data from the interface surfaces
  #Create_Interface_Outputs(outputDir, base_case, case_list, metadata_dict, interface_surfs)

  ## Clean up function to remove un-necessary files
  ## what files can be deleted?
