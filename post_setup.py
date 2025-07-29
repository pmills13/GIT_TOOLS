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
from plotterFunctions import *
from viscous_inviscid_forces import Calculate_Viscous_Inviscid_Forces

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
from scipy.interpolate import interp1d
import json

## setup the argparser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("-forces", help="Optional flag that will report raw Force values in addition to Coefficients. (-forces  (no other input needed)) If not entered, defaults to Coefficients only.\nUsage: -forces", action="store_true", required=False)
parser.add_argument("-model_factor", help="Model factor for the case. \nIf it is a halfspan model, enter 0.5.\nIf it is a fullspan model, enter 1.0.\nIf It is a quarterspan model, enter 0.25.\nIf not entered, defaults to halfspan.\nUsage: -model_factor 0.5", type=float, required = False)
parser.add_argument("-jobID", help="Rescale job ID for the files you want to pull down.\nNOTE: If you do not enter a jobID, the code will search for the jobId from the job_submission_output file created by the journal_setup.py. If you dont enter a jobID AND the job_submission_output file cant be found, the code will error out. .\nUsage: -jobID #####", type=str, required = False)

global MRP_m
global L_ref_m, S_ref_m2, model_factor, isBackPressured
global jobID
global recordForces

###############################
## READ_INPUTS
###############################
def Read_Inputs(currDir):
  global jobID
  global model_factor
  global recordForces

  args = parser.parse_args()
  model_factor = args.model_factor
  recordForces = args.forces
  ID = args.jobID
  if ID == None:
    jobID = []
  else:
    jobID = [ID]

  if model_factor == None:
    model_factor = 0.5

  ## Set the working directory
  working_dir = ''
  rescaleDir = ''
  dirFound = False
  if 'CFD' in currDir:
    working_dir = currDir
    dirFound = True

  else:
    for dir in os.listdir(currDir):
      tempDir = os.path.join(currDir, dir)
      if os.path.isdir(tempDir):
        if dir == 'CFD':
          working_dir = tempDir
          dirFound = True
          break

  if dirFound == False:
    print('!! ERROR !! The \'CFD\' directory could not be found. Please run journal_setup.py BEFORE running post_setup.py. Exiting now')
    sys.exit(0)

  return(working_dir)

###############################
## READ_JOURNAL_OUTPUTS
###############################

def Read_Journal_Outputs(currDir):
  global MRP_m
  global L_ref_m, S_ref_m2, model_factor, isBackPressured
  global jobID

  ## check that the OUTPUTS directory was created and the 3 required files are in it
  metaFileFound = False
  caseFileFound = False
  surfFileFound = False
  jobFileFound = False
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
            elif 'job_submission_output' in file:
              jobFileFound = True
              postFile_Job = os.path.join(tempDir, file)
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
      case 'MRP_m:':
        MRP_m = [float(caseData[entry][1]), float(caseData[entry][2]), float(caseData[entry][3])]
      case 'L_ref:':
        L_ref_m = float(caseData[entry][1])
      case 'S_ref:':
        S_ref_m2 = float(caseData[entry][1])
      case 'model_factor:':
        model_factor = float(caseData[entry][1])
      case 'isBackPressured:':
        BP = caseData[entry][1].lower()
        if BP == 'true':
          isBackPressured = True
        elif BP == 'false':
          isBackPressured = False
        else:
          print('!! ERROR !! Non-boolean value for the variable \'isBackPressured\'. Exiting now.')
          sys.exit(0)
      case _:
        print('>> Warning... Unexpected header (' + str(header) + ') in ' +str(postFile_Case) +'. Skipping this entry.\n')

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

  ## Read the job_submission_output file to get the jobID
  ids = []
  if jobFileFound == True:
    jobData = []
    with open(postFile_Job, 'r') as fin:
      for line in fin:
        jobData.append(line.strip().split())
    fin.close()

    for i_line in range(len(jobData)):
      try:
        if jobData[i_line][3] == 'Job':
          ids.append(jobData[i_line][4].replace(':',''))
      except IndexError:
        doNothing = True

    if len(ids) == 0:
      if len(jobID) == 0:
        print('\n!! ERROR !! The -jobID flag was not passed, and the job_submission_script was found, but no viable jobID was detected in the job_submission_file. Exiting now.\n')
        sys.exit(0)
    else:
      for i_id in range(len(ids)):
        if not ids[i_id] == ids[0]:
          if len(jobID) == 0:
            print('\n!! ERROR !! The -jobID flag was not passed, and the job_submission_script was found, but contains multiple jobIDS. Exiting now\n')
            sys.exit(0)

    if len(jobID) == 0:
      jobID.append(ids[0])
      print('\n>> The jobID was found: ' + str(jobID[0]) + ' ...\n')
    else:
      if not jobID[0] == ids[0]:
        print('\n\n\n' + '*'*100)
        response = input('WARNING... The jobID you entered (' + str(jobID[0]) + ') does NOT match the jobID in the job_submission_output file (' + str(ids[0]) + '). Are you positive you want to continue?\nType y for yes or n for no and press Enter to continue.\n')
        if response.lower() == 'n':
          print('\nStopping the program now...\n')
          sys.exit(0)
        elif response.lower() == 'y':
          print('\n>>>> Continuing post-processing using jobID: ' + str(jobID[0]) + ' ... \n')
        else:
          print('\nInvalid entry of ' + str(response) + '. Continuing post-processing using jobID: ' +str(ids[0]) + ' ...\n')
          jobID = [ids[0]]
  else:
    if len(jobID) == 0:
      print('\n!! ERROR !! The job_submission_output file was not found and the -jobID flag was not passed. The code will not be able to pull anything down from Rescale. Exiting now.\n')
      sys.exit(0)

  return(base_case, journal_name, case_name_list, metadata_dict, surface_list, interface_list, outputDir)

###############################
## CREATE_WALL_OUTPUTS
###############################
def Create_Wall_Outputs(outputDir,rescale_fileDir, rescale_imageDir, base_case, case_name_list, metadata_dict, surface_list):

  global MRP_m
  global L_ref_m, S_ref_m2
  global model_factor
  global isBackPressured
  global recordForces

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
  keys_not_to_scale = ['canard','fin', 'strake'] #outputs containing these strings won't be scaled (use for individual canards, fins, etc.)
  internal_keys_to_scale = ['a_', 'mdot_'] ##all outputs containing these strings will be divided by model factor to arrive at fullspan values

  keys_to_correct = keys_to_scale
  quants_to_sum = keys_to_scale

  dim_keys = ['FX_','FY_','FZ_','MX_','MY_','MZ_'] #restoring correct cases
  nondim_keys = ['CX_','CY_','CZ_','Cl_','Cm_','Cn_']

  for i_case in range(len(case_name_list)):
    outputFile = os.path.join(rescale_fileDir, case_name_list[i_case] + '_reports.out' )
    if os.path.exists(outputFile):
      try:
        output_dict = read_fluent_out_to_dict(outputFile)
      except Exception as error:
        print('\n!! ERROR !! While parsing the file ' + str(outputFile) + ' the following error occured:\n')
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
    metadata_fields = list(tuple(metadata_dict.keys())) ## cast dict.keys-> tuple -> list to make hashable

    for i_field in range(len(metadata_fields)):
      results_dict[metadata_fields[i_field]] = metadata_dict[metadata_fields[i_field]][ind]

    results_dict['L_ref'] = L_ref_m
    results_dict['S_ref'] = S_ref_m2
    results_dict['x_MRC_m'] = MRP_m[0]
    results_dict['y_MRC_m'] = MRP_m[1]
    results_dict['z_MRC_m'] = MRP_m[2]

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

      elif any(internal_key_to_scale in keys[i_key] for internal_key_to_scale in internal_keys_to_scale):
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
  for file in os.listdir(rescale_fileDir):
    _,ext = os.path.splitext(file)
    if ext.lower() == '.jou' and not '.original' in file:
      journal_name = os.path.join(rescale_fileDir, file)
      break

  if not len(master_dict) == 0:
    results_df = pd.DataFrame(master_dict)
    update_master(results_df, file_path=journal_name.replace('.jou','')+'_results.xlsx')
    master_df = pd.read_excel(journal_name.replace('.jou','')+'_results.xlsx').sort_index(axis=1)

    master_df['CLs_wall_total^2'] = master_df['CLs_wall_total']**2

    ## create the plots directory
    plotDir = os.path.join(rescale_fileDir, 'plots')
    if not os.path.exists(plotDir):
      os.mkdir(plotDir)

    # %% plot some stuff
    #plotter = DataFramePlotter(master_df.sort_values(by=['Run','Point']), plotDir)
    #filter_conditions = {'Converged on NS':[True]}#, 'Run_ID': [1, 2]}
    #plotter.filter_data(filter_conditions)
    #title = journal_name.replace('.jou','')
    #plotter.plot_sweep(['Cms_wall_total','CLs_wall_total','CDs_wall_total','L/D_wall_total'], independent_vars='AoA_deg', metadata_columns=['mach_0_BC'], title=title)
    #plotter.plot_sweep(['CDs_wall_total','CDs_wall_total'], independent_vars=['CLs_wall_total','CLs_wall_total^2'], metadata_columns=['mach_0_BC'], title=title)
    #plotter.plot_sweep(['CN_wall_total','CA_wall_total','CLs_wall_total','CDs_wall_total'], independent_vars='AoA_deg', metadata_columns=['mach_0_BC'], title=title)


    ## Create Cane curves, only if back pressured
    if isBackPressured == True:
      ## create the plots directory
      if not os.path.exists(plotDir):
        os.mkdir(plotDir)
      Create_Cane_Curves(plotDir, master_df, interface_surfs)

  else:
    print('!!ERROR!! The master dict is empty, so no result file will be created.')


def Create_Interface_Outputs(outputDir, base_case, case_name_list, metadata_dict, interface_surfs):
  global MRP_m
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
    results_dict['x_MRC_m'] = MRP_m[0]
    results_dict['y_MRC_m'] = MRP_m[1]
    results_dict['z_MRC_m'] = MRP_m[2]

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

def Pull_Down_Rescale_Files(outputDir, journal_name):
  global jobID

  rescale_job_id = jobID
  job_name = [journal_name.replace('.jou','')]
  filetypes = ['*.out','*.jou','*.log','*.png', '.*txt']

  ## create the bat file to pull down files from rescale
  fName = 'pullDownRescale.bat'
  pullFile = os.path.join(outputDir, fName)
  if os.path.exists(pullFile):
    os.remove(pullFile)

  useTempDir = False
  with open(pullFile, 'w') as fin:
    for i_job in range(len(rescale_job_id)): ## Currently only supports a single job ID
      for i_file in range(len(filetypes)):
        if i_file == 0:

          ## maunally create directories instead of having them created via .bat file (cleaner)
          rescale_fileDir = os.path.join(outputDir, job_name[i_job] + '_outputs')
          if not os.path.exists(rescale_fileDir):
            os.mkdir(rescale_fileDir)
          else:
            temp_fileDir = os.path.join(outputDir, job_name[i_job] + '_outputs_TEMP') ## Rescale has weird behavior when downloading files. If the downlaod directory already exists, it doesnt download the files.
            useTempDir = True
            if os.path.exists(temp_fileDir):
              CleanDir(temp_fileDir)
            os.mkdir(temp_fileDir)

          rescale_imageDir = os.path.join(rescale_fileDir, 'images')
          if useTempDir == False:
            if not os.path.exists(rescale_imageDir):
              os.mkdir(rescale_imageDir)
            else:
              temp_imageDir = os.path.join(temp_fileDir, 'images_TEMP')
              if os.path.exists(temp_imageDir):
                CleanDir(temp_imageDir)
              os.mkdir(temp_imageDir)
          else:
            temp_imageDir = os.path.join(temp_fileDir, 'images_TEMP')
            if os.path.exists(temp_imageDir):
              CleanDir(temp_imageDir)
            os.mkdir(temp_imageDir)

          rescale_forceDir = os.path.join(rescale_fileDir, 'forces')
          if useTempDir == False:
            if not os.path.exists(rescale_forceDir):
              os.mkdir(rescale_forceDir)
            else:
              temp_forceDir = os.path.join(temp_fileDir, 'forces_TEMP')
              if os.path.exists(temp_forceDir):
                CleanDir(temp_forceDir)
              os.mkdir(temp_forceDir)
          else:
            temp_forceDir = os.path.join(temp_fileDir, 'forces_TEMP')
            if os.path.exists(temp_forceDir):
              CleanDir(temp_forceDir)
            os.mkdir(temp_forceDir)


          ## rescale requires double backslashes in file paths enclosed by " "
          if useTempDir == False:
            rescale_fileDir = '\"' + rescale_fileDir.replace('\\', '\\\\') +'\"'
            rescale_imageDir = '\"' + rescale_imageDir.replace('\\', '\\\\') + '\"'
            rescale_forceDir = '\"' + rescale_forceDir.replace('\\', '\\\\') + '\"'
          else:
            temp_fileDir = '\"' + temp_fileDir.replace('\\', '\\\\') + '\"'
            temp_imageDir = '\"' + temp_imageDir.replace('\\', '\\\\') + '\"'
            temp_forceDir = '\"' + temp_forceDir.replace('\\', '\\\\') + '\"'

        if 'png' in filetypes[i_file]:
          if useTempDir == False:
            fin.write('rescale-cli -X https://itar.rescale.com sync -j ' + rescale_job_id[i_job] + ' -o ' + rescale_imageDir + ' -f ' + filetypes[i_file] + '\n')
          else:
            fin.write('rescale-cli -X https://itar.rescale.com sync -j ' + rescale_job_id[i_job] + ' -o ' + temp_imageDir+ ' -f ' + filetypes[i_file] + '\n')
        elif 'txt' in filetypes[i_file]:
          if useTempDir == False:
            fin.write('rescale-cli -X https://itar.rescale.com sync -j ' + rescale_job_id[i_job] + ' -o ' + rescale_forceDir + ' -f ' + filetypes[i_file] + '\n')
          else:
            fin.write('rescale-cli -X https://itar.rescale.com sync -j ' + rescale_job_id[i_job] + ' -o ' + temp_forceDir+ ' -f ' + filetypes[i_file] + '\n')
        else:
          if useTempDir == False:
            fin.write('rescale-cli -X https://itar.rescale.com sync -j ' + rescale_job_id[i_job] + ' -o ' + rescale_fileDir + ' -f ' + filetypes[i_file] + '\n')
          else:
            fin.write('rescale-cli -X https://itar.rescale.com sync -j ' + rescale_job_id[i_job] + ' -o ' + temp_fileDir + ' -f ' + filetypes[i_file] + '\n')

    fin.write('\n')
  fin.close()

  ## execute the bat file
  print('>> Attempting to pull down files from rescale...')
  try:
    result = subprocess.run([pullFile], shell=True, check=True, capture_output=True, text=True)
  except subprocess.CalledProcessError as error:
    print('!!ERROR!! Tried to pull down files from rescale resulted in the following error:\n' + error.stderr)
    print('Exiting now...')
    sys.exit(0)

  if result.stderr == '':
    print('>>>> files pulled down from rescale successfully.')

  ## If files were downloaded to the temp directory, copy them to the proper directory.
  if useTempDir == True:
    CopyDownloadedFiles(temp_fileDir, rescale_fileDir)
    CopyDownloadedFiles(temp_imageDir, rescale_imageDir)
    CopyDownloadedFiles(temp_forceDir, rescale_forceDir)
    ## After copying the files, delete the temp_dir
    CleanDir(temp_fileDir)

  ## before returning the new directories, remove the enclosing "" and switch back to single slash for code compatability
  rescale_fileDir = rescale_fileDir.replace('\"','').replace('\\\\','\\')
  rescale_imageDir = rescale_imageDir.replace('\"','').replace('\\\\','\\')
  rescale_forceDir= rescale_forceDir.replace('\"','').replace('\\\\', '\\')

  return(rescale_fileDir, rescale_imageDir, rescale_forceDir)



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

def Clean_Up_Directory(rescale_fileDir):
  files_to_remove = ['aclcleanup', 'ansyscl', 'licdebug']

  for file in os.listdir(rescale_fileDir):
    currFile = os.path.join(rescale_fileDir, file)
    if os.path.isfile(currFile):
      for bad_file_name in files_to_remove:
        if bad_file_name in file:
          os.remove(currFile)
          break


if __name__ == '__main__':

  ## get the working directory
  currDir = os.getcwd()

  ## read command line arguments
  working_dir = Read_Inputs(currDir)

  ## set the currDir to the working dir
  currDir = working_dir

  ## read the data from the journal_setup.py files
  base_case, journal_name, case_list, metadata_dict, surface_list, interface_surfs, outputDir = Read_Journal_Outputs(currDir)

  ## pull down files from Rescale
  rescale_fileDir, rescale_imageDir, rescale_forceDir = Pull_Down_Rescale_Files(outputDir, journal_name)

  ## create the residual plots
  Create_Residual_Plots(rescale_fileDir)

  ## create the viscous/inviscid files
  Calculate_Viscous_Inviscid_Forces(rescale_forceDir, case_list)

  ## collect and output data from the wall surfaces
  Create_Wall_Outputs(outputDir, rescale_fileDir, rescale_imageDir, base_case, case_list, metadata_dict, surface_list)

  ## collect and output data from the interface surfaces
  Create_Interface_Outputs(outputDir, base_case, case_list, metadata_dict, interface_surfs)

  ## Clean up function to remove un-necessary files
  Clean_Up_Directory(rescale_fileDir)

