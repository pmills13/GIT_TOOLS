# -*- coding: utf-8 -*-

## CODE DEPENDANCIES: IF running 2D, you need pyFluent.    ---> pip install ansys-fluent-core

"""
@author: pmills
"""

import sys
import os
import shutil
import subprocess
import argparse

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
from miscellaneous_subroutines import Calculate_Total_Temp
from helper_functions_2D import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
from scipy.interpolate import interp1d
import json

## setup the argparser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("-case", help="Identifying string for the mesh file. If no case is entered, The code will select the most recently created mesh workflow file directory in the current directory.\nNOTE: you do not need to include \'workflow_files\' in the case name.", type=str, required=False)
parser.add_argument("-model_factor", help="Model factor for the case. \nIf it is a halfspan model, enter 0.5.\nIf it is a fullspan model, enter 1.0.\nIf It is a quarterspan model, enter 0.25.\nIf not entered, defaults to halfspan", type=float, required = False)
parser.add_argument("-ext_converge", help="Ignore the presence of inteface planes and converge on my_total or fx_total (depending on AoA). If not entered, will converge on internal planes, if 2+ are present", action="store_true", required=False)
parser.add_argument("-converge_last_station", help="Flag to add ONLY the last (highest numbered) interface plane total pressure recovery to the convergence monitoring. If not entered, will add ALL interface planes total pressure recovery to the convergence monitoring.", action="store_true", required=False)

#PCM: FUTURE WORK-- add ability to print solution(.dat file) multiple times during run at a user specified interval. This would serve as a debugging functionality, allowing users to check the solution (i.e. shock development/placement) as it runs, instead of having to wait for the simulation to finish.
#parser.add_argument("-print_dat", help="Flag to print solution .dat files at intervals of print_dat.\nNOTE: Specifying too small of a value will greatly increase the time it takes to run your simulation AS WELL AS greatly increasing your data footprint.\nRECOMMENDED to not go below 500.\nIf not specified, the .dat file will ONLY be printed at the end of the run.", type=int, required=False)

## Global vars
global MRC_m, D_ref_m
global isBackPressured, backPressure_boundary
isBackPressured = False
backPressure_boundary = ''
global isInternal
global ext_converge
global converge_last_station

#global isDebugSoln, print_soln_iters

def Read_Inputs(currDir):
  global converge_last_station, ext_converge
  global isDebugSoln, print_soln_iters

  args = parser.parse_args()
  model_factor = args.model_factor
  case = args.case
  ext_converge = args.ext_converge
  converge_last_station = args.converge_last_station
  #print_dat = args.print_dat

  #if print_dat == None:
  #  isDebugSoln = False
  #else:
  #  isDebugSoln = True
  #  print_soln_iters = print_dat

  if model_factor == None:
    model_factor = 0.5

  if case == None:
    print('>> The case flag was not set, searching the current directory for the most recent mesh ...')
    mesh_files = []
    for file in os.listdir(currDir):
      temp = os.path.join(currDir, file)
      if temp.split('.')[-1] == 'msh':
        mesh_files.append(os.path.join(currDir,file))

    if not 'CFD' in currDir:
      if len(mesh_files) == 0:
        for dir in os.listdir(currDir):
          tempDir = os.path.join(currDir, dir)
          if os.path.isdir(tempDir):
            if str(dir).upper() == 'CFD':
              for file in os.listdir(tempDir):
                if file.split('.')[-1] == 'msh':
                  mesh_files.append(os.path.join(tempDir, file))

      if len(mesh_files) == 0:
        for dir in os.listdir(currDir):
          tempDir = os.path.join(currDir, dir)
          if os.path.isdir(tempDir):
            if str(dir).upper() == 'CFD':
              for mDir in os.listdir(tempDir):
                temp_mDir = os.path.join(tempDir, mDir)
                if os.path.isdir(temp_mDir):
                  if str(mDir).upper() == 'MESH':
                    for file in os.listdir(temp_mDir):
                      if file.split('.')[-1] == 'msh':
                        mesh_files.append(os.path.join(temp_mDir, file))
    else:
      if len(mesh_files) == 0:
        for dir in os.listdir(currDir):
          tempDir = os.path.join(currDir, dir)
          if os.path.isdir(tempDir):
            if str(dir).upper() == 'MESH':
              for file in os.listdir(tempDir):
                if file.split('.')[-1] == 'msh':
                  mesh_files.append(os.path.join(tempDir, file))

    if len(mesh_files) == 0:
      print('!! ERROR !! no .msh files were found in the current directory (\'' + str(currDir) + '\'). Exiting now.')
      sys.exit(0)
    else:
      case = max(mesh_files, key=os.path.getmtime)
      case_path = os.path.join(currDir, case)
      ## remove the full case path and only get the case name
      case = case.split('\\')[-1]

      case = case.replace('.msh','')
      print('>>>> The case file to be used is :\t' + str(case))

  else:
    caseFound = False
    for dir in os.listdir(currDir):
      temp = os.path.join(currDir, dir)
      if os.path.isdir(temp) and str(case) in temp:
        case_path = temp
        print('\n>>>> The case you entered was found. The code will grab data from within ' +str(case_path))
        caseFound = True

        if '_workflow_files' in case:
          case = case.replace('_workflow_files','')
        break

    if caseFound == False:
      print('!! ERROR !! The case provided was not found in the name of any of the mesh workflow file directories in the current directory. Exiting now.')
      sys.exit(0)


  return(case, case_path, model_factor)

def Get_Named_Selections_Mesh_2D(case_path):

  print('>> Running PyFluent to export Named Selections... (this will launch a seperate window to run fluent in)', flush=True)
  print('... this takes a minute ...', flush=True)

  ## this function will launch an instance of pyFluent to extract boundary conditions
  ## CURRENTLY, THERE IS NOT FUNCTIONALITY TO EXTRACT WHICH WALLS HAVE BOUNDARY LAYER CELLS APPLIED, SO BE PRE-EMPTIVE WHEN NAMING BOUNDARIES IN SPACECLAIM
  from ansys.fluent.core import launch_fluent
  # Launch Fluent in batch mode (no GUI)
  session = launch_fluent(
      dimension=2,
      precision="double",
      processor_count=2,
      start_transcript=False,
      ui_mode="no_gui_or_graphics"
  )

  session.file.read(file_type="mesh", file_name=case_path)

  ## if the mesh_summary file exists, delete it
  summaryFile = os.path.join(currDir, 'mesh_summary.txt')
  if os.path.exists(summaryFile):
    os.remove(summaryFile)

  # Redirect output to a file
  session.execute_tui("report/summary yes mesh_summary.txt")
  # Clean shutdown
  session.exit()

  if os.path.exists(summaryFile):
    print('>>>> PyFluent ran correctly... collecting Named Selections.')
  else:
    print('!! ERROR !! PyFluent did not export the mesh summary file. Exiting now.')
    sys.exit(0)

  data = []
  with open(summaryFile, 'r') as fin:
    for line in fin:
      if line.strip():
        data.append(line.strip().split())
  fin.close()

  ## for general tidiness, delete the mesh_summary file
  os.remove(summaryFile)

  Named_Selections = []
  walls_with_bl_cells = []
  for i in range(len(data)):
    try:
      if data[i][0] == 'Boundary' and data[i][1] == 'Conditions':
        j = i + 5
        while not data[j][0] == 'Setup':
          Named_Selections.append([data[j][0], data[j][2]])
          if data[j][2] == 'wall':
            walls_with_bl_cells.append(data[j][0])
          j += 1
        break

    except IndexError:
      continue

  ## When creating the journals, we dont use the axis, so remove any axis from the Named Selections
  Named_Selections = [bc for bc in Named_Selections if not 'axis' in bc[1]]

  return(Named_Selections, walls_with_bl_cells)

def Get_Aero_Prop_Interface_Surfs_2D(Named_Selections, walls_with_bl_cells):
  global isBackPressured
  global isInternal
  global backPressure_boundary
  global ext_converge

  ## seperate the Named Selections into the aero and prop surfaces
  surface_list_final = []
  interface_surfs = []
  farfield_surfs = []

  farfield_included = False ## for internal cases with set inflow and outflow conditions, there may not be a farfield
  for bc in range(len(Named_Selections)):
    if 'farfield' in Named_Selections[bc][0].lower() or Named_Selections[bc][1] == 'pressure-far-field':
      farfield_included = True
      farfield_surfs.append(Named_Selections[bc][0])
    if 'wall_aero' in Named_Selections[bc][0].lower():
      surface_list_final.append(Named_Selections[bc][0])
    elif 'interface' in Named_Selections[bc][0].lower() or Named_Selections[bc][1] == 'interface':
      if not Named_Selections[bc][0] in walls_with_bl_cells:
        interface_surfs.append(Named_Selections[bc][0])
    elif 'prop' in Named_Selections[bc][0]:
      if 'wall' in Named_Selections[bc][0]:
        surface_list_final.append(Named_Selections[bc][0])
    elif 'outlet' in Named_Selections[bc][0] or Named_Selections[bc][1] == 'pressure-outlet':
      isBackPressured = True
      backPressure_boundary = Named_Selections[bc][0]

  ## Remove any interfaces not of interest for setting up reports
  interface_surfs_for_reports = []
  surfs_to_exclude = ['farfield', 'boi', 'far', 'shock', 'capture']
  for i_surf in range(len(interface_surfs)):
    if not any(bad_surf in interface_surfs[i_surf] for bad_surf in surfs_to_exclude):
      interface_surfs_for_reports.append(interface_surfs[i_surf])

  interface_surfs_final = []
  if len(interface_surfs_for_reports) == 0:
    print('~~ NOTE ~~ the mesh does not have interfaces along station planes, and therefore no station quantities will be reported...')
  else:
    ## in 2D meshing, interfaces will ALWAYS be duplicated. Grab the interface with _2 and remove the interface with _1
    ## 2 has the proper normals, and will give positive mdots
    for i_surf in range(len(interface_surfs_for_reports)):
      interface_name = str(interface_surfs_for_reports[i_surf]).split('_')
      try:
        int_num = int(interface_name[-1])
        if int_num == 2:
          station_num = interface_name[-2]
          if 'station' in station_num:
            station_num = station_num.replace('station','')
          elif 'sta' in station_num:
            station_num = station_num.replace('sta','')
          elif 'STATION' in station_num:
            station_num = station_num.replace('STATION','')
          elif 'STA' in station_num:
            station_num = station_num.replace('STA', '')
          elif 'st' in station_num:
            station_num = station_num.replace('st','')
          elif 'ST' in station_num:
            station_num = station_num.replace('ST','')
          interface_surfs_final.append([interface_surfs_for_reports[i_surf], station_num])
      except ValueError:
        int_num = list(interface_name[-1])[-1]
        if int_num == 2:
          station_num = interface_name[-1][:-1]
          if 'station' in station_num:
            station_num = station_num.replace('station','')
          elif 'sta' in station_num:
            station_num = station_num.replace('sta','')
          elif 'STATION' in station_num:
            station_num = station_num.replace('STATION','')
          elif 'STA' in station_num:
            station_num = station_num.replace('STA', '')
          elif 'st' in station_num:
            station_num = station_num.replace('st','')
          elif 'ST' in station_num:
            station_num = station_num.replace('ST','')
          interface_surfs_final.append([interface_surfs_for_reports[i_surf], station_num])

  ## set the global flags
  isInternal = False
  if isBackPressured == True:
    isInternal = True
    if ext_converge == True:
      print('\n~~ NOTE ~~ You passed the ext_converge flag, but the solution is being back pressured. ext_converge is being set to False. The solution will use mdots through all station planes as the convergence criteria.')
      ext_converge = False
  elif len(interface_surfs_final) > 1 and ext_converge == False:
    isInternal = True

  return(surface_list_final, interface_surfs_final, farfield_included, farfield_surfs)


def Write_Farfield_Reports(reports_journal, farfield_included, farfield_surfs):
  report_names = []

  farfield_reports = ['p_0','pt_0','mach_0','t_0','tt_0','rho_0','q_0','v_mag_0']
  farfield_report_fields = ['pressure','total-pressure','mach-number','temperature','total-temperature','density','dynamic-pressure','velocity-magnitude']

  if farfield_included == True:
    for i_report in range(len(farfield_reports)):
        tui.make_report(reports_journal,report_name=farfield_reports[i_report],report_type='surface-massavg',report_surfaces=[farfield_surfs[0]],report_field=farfield_report_fields[i_report])
        report_names.append(farfield_reports[i_report])
    reports_journal.writelines('\n')

  return(report_names)

def Write_FaM_Reports(reports_journal, surface_list, report_names, interface_surfs):
  global MRC_m
  global ext_converge

  ## entries that span multiple bodies can become duplicated. remove them here.
  surfList_duplicates = set() ## initially a set to ensure no duplicates
  surfList_final = []
  for i_name in range(len(surface_list)):
    if '-' in surface_list[i_name]:
      realName = surface_list[i_name].split('-')[0]
      for i_dup_name in range(len(surface_list)):
        if not i_name == i_dup_name:
          realName_dup = surface_list[i_dup_name].split('-')[0]
          if realName == realName_dup:
            setLength = len(surfList_duplicates)
            surfList_duplicates.add(realName + '*')
            if not len(surfList_duplicates) == setLength:
              surfList_final.append(realName + '*')
    else:
      surfList_duplicates.add(surface_list[i_name])
      surfList_final.append(surface_list[i_name])

  force_vector_list = [[1,0],[0,1]]
  force_vector_names = ['fx','fy']

  moment_vector_list = [[1,0,0],[0,1,0],[0,0,1]]
  moment_vector_names = ['mx','my','mz']

  for i_force in range(len(force_vector_list)):
      for i_surf in range(len(surfList_final)):
          if '*' in surfList_final[i_surf]:
            report_name = force_vector_names[i_force] + '_' + surfList_final[i_surf][:-1]
          else:
            report_name = force_vector_names[i_force] + '_' + surfList_final[i_surf]

          tui.make_force_report(reports_journal,report_name,force_vector_list[i_force],[surfList_final[i_surf]])
          report_names.append(report_name)
          reports_journal.writelines('\n')

          if '*' in surfList_final[i_surf]:
            report_name = moment_vector_names[i_force] + '_' + surfList_final[i_surf][:-1]
          else:
            report_name = moment_vector_names[i_force] + '_' + surfList_final[i_surf]
          tui.make_moment_report(reports_journal,report_name,moment_vector_list[i_force],MRC_m[0:2],[surfList_final[i_surf]])
          report_names.append(report_name)
          reports_journal.writelines('\n')

  return(report_names)


def Write_Station_Reports(reports_journal, interface_surfs, report_names):
  global isInternal

  report_names_int = []

  ## Write Stream thrust reports
  #these are quantities that will never change for thrust drag, note only stations aligned with x work nx/ny/nz coded to 1,0,0

  '''
  stream_thrust_function_names = ['nx','ny','p_minus_p_0x','p_minus_p_0y','u_dot_n','rho_u_dot_n','rho_u_dot_n_ux', 'rho_u_dot_n_uy', 'p_minus_p_0x_plus_rho_u_dot_n_ux', 'p_minus_p_0y_plus_rho_u_dot_n_uy']
  stream_thrust_functions = ['1','0', '(pressure-p_0<rd>)*nx','(pressure-p_0<rd>)*ny','nx*x_velocity+ny*y_velocity',\
                           'density*u_dot_n','rho_u_dot_n*x_velocity','rho_u_dot_n*y_velocity', 'p_minus_p_0x+rho_u_dot_n_ux','p_minus_p_0y+rho_u_dot_n_uy']

  for i_function in range(len(stream_thrust_function_names)):
      tui.make_field_function(reports_journal,stream_thrust_function_names[i_function],stream_thrust_functions[i_function])

  functions_to_integrate = ['p_minus_p_0x_plus_rho_u_dot_n_ux', 'p_minus_p_0y_plus_rho_u_dot_n_uy']

  for i_surface in range(len(interface_surfs)):
      for i_function in range(len(functions_to_integrate)):
          report_name = 'int_' + functions_to_integrate[i_function] + '_' + interface_surfs[i_surface][1]
          tui.make_report(reports_journal,report_name=report_name,report_type='surface-integral',report_surfaces=[interface_surfs[i_surface][0]],report_field=functions_to_integrate[i_function])
          report_names_int.append(report_name)

  reports_journal.writelines('\n')
  '''

  #make station quantities
  station_reports = ['A','mdot','p','pt','mach','t','tt','rho','q','v_mag']
  station_report_fields = ['','','pressure','total-pressure','mach-number','temperature','total-temperature','density','dynamic-pressure','velocity-magnitude']

  mdot_reports = []
  for i_station in range(len(interface_surfs)):
      for i_report in range(len(station_reports)):
          if 'A' == station_reports[i_report]:
              report_type='surface-area'
          elif 'mdot' == station_reports[i_report]:
              report_type='surface-massflowrate'
          else:
              report_type='surface-massavg'

          tui.make_report(reports_journal,report_name=station_reports[i_report]+'_'+interface_surfs[i_station][1],report_type=report_type,report_surfaces=[interface_surfs[i_station][0]],report_field=station_report_fields[i_report])
          report_names.append(station_reports[i_report]+'_'+interface_surfs[i_station][1])
          if 'mdot' == station_reports[i_report]:
            mdot_reports.append(station_reports[i_report]+'_'+interface_surfs[i_station][1])

  reports_journal.writelines('\n')

  return(report_names, report_names_int, mdot_reports)

def Create_Report_Files( reports_journal, report_names, report_names_int):

  reports = ['reports', 'reports1']
  ## create the reports file to be averaged every step

  tui.make_report_file(reports_journal, 1, report_names,'reports')
  reports_journal.writelines('\n')

  ## create the reports file containing integrations
  #tui.make_report_file(reports_journal, 100, report_names_int, 'reports1')
  #reports_journal.writelines('\n')

  return(reports)


#############################################################################
def Write_Sweep_Journal_2D(reports_journal, sweep_journal, currDir, reports, case, sweep_vars, surface_list, interface_surfs, model_factor, mdot_reports, farfield_surfs, reportsDir):

  global isBackPressured, isInternal
  global converge_last_station, ext_converge
  global backPressure_boundary
  global MRC_m
  global D_ref_m
  S_ref_m2 = np.pi*D_ref_m**2/4

  I1_Mach = sweep_vars['Mach']
  I2_h_ft = sweep_vars['h_ft']
  I3_BP_Pa = sweep_vars['BackPress_Pa']

  sweep_type = 'Mach'
  run_start_no = 2
  point_start_no = 1
  run_no = run_start_no
  point_no = point_start_no

  TI = 0.001 #turbulence intensity, fraction
  TVR = 10 #turbulent viscosity ratio, fraction

  Cmb_stop = 0.005
  previous_iterations = 50

  if isBackPressured == True:
    nRamps = 6 ## Depending on the back pressure, this will have to be variable, as well as the number of iterations per step
    backPress_restartFile = ''

  # %%
  case_name_list = []

  tui.set_tui_version(sweep_journal)
  tui.set_batch_options(sweep_journal)
  tui.read_case(sweep_journal,case)

  metadata_dict = {}
  metadata_fields = ['Case file','Sweep','Run','Point','mach_0_BC','h_ft','p_0_BC','q_0_BC','t_0_BC','AoA_deg','AoSw_deg','BackPress_Pa']
  for i_field in range(len(metadata_fields)):
    metadata_dict[metadata_fields[i_field]] = []

  for i_Mach in range(len(I1_Mach)):
    #for i_h in range(len(I2_h_ft)):
    i_h = i_Mach


    #if I1_Mach[i_Mach] == 2.25:
    #  I3_BP_Pa = sweep_vars['BackPress_Pa'][0]
    #elif I1_Mach[i_Mach] == 2.5:
    #  I3_BP_Pa = sweep_vars['BackPress_Pa'][1]


    airspeed = Airspeed(I2_h_ft[i_h],airspeed_value=I1_Mach[i_Mach],airspeed_type='Mach',units='metric').get_results() #get freestream conditions
    p0_Pa = airspeed['pressure']
    T0_K = airspeed['temperature']
    rho_kgpm3 = airspeed['density']
    V_mps = airspeed['trueAirspeed_mps']
    pt0_Pa = airspeed['totalPressure']
    q0_Pa = airspeed['dynamicPressure_Pa']

    for i_BP in range(len(I3_BP_Pa)):

      ## Set the farfield conditions
      wind_vector_c = np.matmul(roty(180),np.matmul(roty(0).transpose(),np.matmul(rotz(0),np.array([-1,0,0]).transpose())))
      for i_farField in range(len(farfield_surfs)):
        tui.set_pressure_farfield(sweep_journal,boundary_name=farfield_surfs[i_farField] ,mode='axisymmetric',turbulence='SST',Mach=I1_Mach[i_Mach],p_Pa=p0_Pa,T_K=T0_K,TI=TI,TVR=TVR,x_component_direction = wind_vector_c[0],y_component_direction=wind_vector_c[1]) #set characteristic boundary

      tui.set_reference_values(sweep_journal,boundary_name=farfield_surfs[0])

      ###########################
      ## Set up the convergence criteria(s). ON first case of sweep, create the convergence criteria.
      ###########################
      if all(ind == 0 for ind in [i_Mach,i_h, i_BP]): #if first case

        ## set the views for the contours
        set_view_for_contours(sweep_journal, currDir)

        if isBackPressured == True or (isInternal == True and ext_converge == False):
          if isBackPressured == True:
            if I3_BP_Pa[i_BP] <= 85000:
              iters_to_ignore = 1050 ## If the solution is back pressured, back pressure/cfl ramping will always start at 550 iters
            else:
              iters_to_ignore = 1550 ## If the solution is back pressured, back pressure/cfl ramping will always start at 550 iters

            ## only if back pressuring, increase the continuity residual slightly
            tui.modify_residual_criterions_2D(sweep_journal, 0.0025, 0.001, 0.001, 0.000001, 0.001, 0.001)

          if converge_last_station == False:
            for i_ptr in (mdot_reports):
              criterion = i_ptr
              tui.make_convergence_criterion(sweep_journal, criterion, criterion, previous_iterations, Cmb_stop, active=True, ignore_iterations_before=iters_to_ignore)
          else:
            criterion = mdot_reports[-1]
            tui.make_convergence_criterion(sweep_journal, criterion, criterion, previous_iterations, Cmb_stop, active=True, ignore_iterations_before=iters_to_ignore)
        else: ## external flow or ext_converge flag is set
          ## initialize both the my_total and the fx_total convergence criteria. The Fx will be set to inactive, until a 0 AoA case is reached.
          fx_criterion = 'fx_total'
          tui.make_convergence_criterion(sweep_journal, fx_criterion, fx_criterion, previous_iterations, Cmb_stop, active=True)
      ###########################
      ## END set convergence metrics
      ###########################

      if isBackPressured:
        case_name_full = (case + '_Mach_' + str(I1_Mach[i_Mach]) + '_h_' + str(I2_h_ft[i_h]) + '_BP_' + str(I3_BP_Pa[i_BP]) + '_Pa').replace('.','p')
      else:
        case_name_full = (case + '_Mach_' + str(I1_Mach[i_Mach]) + '_h_' + str(I2_h_ft[i_h])).replace('.','p')

      ## add the created case to the case list
      case_name_list.append(case_name_full)

      ## update the report file(s)
      #for i_report in range(len(reports)):
      for i_report in range(1):
        tui.change_all_report_file_names(sweep_journal, case_name_full, reports[i_report])

      ## if there is a pressure_outlet, set it up
      if isBackPressured == True:
        ## If gamma != 1.4, pass in gamma as 3rd argument in Calculate_Total_temp
        tui.set_pressure_outlet_2D(sweep_journal, 0 , Calculate_Total_Temp(T0_K, I1_Mach[i_Mach]))

      ## Initialize the solution
      tui.hyb_init(sweep_journal)
      tui.fmg_init(sweep_journal)


      if isBackPressured:

        ## Set CFL to 10, iterate for 50 iters
        tui.set_CFL(sweep_journal, CFL=10)
        tui.solve(sweep_journal,iterations=50)
        ## Ramp CFL up to 40, iterate for rest of simulation
        tui.set_CFL(sweep_journal, CFL=40)
        tui.solve(sweep_journal,iterations=500)
        if I3_BP_Pa[i_BP] <= 50000:
          tui.set_pressure_outlet_2D(sweep_journal, I3_BP_Pa[i_BP] , Calculate_Total_Temp(T0_K, I1_Mach[i_Mach]))
          tui.solve(sweep_journal, iterations = 1000)
        elif I3_BP_Pa[i_BP] <= 85000:
          for i_curr_bp in np.linspace(50000, I3_BP_Pa[i_BP], 5):
            tui.set_pressure_outlet_2D(sweep_journal, round(i_curr_bp,2) , Calculate_Total_Temp(T0_K, I1_Mach[i_Mach]))
            tui.solve(sweep_journal, iterations = 100)
          tui.solve(sweep_journal, iterations = 1000)
        else:
          for i_curr_bp in np.linspace(50000, I3_BP_Pa[i_BP], 10):
            tui.set_pressure_outlet_2D(sweep_journal, round(i_curr_bp,2) , Calculate_Total_Temp(T0_K, I1_Mach[i_Mach]))
            tui.solve(sweep_journal, iterations = 100)
          tui.solve(sweep_journal, iterations = 1500)

      else:
        ## Set CFL to 10, iterate for 50 iters
        tui.set_CFL(sweep_journal, CFL=10)
        tui.solve(sweep_journal,iterations=50)
        ## Ramp CFL up to 40, iterate for rest of simulation
        tui.set_CFL(sweep_journal, CFL=40)
        tui.solve(sweep_journal, iterations = 1950)

      ## save the final solution file
      dataDir = os.path.join(currDir, 'CASE_AND_DATA_FILES')
      if not os.path.exists(dataDir):
        os.mkdir(dataDir)
      caseData = os.path.join(dataDir, case_name_full)
      caseData = caseData.replace('\\', '\\\\')
      caseData = caseData
      tui.write_case_data(sweep_journal,case_name=caseData)

      ## Create the images
      create_contours(sweep_journal, currDir, case_name_full)

      tui.write_converged_reports(sweep_journal, reports[0])

      sweep_journal.writelines('\n')

      metadata_dict['Case file'].append(case_name_full)
      metadata_dict['mach_0_BC'].append(I1_Mach[i_Mach])
      metadata_dict['h_ft'].append(I2_h_ft[i_h])
      metadata_dict['p_0_BC'].append(p0_Pa)
      metadata_dict['q_0_BC'].append(q0_Pa)
      metadata_dict['t_0_BC'].append(T0_K)

      metadata_dict['AoA_deg'].append(0)
      metadata_dict['AoSw_deg'].append(0)

      metadata_dict['BackPress_Pa'].append(I3_BP_Pa[i_BP])

      metadata_dict['Sweep'].append(sweep_type)
      metadata_dict['Run'].append(run_no)
      metadata_dict['Point'].append(point_no)

      if sweep_type == 'BackPress':
        point_no += 1
    if sweep_type == 'Mach':
      point_no+=1
    else:
      run_no+=1

  tui.exit_job(sweep_journal)

  ## Before exiting this function, create the ouput directory and write the required files
  ## create the ouput directory
  outputDirName = os.path.join(currDir, case + '_POST')
  if not os.path.exists(outputDirName):
    os.mkdir(outputDirName)

  ## put all the files to run post in a single folder for general tidiness
  outputDirFiles = os.path.join(outputDirName, 'case_files')
  if not os.path.exists(outputDirFiles):
    os.mkdir(outputDirFiles)

  ## create the metadata file
  metaFile = os.path.join(outputDirFiles, 'METADATA_FOR_POST.json')
  if os.path.exists(metaFile):
    os.remove(metaFile)

  with open(metaFile, 'w') as fout:
    json.dump(metadata_dict, fout)
  fout.close()

  ## create the case file / geometry file
  caseFile = os.path.join(outputDirFiles, 'CASE_DATA_FOR_POST.TXT')
  if os.path.exists(caseFile):
    os.remove(caseFile)

  with open(caseFile, 'w') as fout:
    fout.write('Base_case: ' + str(case) + '\n')
    fout.write('Journal_name:\t' + str(sweep_journal.name).split('\\')[-1] + '\n')
    fout.write('Case_name_list: ')
    for case_name in case_name_list:
      fout.write(case_name + ' ')
    fout.write('\n')
    fout.write('MRC_m: ' + str(MRC_m[0]) + ' ' +str(MRC_m[1]) + ' ' + str(MRC_m[2]) + '\n')
    fout.write('L_ref: ' + str(D_ref_m) + '\n')
    fout.write('S_ref: ' +str(S_ref_m2) + '\n')
    fout.write('model_factor: ' + str(model_factor)+ '\n')
    fout.write('isBackPressured: ' + str(isBackPressured))
  fout.close()

  ## create the surfaces file
  surfFile = os.path.join(outputDirFiles, 'SURFACES_FOR_POST.txt')
  if os.path.exists(surfFile):
    os.remove(surfFile)

  with open(surfFile, 'w') as fout:
    fout.write('walls: ')
    for surf in surface_list:
      fout.write(str(surf) + ' ')
    fout.write('\ninterfaces: ')
    for interface in interface_surfs:
      fout.write(str(interface[0]) + ' ')
  fout.close()

def Create_Reports_Dir(currDir):
  reportsDir = os.path.join(currDir, 'REPORT_FILES')
  if not os.path.exists(reportsDir):
    os.mkdir(reportsDir)

  return(reportsDir)

def Create_Journal(currDir, case, journal_type):

  if journal_type == 'reports':
    case = 'MAKE_REPORTS_' + case
  elif journal_type == 'sweep':
    case = case + '_SWEEP'

  journal_file = os.path.join(currDir, case)
  reports_journal = open(journal_file+'.jou', 'w')
  return(reports_journal)

def Close_Journal(journal):
  journal.close()

if __name__ == '__main__':

  ## VARIABLES TO SET
  ## Geometry variables
  MRC_m = [0.915,0,0]
  D_ref_m = 0.155

  ct_fin_m = 0.1063692
  cr_fin_m = 0.1661163
  full_fin_span_m = 0.175
  LE_sweep_fin_deg = 30
  nose_to_fin_LE_m = 2.609575

  ct_strake_m = 0.6612787
  cr_strake_m = 0.8128001
  full_strake_span_m = 0.3074
  LE_sweep_strake_deg = 10
  nose_to_strake_LE_m = 0.752

  ## Sweep Variables:
  sweep_vars = {}
  sweep_vars['Mach'] = [2.5]
  sweep_vars['h_ft'] = [10000]
  #sweep_vars['BackPress_Pa'] = [list(np.linspace(660000, 665000, 6)), list(np.linspace(950000, 970000, 21))]  ## This will ONLY get applied if you have a boundary called 'pressure_outlet'.
  sweep_vars['BackPress_Pa'] = [971000, 972000, 973000, 974000, 975000, 976000, 978000, 979000, 980000]
  ## END OF VARIABLES TO SET

  ################################
  ## GRAB VALUES FROM MESH FILES
  ################################
  ## get the working directory
  currDir = os.getcwd()

  ## read command line arguments
  case, case_path, model_factor = Read_Inputs(currDir)

  ## grab Named Selections from Mesh
  Named_Selections, walls_with_bl_cells = Get_Named_Selections_Mesh_2D(case_path)

  ## Collect the aero and prop surfaces
  surface_list, interface_surfs, farfield_included, farfield_surfs = Get_Aero_Prop_Interface_Surfs_2D(Named_Selections, walls_with_bl_cells)

  ## create the reports directory where the reports will be printed to
  reportsDir = Create_Reports_Dir(currDir)
  ################################
  ## START REPORTS
  ################################
  ## Create the reports_journal object
  reports_journal = Create_Journal(currDir, case, 'reports')

  ## write the farfield journals
  report_names = Write_Farfield_Reports(reports_journal, farfield_included, farfield_surfs)

  ## write the Force/moment report entries for all the walls
  report_names = Write_FaM_Reports(reports_journal, surface_list, report_names, interface_surfs)

  ## write the stream thrust report entries
  report_names, report_names_int, mdot_reports = Write_Station_Reports(reports_journal, interface_surfs, report_names)

  ## write the total pressure recovery report entries
  ## PCM: NEED TO FIND A WAY TO ADD PTR. Cant be done from TUI, because station_report definitions arent stored in global memory.
  #report_names, report_names_ptr = Write_Total_Pressure_Recovery_Reports(reports_journal, interface_surfs, report_names)

  ## Create the report files
  reports = Create_Report_Files(reports_journal, report_names, report_names_int)

  ################################
  ## START SWEEP
  ################################
  ## Create the sweep journal
  sweep_journal = Create_Journal(currDir, case, 'sweep')

  Write_Sweep_Journal_2D(reports_journal, sweep_journal, currDir, reports, case, sweep_vars, surface_list, interface_surfs, model_factor, mdot_reports, farfield_surfs, reportsDir)

  ## Close the reports journal
  Close_Journal(reports_journal)

  ## Close the sweep journal
  Close_Journal(sweep_journal)




