# -*- coding: utf-8 -*-
"""
@author: pmills
"""

import sys
import os
import numpy as np
import argparse
import subprocess

from airspeed_v1 import Airspeed
from calculateAtmosphere_v1 import *
from unitConverter import *

from writeSpaceClaimFile import Create_SpaceClaim_File

## Add tools folder to system path
pathToEngineeringTools = 'G:\\My Drive\\tools'
if pathToEngineeringTools not in sys.path:
  sys.path.append(pathToEngineeringTools)


########## GLOBAL VARIABLES TO EDIT ##########
global mesh_settings
mesh_settings = {}
global mach, alt
global nPrismCells
global flowthrough
global fidelity

########## DONT EDIT ANYTHING BELOW THIS LINE ##########

## setup the argparser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("-case", help="Name of the SpaceClaim file you want to read. If not specified, will default to most recently created SpaceClaim file in the current working directory.", type=str, required=False)
parser.add_argument("-fidelity", help="desired fidelity level of mesh.\nAcceptable inputs: tiny, coarse, medium, fine, extra_fine, ultra_fine.\nIf not specified, defaults to medium.\n", required=False, action = 'store', default= 'medium')
parser.add_argument("-nPrismCells", help="desired number of prism cells.\nIf not specified, defaults to 45.", type=int, required=False)
parser.add_argument("-mach", help="mach values(s) you want to run at. (Used for the Re calculation -> y+ calculation).\nTo enter multiple mach values: -mach 2.75 3.0 3.5" ,nargs="+", type=float, required=True)
parser.add_argument("-alt", help="altitude (ft) you want to run at. (Used for the Re calculation -> y+ calculation).\nTo enter multiple alt values: -alt 10000 65000 65000." ,nargs="+", type= float, required = True)
parser.add_argument("-flowthrough", help="flag when you are running a full flowthrough model. Will shift down fidelity of cells across AIP and cells across annulus accordingly.\nDefault is NOT flowthrough.", action="store_true", required = False)


def Read_Inputs(currDir):
  global mesh_settings
  global mach, alt
  global nPrismCells
  global flowthrough
  global fidelity

  args = parser.parse_args()
  case = args.case
  fidelity = args.fidelity
  nPrismCells = args.nPrismCells
  mach = args.mach
  alt = args.alt
  flowthrough = args.flowthrough

  if not len(mach) == len(alt):
    if not (len(mach) == 1 or len(alt) == 1):
      print('!! ERROR !! There was ' + str(len(mach)) + ' mach values entered and ' + str(len(alt)) + ' altitude values entered. Mach and alt MUST be the same length, OR 1 mach for n altitude OR 1 altitude for n mach. Exiting now...')
      sys.exit(0)

  if nPrismCells == None:
    nPrismCells = 45


  CADdirFound = {}
  CADdirFound['found'] = False
  CADdirFound['path'] = ''
  if case == None:
    print('>> The case flag was not set, searching the current directory for the most recent .SCDOC file...')
    sc_files = []
    for file in os.listdir(currDir):
      temp = os.path.join(currDir, file)
      if os.path.isfile(temp) and os.path.splitext(temp)[1] == '.scdoc':
        sc_files.append(file)
        CADdirFound['found'] = True
        CADdirFound['path'] = temp

    if len(sc_files) == 0:
      for dir in os.listdir(currDir):
        tempDir = os.path.join(currDir, dir)
        if os.path.isdir(tempDir):
          if str(dir) == 'CAD':
            for file in os.listdir(tempDir):
              tempFile = os.path.join(tempDir, file)
              if os.path.isfile(tempFile) and os.path.splitext(tempFile)[1] == '.scdoc':
                CADdirFound['found'] = True
                CADdirFound['path'] = tempDir
                sc_files.append(os.path.join(tempDir, file))
    if len(sc_files) == 0:
      print('>>>> ERROR... no .SCDOC files found in the currrent directory(\''+str(currDir)+'\'). Exiting now')
      sys.exit(0)
    else:
      case = max(sc_files, key=os.path.getmtime)
      case_path = os.path.join(currDir, case)
      print('>>>> The .SCDOC file to be used is:\t'+str(case))
      case = case.replace('.scdoc', '')


  else:
    caseFound = False
    for dir in os.listdir(currDir):
      temp = os.path.join(currDir, dir)
      if os.path.isfile(temp) and str(case) in temp:
        case_path = temp
        print('\n>> The case you entered was found. The code will grab data from within ' +str(case_path))
        caseFound = True

    if caseFound == False:
      print('!! ERROR. The case provided was not found in the name of any of the .scdoc files in the current directory. Exiting now.')
      sys.exit(0)

  match fidelity:
    case None:
      print('!! The fidelity flag was not passed... setting fidelity to \'medium\'...')
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 1.0
      mesh_settings['curve_normal_angle'] = 7.50
      mesh_settings['surf_mesh_space'] = 75.0
      mesh_settings['growth_rate'] = 1.20
      mesh_settings['num_cells_per_gap'] = 3.0
      mesh_settings['element_type'] = 'tetrahedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 12.5
        mesh_settings['internal_num_cell_per_annulus'] = 4
      else:
        mesh_settings['internal_num_cell_per_diam'] = 25
        mesh_settings['internal_num_cell_per_annulus'] = 7

    case 'tiny':
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 100.0
      mesh_settings['curve_normal_angle'] = 12.50
      mesh_settings['surf_mesh_space'] = 25.0
      mesh_settings['growth_rate'] = 1.20
      mesh_settings['num_cells_per_gap'] = 2.0
      mesh_settings['element_type'] = 'polyhedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 0
        mesh_settings['internal_num_cell_per_annulus'] = 0
      else:
        mesh_settings['internal_num_cell_per_diam'] = 0
        mesh_settings['internal_num_cell_per_annulus'] = 0

    case 'coarse':
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 30.0
      mesh_settings['curve_normal_angle'] = 10.0
      mesh_settings['surf_mesh_space'] = 50.0
      mesh_settings['growth_rate'] = 1.2
      mesh_settings['num_cells_per_gap'] = 2.0
      mesh_settings['element_type'] = 'polyhedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 0
        mesh_settings['internal_num_cell_per_annulus'] = 0
      else:
        mesh_settings['internal_num_cell_per_diam'] = 12.5
        mesh_settings['internal_num_cell_per_annulus'] = 4

    case 'medium':
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 1.0
      mesh_settings['curve_normal_angle'] = 7.5
      mesh_settings['surf_mesh_space'] = 75.0
      mesh_settings['growth_rate'] = 1.20
      mesh_settings['num_cells_per_gap'] = 3.0
      mesh_settings['element_type'] = 'tetrahedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 12.5
        mesh_settings['internal_num_cell_per_annulus'] = 4
      else:
        mesh_settings['internal_num_cell_per_diam'] = 25
        mesh_settings['internal_num_cell_per_annulus'] = 7

    case 'fine':
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 1.0
      mesh_settings['curve_normal_angle'] = 5.0
      mesh_settings['surf_mesh_space'] = 100.0
      mesh_settings['growth_rate'] = 1.15
      mesh_settings['num_cells_per_gap'] = 5.0
      mesh_settings['element_type'] = 'tetrahedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 25
        mesh_settings['internal_num_cell_per_annulus'] = 7
      else:
        mesh_settings['internal_num_cell_per_diam'] = 37.5
        mesh_settings['internal_num_cell_per_annulus'] = 11

    case 'extra_fine':
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 0.87
      mesh_settings['curve_normal_angle'] = 3.75
      mesh_settings['surf_mesh_space'] = 115.0
      mesh_settings['growth_rate'] = 1.10
      mesh_settings['num_cells_per_gap'] = 7.0
      mesh_settings['element_type'] = 'tetrahedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 37.5
        mesh_settings['internal_num_cell_per_annulus'] = 11
      else:
        mesh_settings['internal_num_cell_per_diam'] = 50
        mesh_settings['internal_num_cell_per_annulus'] = 14

    case 'ultra_fine':
      mesh_settings['fidelity'] = fidelity
      mesh_settings['y_plus'] = 0.76
      mesh_settings['curve_normal_angle'] = 2.5
      mesh_settings['surf_mesh_space'] = 132.0
      mesh_settings['growth_rate'] = 1.10
      mesh_settings['num_cells_per_gap'] = 9.0
      mesh_settings['element_type'] = 'tetrahedra'
      if flowthrough == True:
        mesh_settings['internal_num_cell_per_diam'] = 50
        mesh_settings['internal_num_cell_per_annulus'] = 14
      else:
        mesh_settings['internal_num_cell_per_diam'] = 62.5
        mesh_settings['internal_num_cell_per_annulus'] = 18
    case _:
      print('** INVALID ENTRY FOR FIDELITY ** To check acceptable fidelity entries->  mesh_setup.py -h')
      print('Exiting now.')
      sys.exit(0)

  return(case, case_path, CADdirFound)


def Run_Space_Claim(currDir, SC_File):

  ## spaceclaim requires files to have double backslashes and be enclosed by ""
  SC_File = '\"' + SC_File.replace('\\', '\\\\') + '\"'

  print('>> Running SpaceClaim in batch mode to extract geometries...\n')
  print('(this may take a couple of minutes...)\n')

  sc_exe = '\"C:\\Program Files\\ANSYS Inc\\v242\\SCDM\\SpaceClaim.exe\" /RunScript=' + str(SC_File) +' /ScriptAPI=242 /Headless=True /Splash=False /Welcome=False /ExitAfterScript=True'
  process = subprocess.Popen(sc_exe, stdout=subprocess.PIPE)
  output, error = process.communicate()

  ## check for errors
  if error == None:
    print('>>>> SpaceClaim ran successfully...')
  else:
    print('>>>> SpaceClaim failed with error:\n')
    print(error)
    print('\n\n Exiting now.')
    sys.exit(0)

  geomFileFound = False
  meshFileFound = False
  cfd_dir = os.path.join(currDir, 'CFD')
  if os.path.exists(cfd_dir):
    sc_out_dir = os.path.join(cfd_dir, 'MESH_INPUTS')
    if os.path.exists(sc_out_dir):
      for file in os.listdir(sc_out_dir):
        if file == 'mesh_vars.txt':
          meshFileFound = True
          meshFile = os.path.join(sc_out_dir, 'mesh_vars.txt')
        elif file == 'geometry_values.txt':
          geomFileFound = True
          geomFile = os.path.join(sc_out_dir, 'geometry_values.txt')
    else:
      print('!! ERROR !! The directory ' + str(sc_out_dir) + ' was not found. Exiting now.')
      sys.exit(0)
  else:
    print('!! ERROR !! The directory ' + str(cfd_dir) + ' was not found. Exiting now.')
    sys.exit(0)

  if geomFileFound == False:
    print('!! ERROR !! The geometry_values file was not found. Exiting now.')
    sys.exit(0)
  if meshFileFound == False:
    print('!! ERROR !! The mesh_vars files was not found. Exiting now.')
    sys.exit(0)

  return(meshFile, geomFile, sc_out_dir)

def Collect_SpaceClaim_Data(meshFile, geomFile):

  ## read the mesh_vars file
  data = []
  with open(meshFile, 'r') as fin:
    for line in fin:
      data.append(line.strip().split())
  fin.close()

  INTERFACES = []
  interfaceKeys = []
  BOIS = []
  interfaceFound = False
  for i in range(len(data)):
    interface = {}
    try:
      if data[i][0] == 'Body_of_Influence:':
        BOIS.append(data[i][1])
      elif data[i][0] == 'INTERFACE:':
        interface['interface'] = data[i][1]
        interface['diameter'] = float(data[i+1][1])
        interface['inner_diameter'] = float(data[i+2][1])
        INTERFACES.append(interface)

        interfaceKeys.append(data[i][1])
    except IndexError:
      continue

  if len(interfaceKeys) > 0:
    interfaceFound = True

  if interfaceFound == True:
    ## sort the interfaces in increasing order by station name.
    sortOrder = []
    for key in interfaceKeys:
      stationNum = ''
      if 'aip' == key.lower():
        stationNum = 999999.0
      elif 'station' in key.lower():
        stationNum = float(key.replace('station', '').replace('p','.').replace('_',''))
      elif 'sta' in key.lower():
        stationNum = float(key.replace('sta', '').replace('p','.').replace('_',''))
      elif 'st' in key.lower():
        stationNum = float(key.replace('st', '').replace('p','.').replace('_',''))

      if not stationNum == '':
        sortOrder.append(stationNum)
      else:
        print('~ NOTE ~ Station \'' + str(key) + '\' does not meet the naming conventions and will not be added to the sorting of the stations')

    sort_lists = sorted(zip(sortOrder, interfaceKeys))

    sortOrder, interfaceKeys = zip(*sort_lists)
    sortOrder = list(sortOrder)
    interfaceKeys = list(interfaceKeys)

    if sortOrder[-1] == 999999.0:
      sortOrder.pop(-1)
      aip = interfaceKeys.pop(-1)
      for i_sta in range(len(sortOrder)):
        if sortOrder[i_sta] == 1.0:
          sortOrder.insert(i_sta+1, 999999.0)
          interfaceKeys.insert(i_sta+1, aip)
          break

  ## read the geometry_values file
  data = []
  with open(geomFile, 'r') as fin:
    for line in fin:
      data.append(line.strip().split())
  fin.close()

  GEOM_VARS = []
  for i in range(len(data)):
    geom = {}
    try:
      if data[i][1] == 'canard' or 'fin' in data[i][1] or data[i][1] == 'strake':
        geom['surface'] = data[i][1]
        geom['root_chord'] = float(data[i+1][1])
        geom['tip_chord'] = float(data[i+2][1])
        geom['span'] = float(data[i+3][1])
        geom['LE_sweep_angle'] = float(data[i+4][1])
        geom['LE_distance_from_nose'] = float(data[i+5][1])
        GEOM_VARS.append(geom)
      elif data[i][1] == 'body':
        geom['surface'] = data[i][1]
        geom['diameter'] = float(data[i+1][1])
        geom['length'] = float(data[i+2][1])
        GEOM_VARS.append(geom)
    except IndexError:
      continue

  return(BOIS, INTERFACES, GEOM_VARS, interfaceFound, interfaceKeys)

def Calculate_Cbar(surf):
  root_chord = surf['root_chord']
  tip_chord = surf['tip_chord']
  taper_ratio = tip_chord/root_chord
  cbar = (2/3)*root_chord*(1+taper_ratio+taper_ratio**2)/(1+taper_ratio)

  return(cbar)


def Calculate_Mesh_Parameters(fout, bois, interfaces, geom_vars, interfaceFound, sortedInterfaceKeys):
  global mesh_settings
  global mach, alt

  s_surface = []
  x_surface = []
  min_yplus = []

  body_size_min = False

  ## check geom_vars for control surfs. IF present, calculate cbar
  for i_surf in range(len(geom_vars)):
    if geom_vars[i_surf]['surface'] == 'canard' or 'fin' in geom_vars[i_surf]['surface'] or geom_vars[i_surf]['surface'] == 'strake':
      cbar = Calculate_Cbar(geom_vars[i_surf])
      x_surface.append([geom_vars[i_surf]['surface'], cbar])
      s_surface.append([geom_vars[i_surf]['surface'], cbar/mesh_settings['surf_mesh_space']])
    elif geom_vars[i_surf]['surface'] == 'body':
      x_surface.append([geom_vars[i_surf]['surface'], geom_vars[i_surf]['length']])
      ## check curvature normal angle spacing against D_ref/no_cells, determine which is driving
      circ_body = geom_vars[i_surf]['diameter']*np.pi / (360/mesh_settings['curve_normal_angle'])
      s_surf_min = min(geom_vars[i_surf]['length']/mesh_settings['surf_mesh_space'], circ_body)
      body_size_min = s_surf_min
      s_surface.append([geom_vars[i_surf]['surface'], s_surf_min])

  fout.write('Mesh input variables for fidelity level:\t' + mesh_settings['fidelity'].upper() + '\n')
  fout.write('~'*30 + '\n')
  fout.write('y_plus:\t' + str(mesh_settings['y_plus']) + '\n')
  fout.write('Curvature normal angle:\t' + str(mesh_settings['curve_normal_angle']) + '\n')
  fout.write('Number of cells across trailing edges (cells per gap):\t' + str(mesh_settings['num_cells_per_gap']) + '\n')
  fout.write('Surface mesh spacing (num cells across Lref):\t' + str(mesh_settings['surf_mesh_space']) + '\n')
  fout.write('Growth rate:\t' + str(mesh_settings['growth_rate']) + '\n')
  fout.write('Propulsion flowpath volume refinement (num. cells across AIP diameter):\t' + str(mesh_settings['internal_num_cell_per_diam']) + '\n')
  fout.write('Element type:\t' + str(mesh_settings['element_type']) + '\n')
  fout.write('~'*30 + '\n\n')

  writeFlag = True
  for i_mach in range(len(mach)):

    if len(alt) == len(mach):
      fout.write('## MESH SETTINGS FOR MACH ' +str(mach[i_mach]) + ' ALTITUDE ' + str(alt[i_mach]) + ' (ft) .....................\n')
      min_y_inner = Calculate_Min_Y_Inner(meshFile, mach[i_mach], alt[i_mach], x_surface, s_surface, interfaces, sortedInterfaceKeys)
      for i_surf in range(len(x_surface)):
        min_yplus.append(Calculate_Mesh_Settings(fout, writeFlag, mach[i_mach], alt[i_mach], x_surface[i_surf], s_surface[i_surf], min_y_inner))
      if interfaceFound == True:
        min_yplus.append(Calculate_Internal_Mesh_Settings(fout, writeFlag, mach[i_mach], alt[i_mach], interfaces, sortedInterfaceKeys, min_y_inner))

      Write_BOI_Sizing(fout, body_size_min, bois)
      fout.write('## END MESH SETTINGS FOR MACH ' +str(mach[i_mach]) + ' ALTITUDE ' + str(alt[i_mach]) + ' (ft) .....................\n\n')

    elif len(alt) == 1:
      fout.write('## MESH SETTINGS FOR MACH ' +str(mach[i_mach]) + ' ALTITUDE ' + str(alt[0]) + ' (ft) .....................\n')
      min_y_inner = Calculate_Min_Y_Inner(meshFile, mach[i_mach], alt[0], x_surface, s_surface, interfaces, sortedInterfaceKeys)
      for i_surf in range(len(x_surface)):
        min_yplus.append(Calculate_Mesh_Settings(fout, writeFlag, mach[i_mach], alt[i_mach], x_surface[i_surf], s_surface[i_surf], min_y_inner))
      if interfaceFound == True:
        min_yplus.append(Calculate_Internal_Mesh_Settings(fout, writeFlag, mach[i_mach], alt[i_mach], interfaces, sortedInterfaceKeys, min_y_inner))

      Write_BOI_Sizing(fout, body_size_min, bois)
      fout.write('## END MESH SETTINGS FOR MACH ' +str(mach[i_mach]) + ' ALTITUDE ' + str(alt[0]) + ' (ft) .....................\n\n')

    else:
      for i_alt in range(len(alt)):
        fout.write('## MESH SETTINGS FOR MACH ' +str(mach[i_mach]) + ' ALTITUDE ' + str(alt[i_alt]) + ' (ft) .....................\n')
        min_y_inner = Calculate_Min_Y_Inner(meshFile, mach[i_mach], alt[i_alt], x_surface, s_surface, interfaces, sortedInterfaceKeys)
        for i_surf in range(len(x_surface)):
          min_yplus.append(Calculate_Mesh_Settings(fout, writeFlag, mach[i_mach], alt[i_alt], x_surface[i_surf], s_surface[i_surf], min_y_inner))
        if interfaceFound == True:
          min_yplus.append(Calculate_Internal_Mesh_Settings(fout, writeFlag, mach[i_mach], alt[i_alt], interfaces, sortedInterfaceKeys, min_y_inner))

        Write_BOI_Sizing(fout, body_size_min, bois)
        fout.write('## END MESH SETTINGS FOR MACH ' +str(mach[i_mach]) + ' ALTITUDE ' + str(alt[i_alt]) + ' (ft) .....................\n\n')

def Write_BOI_Sizing(fout, body_size_min, bois):
  bodyFound = False
  for i_boi in bois:
    if 'body' in i_boi:
      body = i_boi
      bodyFound = True
      break

  if bodyFound == True:
    if not body_size_min == False:
      fout.write('Surface:\t' + str(body) + '\n')
      fout.write('\tBody of Influence Size:\t' + str(2*body_size_min) + 'm (' + str(2*body_size_min*1000) +' mm) \n\n')
    else:
      fout.write('NOTE: A BOI with \'body\' in the name was found, but no boundary name with \'body\' was found...\n\n')


def Calculate_Internal_Mesh_Settings(fout, writeFlag, mach, alt, interfaces, sortedInterfaceKeys, min_y_inner):
  global mesh_settings
  global flowthrough
  global fidelity
  global nPrismCells

  y_inner_m = 0
  y_inner_prop = 0

  ## determine if there is an annulus present
  isAnnulus = False
  annulusStartIdx = -1
  annulusEndIdx = -1
  minAnnulusSize = float('inf')
  l_ref = 0
  for i_int in range(len(sortedInterfaceKeys)):
    ## Get the index in interfaces that corresponds to the current sorted interface
    for i_orig in range(len(interfaces)):
      if interfaces[i_orig]['interface'] == sortedInterfaceKeys[i_int]:
        currIdx = i_orig
        break

    if interfaces[currIdx]['diameter'] > interfaces[currIdx]['inner_diameter']:
      isAnnulus = True

      if annulusStartIdx == -1:
        annulusStartIdx = i_int

      if interfaces[currIdx]['inner_diameter'] < minAnnulusSize:
        minAnnulusSize = interfaces[currIdx]['inner_diameter']
        l_ref = interfaces[currIdx]['diameter']*np.pi

    else:
      if isAnnulus == True:
        if annulusEndIdx == -1:
          annulusEndIdx = i_int

  if isAnnulus == True and annulusEndIdx == -1: ## if the whole thing is an annulus
    annulusEndIdx = i_int

  if isAnnulus == True:

    refineAnnulus = False
    if not mesh_settings['internal_num_cell_per_diam'] == 0:
      refineAnnulus = True

    if writeFlag == True:
      fout.write('Surface:\tAnnulus\n')
      fout.write('\tAn annulus was detected between ' + str(sortedInterfaceKeys[annulusStartIdx]) + ' and ' + str(sortedInterfaceKeys[annulusEndIdx]) +' ...\n')
      if refineAnnulus == False:
        fout.write('\tAs per the mesh sizing guidelines, With the fidelity set to ' + str(fidelity) + ' and flowthrough = ' + str(flowthrough) + ' , the annulus does NOT need any local sizings.\n')

    if refineAnnulus == True:
      annulus_s = minAnnulusSize / mesh_settings['internal_num_cell_per_annulus']
      if writeFlag == True:
        fout.write('\tsurface mesh size: ' + str(annulus_s) + ' m  (' + str(annulus_s*1000) + ' mm)\n\n')

    ## calculate the air speed at the given mach/alt combo
    airspeed = Airspeed(alt,mach,airspeed_type='Mach').get_results()

    ## calculate the y_inner_m
    Re_x = airspeed['ReynoldsNumber']*l_ref/10
    Cf = (2*np.log10(Re_x)-0.65)**-2.3
    tau_wall_Pa = Cf*airspeed['dynamicPressure_Pa']
    ustar_mps = np.sqrt(tau_wall_Pa/airspeed['density'])
    y_plus = mesh_settings['y_plus']
    y_inner_m = round((y_plus*airspeed['dynamicViscosity_kgpms']/(airspeed['density']*ustar_mps)), 8)  # y+ = 1 for 10% chord resolution
    if refineAnnulus == True:
      if writeFlag == True:
        fout.write('\tTo achieve a y+ of ' +str(y_plus)+':\tFirst cell height:\t'+str(y_inner_m)+' m ('+str(y_inner_m*1000)+ ' mm)\n')
    if y_inner_m < 1e-6:
      y_inner_m = 1e-6
      if refineAnnulus == True:
        if writeFlag == True:
          fout.write('\tRECOMMENDED first cell height (based on company standards):  1e-6m (.001 mm)\n')
          fout.write('\tNOTE: The number of prism layers is calculated based on y+ = 1e-6m ...\n')

    ## calculate max number of prism layers
    bl_height = y_inner_m
    max_bl_height = minAnnulusSize/2
    y_inner = y_inner_m
    numLayers = 1
    growth_rate = mesh_settings['growth_rate']
    while bl_height < max_bl_height:
      y_inner *= growth_rate
      bl_height += y_inner
      numLayers += 1

    if writeFlag == True:
      fout.write('\t\tThe MAXIMUM number of prism cells == ' + str(numLayers-1) + '\n')
      fout.write('\t\t(This is the absolute maximum number of layers that will fit in the smallest portion of the annulus.)\n')
      fout.write('\t\t(Anything more than this will result in Fluent not putting any BL cells in the annulus.)\n')

    ## calculate max number of prism cells based on min_y_inner
    if not y_inner_m == min_y_inner:
      if writeFlag == True:
        bl_height = min_y_inner
        y_inner = min_y_inner
        max_bl_height = minAnnulusSize/2
        numLayers = 1
        growth_rate = mesh_settings['growth_rate']
        while bl_height < max_bl_height:
          y_inner *= growth_rate
          bl_height += y_inner
          numLayers += 1

      if writeFlag == True:
        y_plus_new = (min_y_inner*(airspeed['density']*ustar_mps))/airspeed['dynamicViscosity_kgpms']
        fout.write('\n\tUsing the smallest first cell height of ' + str(min_y_inner) + ' m (' + str(min_y_inner*1000) + ' mm) results in a y+ of ' + str(round(y_plus_new, 4)) +'\n')
        fout.write('\t\tThe MAXIMUM number of prism cells == ' + str(numLayers-1) + '\n')
        fout.write('\t\t(This is the absolute maximum number of layers that will fit in the smallest portion of the annulus.)\n')
        fout.write('\t\t(Anything more than this will result in Fluent not putting any BL cells in the annulus.)\n\n')
    else:
      if writeFlag == True:
        fout.write('\n\tThis surface has the smallest first cell height\n\n')

    ## INTERNAL SIZING ##
    if not annulusEndIdx == len(sortedInterfaceKeys)-1:
      for i_orig in range(len(interfaces)):
        if interfaces[i_orig]['interface'] == sortedInterfaceKeys[annulusEndIdx]:
          l_ref = interfaces[i_orig]['diameter']
          if refineAnnulus == True:
            prop_s = l_ref / mesh_settings['internal_num_cell_per_diam']
            y_inner_prop = Calculate_Mesh_Settings(fout, writeFlag, mach, alt, ['All surfaces from ' + str(sortedInterfaceKeys[annulusEndIdx]) + ' to ' + str(sortedInterfaceKeys[-1]), l_ref], ['', prop_s], min_y_inner)
          else:
            if writeFlag == True:
              fout.write('Surface:\tPropulsion flowpath\n')
              fout.write('\tAs per the mesh sizing guidelines, With the fidelity set to ' + str(fidelity) + ' and flowthrough = ' + str(flowthrough) + ' , the propulsion surfaces do NOT need any local sizings.\n')
          break

  ##If no annulus is detected
  else:
    refinePropSurfs = False
    if not mesh_settings['internal_num_cell_per_diam'] == 0:
      refinePropSurfs = True

    for i_orig in range(len(interfaces)):
      if interfaces[i_orig]['interface'].lower() == 'aip':
        l_ref = interfaces[i_orig]['diameter']
        if refinePropSurfs == True:
          prop_s = l_ref / mesh_settings['internal_num_cell_per_diam']
          y_inner_prop = Calculate_Mesh_Settings(fout, mach, alt, ['Propulsion surfaces', l_ref], ['', prop_s], min_y_inner)
        else:
          if writeFlag == True:
            fout.write('Surface:\tPropulsion flowpath\n')
            fout.write('\tAs per the mesh sizing guidelines, With the fidelity set to ' + str(fidelity) + ' and flowthrough = ' + str(flowthrough) + ' , the propulsion surfaces do NOT need any local sizings.\n')
        break

  if y_inner_m == 0 and y_inner_prop == 0:
    return(-1)
  elif y_inner_prop == 0:
    return(y_inner_m)
  elif y_inner_m == 0:
    return(y_inner_prop)
  else:
    return(min(y_inner_m, y_inner_prop))

def Calculate_Min_Y_Inner(meshFile, mach, alt, x_surface, s_surface, interfaces, sortedInterfaceKeys):
  min_yplus_values = []
  writeFlag = False
  for i_surf in range(len(x_surface)):
    min_yplus_values.append(Calculate_Mesh_Settings(fout, writeFlag, mach, alt, x_surface[i_surf], s_surface[i_surf], -1))
  if interfaceFound == True:
    min_yplus_values.append(Calculate_Internal_Mesh_Settings(fout, writeFlag, mach, alt, interfaces, sortedInterfaceKeys, -1))

  min_yplus = [i_yp for i_yp in min_yplus_values if i_yp > 0]
  min_yplus = min(min_yplus)

  if min_yplus < 1e-6:
    min_yplus = 1e-6

  return(min_yplus)

def Calculate_Num_Prism_Layers(meshFile, writeFlag, s_surf_full, y_inner_m):
  global mesh_settings
  global nPrismCells

  ## calculate the number of prism cells
  no_layers_field = 1
  no_layers_prisms = 1

  transition_ratio = 0.6
  tol = 0.1
  max_layers = nPrismCells

  growth_rate = mesh_settings['growth_rate']

  s_surface_m = s_surf_full[1]

  converged_flag = False
  results = []

  if writeFlag == True:
    meshFile.write('\tTarget number of prism cells == ' + str(nPrismCells) + '\n')

  for i_field in range(100):
    no_layers_prisms = 1
    for i_prism in range(max_layers):
      y_outer_prisms_m = y_inner_m*growth_rate**(no_layers_prisms-1)
      y_outer_field_m = s_surface_m*growth_rate**(no_layers_field-1)
      y_outer_target_prisms_m = y_outer_field_m*transition_ratio #get close to 1.2 growth rate by matching y_total_field/prisms and y_outer_prism/prism_target
      y_total_prisms_m = y_inner_m*(1-growth_rate**no_layers_prisms)/(1-growth_rate)
      y_total_field_m = s_surface_m*(1-growth_rate**no_layers_field)/(1-growth_rate)

      if y_total_prisms_m > y_total_field_m:
        no_layers_field+=1
        no_layers_prisms=1
        break

      if abs(y_total_prisms_m-y_total_field_m)/y_total_field_m < tol and abs(y_outer_prisms_m-y_outer_target_prisms_m)/y_outer_target_prisms_m < tol:
        converged_flag = True
        results.append([no_layers_prisms,no_layers_field,abs(y_total_prisms_m-y_total_field_m)/y_total_field_m,abs(y_outer_prisms_m-y_outer_target_prisms_m)/y_outer_target_prisms_m])
        if writeFlag == True:
          meshFile.write('\t\tSolution found, no. prisms = ' + str(no_layers_prisms) + '\n')
        break
      else:
        total_error = abs(y_total_prisms_m-y_total_field_m)/y_total_field_m
        transition_error = abs(y_outer_prisms_m-y_outer_target_prisms_m)/y_outer_target_prisms_m
        results.append([no_layers_prisms,no_layers_field,total_error,transition_error,np.sqrt(total_error**2+transition_error**2)])

        no_layers_prisms+=1
      if converged_flag:
        break

  if converged_flag == False:
    results = np.array(results)
    i_results = np.argmin(np.array(results[:,4]))
    if writeFlag == True:
      meshFile.write('\t\tNo solution within tolerance (0.1)\n')
      meshFile.write('\t\tBest you can do is:\n')
      meshFile.write('\t\t\tnumber of prism layers:\t' + str(results[i_results][0]) + '\n')
      meshFile.write('\t\t\tnumber of layers field:\t' + str(results[i_results][1]) + '\n')
      meshFile.write('\t\t\ttotal error:\t' + str(results[i_results][2]) + '\n')
      meshFile.write('\t\t\ttransition error:\t' + str(results[i_results][3]) + '\n')
      meshFile.write('\t\t\tsqrt(total_error**2 + transition_error**2):\t' + str(results[i_results][4]) + '\n')


def Calculate_Mesh_Settings(meshFile, writeFlag, mach, alt, l_ref_full, s_surf_full, min_y_inner):
  global mesh_settings
  global nPrismCells
  y_plus = mesh_settings['y_plus']

  if writeFlag == True:
    meshFile.write('Surface:\t' + l_ref_full[0] + '\n')
  l_ref = l_ref_full[1]

  if writeFlag == True:
    meshFile.write('\tsurface mesh size: ' + str(s_surf_full[1]) +' m  (' + str(s_surf_full[1]*1000) + ' mm)\n\n')

  ## calculate the air speed at the given mach/alt combo
  airspeed = Airspeed(alt,mach,airspeed_type='Mach').get_results()

  ## calculate the y_inner_m
  Re_x = airspeed['ReynoldsNumber']*l_ref/10 # y+ = 1 for 10% chord resolution
  Cf = (2*np.log10(Re_x)-0.65)**-2.3
  tau_wall_Pa = Cf*airspeed['dynamicPressure_Pa']
  ustar_mps = np.sqrt(tau_wall_Pa/airspeed['density'])
  y_inner_m = round((y_plus*airspeed['dynamicViscosity_kgpms']/(airspeed['density']*ustar_mps)), 8)
  if writeFlag == True:
    meshFile.write('\tTo achieve a y+ of ' +str(y_plus)+':\tFirst cell height:\t'+str(y_inner_m)+' m ('+str(y_inner_m*1000)+ 'mm)\n')
  if y_inner_m < 1e-6:
    y_inner_m = 1e-6
    if writeFlag == True:
      meshFile.write('\tRECOMMENDED first cell height (based on company standards):  1e-6m (.001 mm)\n')
      meshFile.write('\tNOTE: The number of prism layers is calculated based on y+ = 1e-6m ...\n')


  Calculate_Num_Prism_Layers(meshFile, writeFlag, s_surf_full, y_inner_m)

  if not y_inner_m == min_y_inner:
    y_plus_new = (min_y_inner*(airspeed['density']*ustar_mps))/airspeed['dynamicViscosity_kgpms']
    if writeFlag == True:
      meshFile.write('\n\tUsing the smallest first cell height of ' + str(min_y_inner) + ' m (' + str(min_y_inner*1000) + ' mm) results in a y+ of ' + str(round(y_plus_new, 4)) +'\n')
      Calculate_Num_Prism_Layers(meshFile, writeFlag, s_surf_full, min_y_inner)
      meshFile.write('\n')
  else:
    if writeFlag == True:
      meshFile.write('\n\tThis surface has the smallest first cell height.\n\n')

  return(y_inner_m)

def Create_Output_File(sc_out_dir):

  outFile = os.path.join(sc_out_dir, 'MESH_VALUES.txt')
  if os.path.exists(outFile):
    os.remove(outFile)

  meshFile = open(outFile, 'w')
  return(meshFile)

def Close_Output_File(meshFile):
  meshFile.close()

## main
if __name__ == '__main__':

  currDir = os.getcwd()

  ## read command line argument and set the appropriate mesh parameters
  case, case_path, CAD_dir = Read_Inputs(currDir)

  ## create the python file to be passed to spaceclaim in batch
  SC_File = Create_SpaceClaim_File(currDir, case_path, CAD_dir)

  ## run SpaceClaim in batch
  meshFile, geomFile, sc_out_dir = Run_Space_Claim(currDir, SC_File)

  ## collect outputs from SpaceClaim
  bois, interfaces, geom_vars, interfaceFound, sortedInterfaceKeys = Collect_SpaceClaim_Data(meshFile, geomFile)

  ## Create the Mesh values file
  fout = Create_Output_File(sc_out_dir)

  ## calculate mesh parameters
  Calculate_Mesh_Parameters(fout, bois, interfaces, geom_vars, interfaceFound, sortedInterfaceKeys)

  ## Close the Mesh values file
  Close_Output_File(fout)
