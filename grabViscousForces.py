# -*- coding: utf-8 -*-
"""
@author: pmills
"""

import os
import sys
import math


def Find_Files(currDir):

  xFile = ''
  yFile = ''
  zFile = ''

  xFound = False
  yFound = False
  zFound = False

  for file in os.listdir(currDir):
    temp = os.path.join(currDir, file)
    if os.path.isfile(temp):
      if 'viscous_vs_inviscid' in file:
        if 'x_component' in file:
          xFile = temp
          xFound = True
        elif 'y_component' in file:
          yFile = temp
          yFound = True
        elif 'z_component' in file:
          zFile = temp
          zFound = True

  if xFound == False:
    print('\nThe x-component file of the forces breakdown was not found...\n')
  if yFound == False:
    print('\nThe y-component file of the forces breakdown was not found...\n')
  if zFound == False:
    print('\nThe z-component file of the forces breakdown was not found...\n')
  if xFound == False and yFound == False and zFound == False:
    print('\n!!ERROR!! None of the viscous_vs_inviscid files were found. Exiting now...\n')
    sys.exit(0)


  return(xFile, yFile, zFile)


def Read_Files(xFile, yFile, zFile):

  data_x = []
  with open(xFile, 'r') as fin:
    for line in fin:
      data_x.append(line.strip().split())
  fin.close()
  if len(data_x) == 0:
    print('\n!!ERROR!! the x component of the forces is empty... Exiting now.\n')
    sys.exit(0)

  ## collect the data
  forces_init = {}
  for line in range(len(data_x)):
    try:
      if 'wall' in data_x[line][0]:
        bc = data_x[line][0]
        forces_init[bc] = {}
        forces_init[bc]['inviscid'] = []
        forces_init[bc]['inviscid'].append(round(float(data_x[line][1].replace('(','')), 8))
        forces_init[bc]['inviscid'].append(round(float(data_x[line][2]),8))
        forces_init[bc]['inviscid'].append(round(float(data_x[line][3].replace(')','')), 8))

        forces_init[bc]['viscous'] = []
        forces_init[bc]['viscous'].append(round(float(data_x[line][4].replace('(','')), 8))
        forces_init[bc]['viscous'].append(round(float(data_x[line][5]), 8))
        forces_init[bc]['viscous'].append(round(float(data_x[line][6].replace(')','')), 8))

      elif '-----' in data_x[line][0]:
        break

    except IndexError:
      continue

  ## combine bcs that are split across an enclosure
  bc_names = list(forces_init.keys())
  forces = {}
  for bc in bc_names:
    if '.' in bc:
      bc_orig = bc.split('.')[0]
      forces[bc_orig]['inviscid'][0] = round(forces[bc_orig]['inviscid'][0] + forces_init[bc]['inviscid'][0], 8)
      forces[bc_orig]['inviscid'][1] = round(forces[bc_orig]['inviscid'][1] + forces_init[bc]['inviscid'][1], 8)
      forces[bc_orig]['inviscid'][2] = round(forces[bc_orig]['inviscid'][2] + forces_init[bc]['inviscid'][2], 8)

      forces[bc_orig]['viscous'][0] = round(forces[bc_orig]['viscous'][0] + forces_init[bc]['viscous'][0], 8)
      forces[bc_orig]['viscous'][1] = round(forces[bc_orig]['viscous'][1] + forces_init[bc]['viscous'][1], 8)
      forces[bc_orig]['viscous'][2] = round(forces[bc_orig]['viscous'][2] + forces_init[bc]['viscous'][2], 8)
    else:
      forces[bc] = {}
      forces[bc] = forces_init[bc]

  ## compute other quantities
  for bc in list(forces.keys()):
    forces[bc]['total_inv'] = {}
    forces[bc]['total_vis'] = {}
    forces[bc]['total'] = {}
    forces[bc]['per_inv'] = {}
    forces[bc]['per_vis'] = {}
    inv = round(math.sqrt(forces[bc]['inviscid'][0]**2 + forces[bc]['inviscid'][1]**2 + forces[bc]['inviscid'][2]**2), 8)
    vis = round(math.sqrt(forces[bc]['viscous'][0]**2 + forces[bc]['viscous'][1]**2 + forces[bc]['viscous'][2]**2), 8)
    tot = round(inv + vis, 8)
    forces[bc]['total_inv'] = inv
    forces[bc]['total_vis'] = vis
    forces[bc]['total'] = tot
    forces[bc]['per_inv'] = round((inv/tot)*100, 8)
    forces[bc]['per_vis'] = round((vis/tot)*100, 8)

  return(forces)

def Write_Data(currDir, forces):


  ## formatting for the output
  max_bc = len(max(list(forces.keys()), key=len)) + 2

  fName = os.path.join(currDir, 'viscous_vs_inviscid forces.txt')
  if os.path.exists(fName):
    os.remove(fName)

  headers = ['Boundary'.ljust(max_bc), 'Total Force(N)  ', 'Total Inviscid  ', 'Total Viscous  ', 'Inviscid(x)   ', 'Inviscid(y)   ', 'Inviscid(z)   ', 'Viscous(x)   ', 'Viscous(y)   ', 'Viscous(z)   ', '% Inviscid  ', '% Viscous']

  with open(fName, 'w') as f:
    for i_header in headers:
      f.write(i_header)
    f.write('\n')
    f.write('-'*175+'\n')
    for bc in list(forces.keys()):
      f.write(str(bc).ljust(len(headers[0])) + str(forces[bc]['total']).ljust(len(headers[1])) +str(forces[bc]['total_inv']).ljust(len(headers[2])) + str(forces[bc]['total_vis']).ljust(len(headers[3])) + str(forces[bc]['inviscid'][0]).ljust(len(headers[4])) + str(forces[bc]['inviscid'][1]).ljust(len(headers[5])) + str(forces[bc]['inviscid'][2]).ljust(len(headers[6])) + str(forces[bc]['viscous'][0]).ljust(len(headers[7])) + str(forces[bc]['viscous'][1]).ljust(len(headers[8])) + str(forces[bc]['viscous'][2]).ljust(len(headers[9])) + str(forces[bc]['per_inv']).ljust(len(headers[10])) + str(forces[bc]['per_vis']).ljust(len(headers[11])) +'\n\n')

  f.close()

if __name__ == '__main__':
  currDir = os.getcwd()

  ## find the files
  xFile, yFile, zFile = Find_Files(currDir)

  ## Read the files
  forces = Read_Files(xFile, yFile, zFile)

  ## print data to file
  Write_Data(currDir, forces)