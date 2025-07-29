import os

def calculate_iters_to_ignore(idx, curr_BP, prev_BP):
  iters_to_ignore = 0
  if idx == 0:
    if curr_BP <= 50000:
      iters_to_ignore = 650 # 50@CFL=10, 500@CFL=40, 100@BackPress
    elif curr_BP <= 85000:
      iters_to_ignore = 1050 # 50@CFL=10, 500@CFL=40, 5*100@BackPress-ramping
    else:
      iters_to_ignore = 1550 # 50@CFL=10, 500@CFL=40, 10*100@BackPress-ramping
  else:
    ## reading from restarts
    if abs(curr_BP - prev_BP) <= 2500:
      iters_to_ignore = 101 # 1@restart, 100 @ BackPress
    elif abs(curr_BP - prev_BP) <= 5000:
      iters_to_ignore = 501 # 1@restart, 5*100@BackPress-ramping
    else:
      iters_to_ignore = 1001 # 1@restart, 10*100@BackPress-ramping
  return(iters_to_ignore)


def set_view_for_contours(f, currDir):
  ## create the parent directory for the images
  parentImageDir = os.path.join(currDir, 'IMAGES')
  if not os.path.exists(parentImageDir):
    os.mkdir(parentImageDir)

  ## create the view file for the contours
  viewFile = os.path.join(parentImageDir, 'inlet_2d.vw')
  if not os.path.exists(viewFile):
    with open(viewFile, 'w') as fout:
      fout.write('(38 ((\n')
      fout.write('(view-list (\n')
      fout.write('(inlet_2d ((0.1736933290958405 0.08735045790672302 0.6793386340141296) (0.1736933290958405 0.08735045790672302 1.787916801276879e-07) (-9.317591320723295e-08 0.9999998807907104 4.656612873077393e-10) 0.1921459138393402 0.1921459138393402 "perspective") #(1. 0. 0. 0. 0. 1. 0. 0. 0. 0. 1. 0. 0. 0. 0. 1.))\n')
      fout.write(')))))\n')
    fout.close()

  f.writelines('/preferences/graphics/colormap-settings/number-format-precision 2\n')
  f.writelines('/preferences/graphics/colormap-settings/number-format-type float\n')
  f.writelines('/preferences/graphics/colormap-settings/text-font-automatic-vertical-size 0.125\n')
  f.writelines('/preferences/graphics/colormap-settings/text-font-name \"Helvetica-Bold\"\n')
  f.writelines('/preferences/graphics/lighting/headlight \"On\"\n')
  f.writelines('display/views/read-views \"' + str(viewFile) + '\"\n')


def create_contours(f, currDir, case_name_full):

  ## create the parent directory for the images
  parentImageDir = os.path.join(currDir, 'IMAGES')
  if not os.path.exists(parentImageDir):
    os.mkdir(parentImageDir)


  ## dispaly and save every contour
  contours = ["static_pressure", "mach", "total_pressure", "static_temp"]
  for i_contour in range(len(contours)):
    imageDir = os.path.join(parentImageDir, contours[i_contour])

    if not os.path.exists(imageDir):
      os.mkdir(imageDir)
    imageDir = imageDir.replace('\\','\\\\' )

    f.writelines('display/objects/display ' + str(contours[i_contour]) + '\n')
    f.writelines('display/views/restore-view \"inlet_2d\" quit quit\n')
    saveName = os.path.join(imageDir, case_name_full + '_' + contours[i_contour] + '.png')
    saveName = '\"' + saveName + '\"'
    f.writelines('display/save-picture ' + saveName + '\n')



# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 08:36:24 2025

@author: pmills
"""

