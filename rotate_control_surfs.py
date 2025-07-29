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


## GLOBAL VARS ##
global surf
global deg
global individual_surf
global deflection


parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument("-case", help="Name of the SpaceClaim file you want to read. If not specified, will default to most recently created SpaceClaim file in the current working directory.", type=str, required=False)
parser.add_argument("-surf", help="Identifying string to tell code which control surface you want to rotate.\nONLY valid option are fin or canard.\nRight now, can only enter one at a time.", type=str, required=True)
parser.add_argument("-deg", help="DEGREES you want to rotate the control surface.\nCan be positive or negative.\nCan enter multiple degrees: -deg 5 10 15 -5", nargs = "+", type = float, required = True)
parser.add_argument("-deflection_type", help="Elevator deflection type you want.\nOptions are: roll, pitch, yaw,\nIf not entered, will default to pitch.", type=str, required=False, default='pitch')
parser.add_argument("-individual_surf", help="Flag to identify if you only want to rotate a specific surface, and not all of them.\nEnter the INTEGER number of the specific surface you want to rotate.\nCan enter multiple numbers: -individual_surf 2 3\nIf this flag is not set, default will be all \'surfs\'.NOTE: If you enter a number that is not is the geometry, it will simply be skipped.", nargs = "+", type=int, required=False)


def Read_Inputs(currDir):
  global surf
  global deg
  global individual_surf
  global deflection

  args = parser.parse_args()
  case = args.case
  surf = args.surf
  deg = args.deg
  individual_surf = args.individual_surf
  deflection = args.deflection_type

  if not (surf == 'fin' or surf == 'canard'):
    print('!! ERROR !! Invalid surf type entered, ' + str(surf) + '. Only allowable surfs are fin or canard. Exiting now.')
    sys.exit(0)

  if not (deflection == 'roll' or deflection == 'pitch' or deflection == 'yaw'):
    print('!! ERROR !! Invlaid deflection_type entered, ' + str(deflection) + '. Only allowable deflection types are roll, pitch, or yaw. Exiting now')
    sys.exit(0)

  ## SpaceClaims roatation scheme is opposite ours.
  for i_deg in range(len(deg)):
    deg[i_deg] *= -1

  CADdir = ''
  if case == None:
    print('>> The case flag was not set, searching the current directory for the most recent .SCDOC file...')
    sc_files = []
    for file in os.listdir(currDir):
      temp = os.path.join(currDir, file)
      if os.path.isfile(temp) and os.path.splitext(temp)[1] == '.scdoc':
        if not ('ROTATE' in temp or 'DUMMY' in temp):
          sc_files.append(file)
          CADdir = currDir

    if len(sc_files) == 0:
      for dir in os.listdir(currDir):
        tempDir = os.path.join(currDir, dir)
        if os.path.isdir(tempDir):
          if str(dir) == 'CAD':
            for file in os.listdir(tempDir):
              tempFile = os.path.join(tempDir, file)
              if os.path.isfile(tempFile) and os.path.splitext(tempFile)[1] == '.scdoc':
                if not ('ROTATE' in tempFile or 'DUMMY' in tempFile):
                  CADdir = tempDir
                  sc_files.append(os.path.join(tempDir, file))
    if len(sc_files) == 0:
      print('>>>> ERROR... no .SCDOC files found in the currrent directory(\''+str(currDir)+'\'). Exiting now')
      sys.exit(0)
    else:
      case = max(sc_files, key=os.path.getmtime)
      case_path = os.path.join(currDir, case)
      print('>>>> The .SCDOC file to be used is:\t'+str(case))
      case = case.replace('.scdoc', '')

  return(case, case_path, CADdir)

def Create_Rotation_File(currDir, case_path):
  global surf
  global deg
  global individual_surf

  SC_FILE = os.path.join(currDir, 'Rotate_Surfs.py')
  if os.path.exists(SC_FILE):
    os.remove(SC_FILE)

  ## spaceclaim requires files to have double backslashes and be enclosed by ""
  case_path = '\"' + case_path.replace('\\','\\\\') + '\"'

  with open(SC_FILE, 'w') as f:
    f.write('# Python Script, API Version = V242\n\n')
    f.write('from SpaceClaim.Api.V242 import *\n')
    f.write('from System import Array\n')
    f.write('from SpaceClaim.Api.V242.Geometry import IShape\n')
    f.write('from SpaceClaim.Api.V242.Geometry import Matrix, Frame\n')
    f.write('from SpaceClaim.Api.V242.Geometry import PointUV\n')
    f.write('import math\n\n')

    f.write('DocumentOpen.Execute(' + case_path + ')\n\n')
    f.write('doc = Window.ActiveWindow.Document\n')
    f.write('main_part = doc.MainPart\n\n')
    f.write('identity_matrix = Matrix.Identity\n\n')
    f.write('global SURFS\n')
    f.write('global individual_surf\n')
    f.write('global surf\n')
    f.write('global degrees\n')
    f.write('global deflection\n')
    f.write('SURFS = {}\n')
    f.write('surf = \'' + str(surf) + '\'\n')
    if individual_surf == None:
      f.write('individual_surf = None\n\n')
    else:
      f.write('individual_surf = [')
      if len(individual_surf) == 1:
        f.write(str(individual_surf[0]) + ']\n')
      else:
        for i_surf in range(len(individual_surf)-1):
          f.write(str(individual_surf[i_surf]) + ',')
        f.write(str(individual_surf[-1]) + ']\n')

    f.write('degrees = [')
    if len(deg) == 1:
      f.write(str(float(deg[0])) + ']\n')
    else:
      for i_deg in range(len(deg)-1):
        f.write(str(float(deg[i_deg])) + ',')
      f.write(str(float(deg[-1])) + ']\n')

    f.write('deflection = \'' + str(deflection) + '\'\n\n')

    f.write('names_to_exclude = [\'farfield\', \'interface\', \'outlet\', \'symmetry\']\n')
    f.write('goodNames = []\n')
    f.write('badNames = []\n')
    f.write('badAreas = []\n')
    f.write('BCs = NamedSelection.GetGroups()\n')
    f.write('for name in BCs:\n')
    f.write('\tbadNameFound = False\n')
    f.write('\tfor i_name in names_to_exclude:\n')
    f.write('\t\tif i_name in name.GetName().ToString().ToLower():\n')
    f.write('\t\t\tbadNameFound = True\n')
    f.write('\t\t\tbadNames.append(name.GetName().ToString())\n')
    f.write('\t\t\tbadFaces = FaceSelection.CreateByGroups(name.GetName().ToString()).Items\n')
    f.write('\t\t\tfor i_bad_face in badFaces:\n')
    f.write('\t\t\t\tbadAreas.append(round(i_bad_face.Area, 8))\n')
    f.write('\tif badNameFound == False:\n')
    f.write('\t\tgoodNames.append(name.GetName().ToString())\n\n')
    f.write('original_bodies = []\n')
    f.write('for i_body in main_part.Bodies:\n')
    f.write('\toriginal_bodies.append(i_body.Name)\n\n')
    f.write('original_enclosures = []\n')
    f.write('for part in doc.Parts:\n')
    f.write('\tfor body in part.Bodies:\n')
    f.write('\t\tif \'enclosure\' in part.Name.ToLower() or \'enclosure\' in body.Name.ToLower():\n')
    f.write('\t\t\toriginal_enclosures.append(body.Name)\n\n')

    f.write('bodyNames = []\n')
    f.write('bodyCheckNames = [\'body\', \'boattail\', \'boat_tail\']\n')
    f.write('for i_goodName in range(len(goodNames)):\n')
    f.write('\tif any(i_bodyName in goodNames[i_goodName] for i_bodyName in bodyCheckNames):\n')
    f.write('\t\tbodyNames.append(goodNames[i_goodName])\n')
    f.write('for i_bodyName in bodyNames:\n')
    f.write('\tgoodNames.remove(i_bodyName)\n')
    f.write('\tgoodNames.append(i_bodyName)\n\n')

    f.write('created_enclosures = []\n')
    f.write('created_physical_enclosures = []\n')
    f.write('for i_name in range(len(goodNames)):\n')
    f.write('\tselFaces = FaceSelection.CreateByGroups(goodNames[i_name])\n')
    f.write('\tresult = DetachFaces.Execute(selFaces)\n\n')
    f.write('\tnew_enclosures = []\n')
    f.write('\tnew_enclosures_names = []\n')
    f.write('\tfor part in doc.Parts:\n')
    f.write('\t\tfor body in part.Bodies:\n')
    f.write('\t\t\tnew_enclosures.append(body)\n')
    f.write('\t\t\tnew_enclosures_names.append(body.Name)\n\n')
    f.write('\tcount = 2\n')
    f.write('\tfor i_enclosure in new_enclosures:\n')
    f.write('\t\tif not i_enclosure.Name in original_bodies:\n')
    f.write('\t\t\tif not i_enclosure.Name in original_enclosures:\n')
    f.write('\t\t\t\tselection = Selection.Create(i_enclosure)\n')
    f.write('\t\t\t\tif \'CREATED_\' + str(goodNames[i_name]) in original_enclosures:\n')
    f.write('\t\t\t\t\tresult = RenameObject.Execute(selection, \'CREATED_\' + str(count) + \'_\' + goodNames[i_name])\n')
    f.write('\t\t\t\t\tcount += 1\n')
    f.write('\t\t\t\telse:\n')
    f.write('\t\t\t\t\tresult = RenameObject.Execute(selection, \'CREATED_\' + goodNames[i_name])\n\n')
    f.write('\t\t\t\tcurr_enclosure_faces = i_enclosure.Shape.Faces\n')
    f.write('\t\t\t\tbadFace_idx = []\n')
    f.write('\t\t\t\tfor i_face in curr_enclosure_faces:\n')
    f.write('\t\t\t\t\tif round(i_face.Area,8) in badAreas:\n')
    f.write('\t\t\t\t\t\tbadFace_idx.append(i_face)\n\n')
    f.write('\t\t\t\tisUpdateEnclosures = True\n')
    f.write('\t\t\t\tif not len(badFace_idx) == 0:\n')
    f.write('\t\t\t\t\tif len(badFace_idx) == len(curr_enclosure_faces):\n')
    f.write('\t\t\t\t\t\tresult = Delete.Execute(selection)\n')
    f.write('\t\t\t\t\t\tisUpdateEnclosures = False\n\n')
    f.write('\t\t\t\tif isUpdateEnclosures == True:\n')
    f.write('\t\t\t\t\tif count == 2:\n')
    f.write('\t\t\t\t\t\toriginal_enclosures.append(\'CREATED_\' + goodNames[i_name])\n')
    f.write('\t\t\t\t\t\tcreated_enclosures.append(\'CREATED_\' + goodNames[i_name])\n')
    f.write('\t\t\t\t\telse:\n')
    f.write('\t\t\t\t\t\toriginal_enclosures.append(\'CREATED_\' + str(count-1) + \'_\' + goodNames[i_name])\n')
    f.write('\t\t\t\t\t\tcreated_enclosures.append(\'CREATED_\' + str(count-1) + \'_\' + goodNames[i_name])\n')
    f.write('\t\t\t\t\tcreated_physical_enclosures.append(i_enclosure)\n\n')

    f.write('enclosures_to_delete = set()\n')
    f.write('final_enclosures = []\n')
    f.write('for i_enc in created_physical_enclosures:\n')
    f.write('\tarea = 0\n')
    f.write('\tfor i_face in i_enc.Shape.Faces:\n')
    f.write('\t\tarea += i_face.Area\n')
    f.write('\tarea = round(area,8)\n\n')
    f.write('\tbadAreaFound = False\n')
    f.write('\tfor i_area in badAreas:\n')
    f.write('\t\tif abs(i_area - area) <= 0.0000001:\n')
    f.write('\t\t\tbadAreaFound = True\n')
    f.write('\t\t\tenclosures_to_delete.add(i_enc)\n')
    f.write('\t\t\tbreak\n')
    f.write('\tif badAreaFound == False:\n')
    f.write('\t\tfinal_enclosures.append(i_enc)\n\n')
    f.write('for i in enclosures_to_delete:\n')
    f.write('\tselection = Selection.Create(i)\n')
    f.write('\tresult = Delete.Execute(selection)\n\n')

    f.write('def Calculate_Midpoint(pt1, pt2):\n')
    f.write('\tmidpoint = [((pt1[0] + pt2[0])/2), ((pt1[1] + pt2[1])/2), ((pt1[2] + pt2[2])/2)]\n')
    f.write('\treturn(midpoint)\n\n')

    f.write('def Calculate_Control_Surface_MRP(enclosure):\n')
    f.write('\tmax_area = float(\'-inf\')\n')
    f.write('\tmax_area_face = \'\'\n')
    f.write('\tfor i_face in i_enclosure.Faces:\n')
    f.write('\t\tarea = i_face.Shape.Area\n')
    f.write('\t\tif area > max_area:\n')
    f.write('\t\t\tmax_area = area\n')
    f.write('\t\t\tmax_area_face = i_face\n\n')
    f.write('\twidth_faces = []\n')
    f.write('\tfor i_face in i_enclosure.Faces:\n')
    f.write('\t\tarea = i_face.Shape.Area\n')
    f.write('\t\tif abs(area - max_area) <= 0.0001:\n')
    f.write('\t\t\twidth_faces.append(i_face)\n\n')
    f.write('\tadjacentFaces = set()\n')
    f.write('\tmax_area_face = width_faces[0]\n')
    f.write('\tbox = max_area_face.Shape.GetBoundingBox(identity_matrix)\n')
    f.write('\tadjacentFaces.add(tuple([max_area_face.Shape, box.MinCorner.X]))\n')
    f.write('\tfor i_edge in max_area_face.Shape.Edges:\n')
    f.write('\t\tfor i_adjFace in i_edge.Faces:\n')
    f.write('\t\t\tif i_adjFace is not max_area_face:\n')
    f.write('\t\t\t\tif isinstance(i_adjFace.Geometry, Plane):\n')
    f.write('\t\t\t\t\tbox = i_adjFace.GetBoundingBox(identity_matrix)\n')
    f.write('\t\t\t\t\tadjacentFaces.add(tuple([i_adjFace, box.MinCorner.X]))\n\n')
    f.write('\tadjacentFaces = [list(x) for x in adjacentFaces]\n')
    f.write('\tadjacentFaces = sorted(adjacentFaces, key = lambda x:(x[1]))\n')
    f.write('\tfacePoints_min_A = set()\n')
    f.write('\tfor i_loop in adjacentFaces[0][0].Loops:\n')
    f.write('\t\tfor i_edge in i_loop.Edges:\n')
    f.write('\t\t\tpt = i_edge.StartPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_min_A.add(tuple(pt))\n')
    f.write('\t\t\tpt = i_edge.EndPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_min_A.add(tuple(pt))\n')
    f.write('\tfacePoints_min_A = [list(i) for i in facePoints_min_A]\n')
    f.write('\tmin_point_root_A =  sorted(facePoints_min_A, key = lambda x:(x[0]))[0]\n\n')
    f.write('\tif not abs(facePoints_min_A[2][0]) == abs(facePoints_min_A[2][-1]):\n')
    f.write('\t\tmax_y = max(abs(p[1]) for p in facePoints_min_A)\n')
    f.write('\t\ttop_points = [p for p in facePoints_min_A if abs(abs(p[1]) - max_y) <= 1e-2]\n')
    f.write('\t\tmin_point_tip_A = sorted(top_points, key = lambda x:(x[0]))[0]\n')
    f.write('\telse:\n')
    f.write('\t\tmax_y = max(abs(p[2]) for p in facePoints_min_A)\n')
    f.write('\t\ttop_points = [p for p in facePoints_min_A if abs(abs(p[2]) - max_y) <= 1e-2]\n')
    f.write('\t\tmin_point_tip_A = sorted(top_points, key = lambda x:(x[0]))[0]\n')
    f.write('\tfacePoints_max_A = set()\n')
    f.write('\tfor i_loop in adjacentFaces[-1][0].Loops:\n')
    f.write('\t\tfor i_edge in i_loop.Edges:\n')
    f.write('\t\t\tpt = i_edge.StartPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_max_A.add(tuple(pt))\n')
    f.write('\t\t\tpt = i_edge.EndPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_max_A.add(tuple(pt))\n\n')
    f.write('\tfacePoints_max_A = [list(i) for i in facePoints_max_A]\n')
    f.write('\troot_max_x = max(p[0] for p in facePoints_max_A)\n')
    f.write('\tright_points = [p for p in facePoints_max_A if abs(p[0] - root_max_x) <= 1e-2]\n')
    f.write('\tmax_point_root_A = sorted(right_points, key = lambda x:(abs(x[2])))[0]\n')
    f.write('\tif not abs(facePoints_min_A[2][0]) == abs(facePoints_min_A[2][-1]):\n')
    f.write('\t\tmax_y = max(abs(p[1]) for p in facePoints_max_A)\n')
    f.write('\t\ttop_points = [p for p in facePoints_max_A if abs(abs(p[1]) -max_y) <= 1e-2]\n')
    f.write('\t\tmax_point_tip_A = sorted(top_points, key = lambda x:(x[0]))[-1]\n')
    f.write('\telse:\n')
    f.write('\t\tmax_y = max(abs(p[2]) for p in facePoints_max_A)\n')
    f.write('\t\ttop_points = [p for p in facePoints_max_A if abs(abs(p[2]) -max_y) <= 1e-2]\n')
    f.write('\t\tmax_point_tip_A = sorted(top_points, key = lambda x:(x[0]))[-1]\n\n')
    f.write('\tadjacentFaces = set()\n')
    f.write('\tmax_Area_face = width_faces[-1]\n')
    f.write('\tbox = max_Area_face.Shape.GetBoundingBox(identity_matrix)\n')
    f.write('\tadjacentFaces.add(tuple([max_Area_face.Shape, box.MinCorner.X]))\n')
    f.write('\tfor i_edge in max_Area_face.Shape.Edges:\n')
    f.write('\t\tfor i_BdjFace in i_edge.Faces:\n')
    f.write('\t\t\tif i_BdjFace is not max_Area_face:\n')
    f.write('\t\t\t\tif isinstance(i_BdjFace.Geometry, Plane):\n')
    f.write('\t\t\t\t\tbox = i_BdjFace.GetBoundingBox(identity_matrix)\n')
    f.write('\t\t\t\t\tadjacentFaces.add(tuple([i_BdjFace, box.MinCorner.X]))\n\n')
    f.write('\tadjacentFaces = [list(x) for x in adjacentFaces]\n')
    f.write('\tadjacentFaces = sorted(adjacentFaces, key = lambda x:(x[1]))\n')
    f.write('\tfacePoints_min_B = set()\n')
    f.write('\tfor i_loop in adjacentFaces[0][0].Loops:\n')
    f.write('\t\tfor i_edge in i_loop.Edges:\n')
    f.write('\t\t\tpt = i_edge.StartPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_min_B.add(tuple(pt))\n')
    f.write('\t\t\tpt = i_edge.EndPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_min_B.add(tuple(pt))\n')
    f.write('\tfacePoints_min_B = [list(i) for i in facePoints_min_B]\n')
    f.write('\tmin_point_root_B =  sorted(facePoints_min_B, key = lambda x:(x[0]))[0]\n')
    f.write('\tif not abs(facePoints_min_B[2][0]) == abs(facePoints_min_B[2][-1]):\n')
    f.write('\t\tmax_y = max(abs(p[1]) for p in facePoints_min_B)\n')
    f.write('\t\ttop_points = [p for p in facePoints_min_B if abs(abs(p[1]) - max_y) <= 1e-2]\n')
    f.write('\t\tmin_point_tip_B = sorted(top_points, key = lambda x:(x[0]))[0]\n')
    f.write('\telse:\n')
    f.write('\t\tmax_y = max(abs(p[2]) for p in facePoints_min_B)\n')
    f.write('\t\ttop_points = [p for p in facePoints_min_B if abs(abs(p[2]) - max_y) <= 1e-2]\n')
    f.write('\t\tmin_point_tip_B = sorted(top_points, key = lambda x:(x[0]))[0]\n\n')
    f.write('\tfacePoints_max_B = set()\n')
    f.write('\tfor i_loop in adjacentFaces[-1][0].Loops:\n')
    f.write('\t\tfor i_edge in i_loop.Edges:\n')
    f.write('\t\t\tpt = i_edge.StartPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_max_B.add(tuple(pt))\n')
    f.write('\t\t\tpt = i_edge.EndPoint\n')
    f.write('\t\t\tpt = [round(i,8) for i in [pt.X, pt.Y, pt.Z]]\n')
    f.write('\t\t\tfacePoints_max_B.add(tuple(pt))\n\n')
    f.write('\tfacePoints_max_B = [list(i) for i in facePoints_max_B]\n')
    f.write('\troot_max_x = max(p[0] for p in facePoints_max_B)\n')
    f.write('\tright_points = [p for p in facePoints_max_B if abs(p[0] - root_max_x) <= 1e-2]\n')
    f.write('\tmax_point_root_B = sorted(right_points, key = lambda x:(abs(x[2])))[0]\n\n')
    f.write('\tif not abs(facePoints_min_B[2][0]) == abs(facePoints_min_B[2][-1]):\n')
    f.write('\t\tmax_y = max(abs(p[1]) for p in facePoints_max_B)\n')
    f.write('\t\ttop_points = [p for p in facePoints_max_B if abs(abs(p[1]) -max_y) <= 1e-2]\n')
    f.write('\t\tmax_point_tip_B = sorted(top_points, key = lambda x:(x[0]))[-1]\n')
    f.write('\telse:\n')
    f.write('\t\tmax_y = max(abs(p[2]) for p in facePoints_max_B)\n')
    f.write('\t\ttop_points = [p for p in facePoints_max_B if abs(abs(p[2]) -max_y) <= 1e-2]\n')
    f.write('\t\tmax_point_tip_B = sorted(top_points, key = lambda x:(x[0]))[-1]\n\n')
    f.write('\tmid_A = Calculate_Midpoint(min_point_root_A, max_point_root_A)\n')
    f.write('\tmid_B = Calculate_Midpoint(min_point_root_B, max_point_root_B)\n')
    f.write('\tmrp_x = Calculate_Midpoint([mid_A[0],0,0], [mid_B[0],0,0])\n')
    f.write('\tmrp_y = Calculate_Midpoint([0,mid_A[1],0], [0,mid_B[1],0])\n')
    f.write('\tmrp_z = Calculate_Midpoint([0,0,mid_A[2]], [0,0,mid_B[2]])\n')

    f.write('\tmid_A_tip = Calculate_Midpoint(min_point_tip_A, max_point_tip_A)\n')
    f.write('\tmid_B_tip = Calculate_Midpoint(min_point_tip_B, max_point_tip_B)\n')
    f.write('\tmrp_x_tip = Calculate_Midpoint([mid_A_tip[0],0,0], [mid_B_tip[0],0,0])\n')
    f.write('\tmrp_y_tip = Calculate_Midpoint([0,mid_A_tip[1],0], [0,mid_B_tip[1],0])\n')
    f.write('\tmrp_z_tip = Calculate_Midpoint([0,0,mid_A_tip[2]], [0,0,mid_B_tip[2]])\n')
    f.write('\tenclosureName = enclosure.Name.replace(\'CREATED_\', \'\')\n')
    f.write('\tSURFS[enclosureName] = {}\n')
    f.write('\tSURFS[enclosureName][\'mrp_root\'] = [mrp_x[0], mrp_y[1], mrp_z[2]]\n')
    f.write('\tSURFS[enclosureName][\'mrp_tip\'] = [mrp_x_tip[0], mrp_y_tip[1], mrp_z_tip[2]]\n\n')

    f.write('control_surfaces = [\'canard\', \'fin\', \'tailfin\']\n')
    f.write('for i_enclosure in final_enclosures:\n')
    f.write('\ttry:\n')
    f.write('\t\tif len(i_enclosure) > 1:\n')
    f.write('\t\t\tmin_x = min_y = min_z = float(\'inf\')\n')
    f.write('\t\t\tmax_x = max_y = max_z = float(\'-inf\')\n')
    f.write('\t\t\tfor split in i_enclosure:\n')
    f.write('\t\t\t\tbox = split.Shape.GetBoundingBox(identity_matrix)\n')
    f.write('\t\t\t\tmin_x = min(min_x, box.MinCorner.X)\n')
    f.write('\t\t\t\tmin_y = min(min_y, box.MinCorner.Y)\n')
    f.write('\t\t\t\tmin_z = min(min_z, box.MinCorner.Z)\n')
    f.write('\t\t\t\tmax_x = max(max_x, box.MaxCorner.X)\n')
    f.write('\t\t\t\tmax_y = max(max_y, box.MaxCorner.Y)\n')
    f.write('\t\t\t\tmax_z = max(max_z, box.MaxCorner.Z)\n')
    f.write('\t\t\tlength_x = max_x - min_x\n')
    f.write('\t\t\tlength_y = max_y - min_y\n')
    f.write('\t\t\tlength_z = max_z - min_z\n')
    f.write('\texcept TypeError:\n')
    f.write('\t\tif any(i_name in i_enclosure.Name.ToLower() for i_name in control_surfaces):\n')
    f.write('\t\t\tCalculate_Control_Surface_MRP(i_enclosure)\n\n')

    f.write('fName = \'' + str(case_path.replace('.scdoc"', '').replace('\"','')) + '\'\n')
    f.write('saveNameDummy = str(fName) + \'_DUMMY\'+ \'.scdoc\'\n')
    f.write('options = ExportOptions.Create()\n')
    f.write('DocumentSave.Execute(saveNameDummy, options)\n\n')

    f.write('def Determine_Degree(surf, deg):\n')
    f.write('\tglobal deflection\n')
    f.write('\ttry:\n')
    f.write('\t\tsurfNum = str(int(surf[-1]))\n')
    f.write('\texcept ValueError:\n')
    f.write('\t\traise Exception(\'surface: \' + str(surf) + \' does not follow the naming convention (last character should be a number.) Exiting now.\')\n\n')
    f.write('\tif deflection == \'roll\':\n')
    f.write('\t\tdeg = deg\n')
    f.write('\telif deflection == \'pitch\':\n')
    f.write('\t\tif surfNum == \'1\' or surfNum == \'2\':\n')
    f.write('\t\t\tdeg *= -1\n')
    f.write('\t\telse:\n')
    f.write('\t\t\tdeg = deg\n')
    f.write('\telif deflection == \'yaw\':\n')
    f.write('\t\tif surfNum == \'2\' or surfNum == \'3\':\n')
    f.write('\t\t\tdeg *= -1\n')
    f.write('\t\telse:\n')
    f.write('\t\t\tdeg = deg\n\n')
    f.write('\treturn(deg)\n\n')


    f.write('def sign(val, ref):\n')
    f.write('\tif ref > 0:\n')
    f.write('\t\treturn(val)\n')
    f.write('\telse:\n')
    f.write('\t\treturn(val*-1)\n\n')
    ## Create a plane on each on the surfs
    f.write('plane_surfs = []\n')
    f.write('if not individual_surf == None:\n')
    f.write('\tfor i_ind in individual_surf:\n')
    f.write('\t\tfor key in SURFS:\n')
    f.write('\t\t\tif str(surf)+str(i_ind) in key or str(surf)+\'_\'+str(i_ind) in key:\n')
    f.write('\t\t\t\tplane_surfs.append(key)\n')
    f.write('else:\n')
    f.write('\tfor key in SURFS:\n')
    f.write('\t\tif surf in key:\n')
    f.write('\t\t\tplane_surfs.append(key)\n\n')

    f.write('for i_deg in degrees:\n')
    f.write('\tWindow.ActiveWindow.Close()\n')
    f.write('\tDocumentOpen.Execute(' + case_path + ')\n\n')
    f.write('\tfor i_surf in plane_surfs:\n')
    f.write('\t\tif not individual_surf == None:\n')
    f.write('\t\t\tWindow.ActiveWindow.Close()\n')
    f.write('\t\t\tDocumentOpen.Execute(' + case_path + ')\n\n')
    f.write('\t\tselection = Selection.Create(GetRootPart().CoordinateSystems[0].Axes[0])\n')
    f.write('\t\tresult = DatumPlaneCreator.Create(selection, False, None)\n')
    f.write('\t\tselection = Selection.Create(GetRootPart().DatumPlanes[0])\n')
    f.write('\t\tdirection = Direction.DirX\n')
    f.write('\t\toptions = MoveOptions()\n')
    f.write('\t\tx = SURFS[i_surf][\'mrp_root\'][0]*1000\n')
    f.write('\t\ty1 = SURFS[i_surf][\'mrp_root\'][1]*1000\n')
    f.write('\t\tz1 = SURFS[i_surf][\'mrp_root\'][2]*1000\n')
    f.write('\t\ty2 = SURFS[i_surf][\'mrp_tip\'][1]*1000\n')
    f.write('\t\tz2 = SURFS[i_surf][\'mrp_tip\'][2]*1000\n')

    f.write('\t\tresult = Move.Translate(selection, direction, MM(x), options)\n')
    f.write('\t\tsectionPlane = Plane.Create(Frame.Create(Point.Create(MM(x), MM(0), MM(0)),Direction.DirY, Direction.DirZ))\n')
    f.write('\t\tresult = ViewHelper.SetSketchPlane(sectionPlane, None)\n')
    f.write('\t\tstart = Point2D.Create(MM(y1), MM(z1))\n')
    f.write('\t\tend = Point2D.Create(MM(y2), MM(z2))\n')
    f.write('\t\tresult = SketchLine.Create(start, end)\n')
    f.write('\t\tbaseSel = SelectionPoint.Create(GetRootPart().DatumPlanes[0].Curves[0])\n')
    f.write('\t\ttargetSel = SelectionPoint.Create(GetRootPart().DatumPlanes[0])\n')
    f.write('\t\tmode = InteractionMode.Solid\n')
    f.write('\t\tresult = ViewHelper.SetViewMode(mode, None)\n')
    f.write('\t\tselection = Selection.Create(GetRootPart().DatumPlanes[0])\n')
    f.write('\t\tresult = Delete.Execute(selection)\n')

    f.write('\t\tfor i_name in range(len(goodNames)):\n')
    f.write('\t\t\tif goodNames[i_name] == i_surf:\n')
    f.write('\t\t\t\tselFaces = FaceSelection.CreateByGroups(goodNames[i_name])\n')
    f.write('\t\t\t\tanchor = Selection.Create(GetRootPart().Curves[0])\n')
    f.write('\t\t\t\tanchorPoint = Move.GetAnchorPoint(anchor)\n')
    f.write('\t\t\t\tdeg = 0.707106781186548\n')
    f.write('\t\t\t\taxis = Line.Create(anchorPoint, Direction.Create(0, sign(deg, y1), sign(deg, z1)))\n')
    f.write('\t\t\t\toptions = MoveOptions()\n')
    f.write('\t\t\t\ti_surfDeg = Determine_Degree(i_surf, DEG(i_deg))\n')
    f.write('\t\t\t\tresult = Move.Rotate(selFaces, axis, DEG(i_deg), options)\n')
    f.write('\t\tselection = Selection.Create(GetRootPart().Curves[0])\n')
    f.write('\t\tresult = Delete.Execute(selection)\n\n')
    f.write('\t\tif not individual_surf == None:\n')
    f.write('\t\t\tif i_deg < 0:\n')
    f.write('\t\t\t\tfName = \'' + str(case_path.replace('.scdoc"', '').replace('\"','')) + '\'\n')
    f.write('\t\t\t\tsaveName = str(fName) + \'_ROTATE_\' + str(deflection.upper()) + \'_\' + str(i_surf.replace(\'wall_aero\',\'\').upper()) + \'_POS_\' + str(i_deg).replace(\'-\',\'\').replace(\'.\', \'p\') + \'_deg.scdoc\'\n')
    f.write('\t\t\telse:\n')
    f.write('\t\t\t\tfName = \'' + str(case_path.replace('.scdoc"', '').replace('\"','')) + '\'\n')
    f.write('\t\t\t\tsaveName = str(fName) + \'_ROTATE_\' + str(deflection.upper()) + \'_\' + str(i_surf.replace(\'wall_aero\',\'\').upper()) + \'_NEG_\' + str(i_deg).replace(\'.\', \'p\') + \'_deg.scdoc\'\n')
    f.write('\t\t\toptions = ExportOptions.Create()\n')
    f.write('\t\t\tDocumentSave.Execute(saveName, options)\n\n')
    f.write('\tif individual_surf == None:\n')
    f.write('\t\tif i_deg < 0:\n')
    f.write('\t\t\tfName = \'' + str(case_path.replace('.scdoc"', '').replace('\"','')) + '\'\n')
    f.write('\t\t\tsaveName = str(fName) + \'_ROTATE_\' + str(deflection.upper()) + \'_\' + str(surf.upper()) + \'S_POS_\' + str(i_deg).replace(\'-\',\'\').replace(\'.\', \'p\') + \'_deg.scdoc\'\n')
    f.write('\t\telse:\n')
    f.write('\t\t\tfName = \'' + str(case_path.replace('.scdoc"', '').replace('\"','')) + '\'\n')
    f.write('\t\t\tsaveName = str(fName) + \'_ROTATE_\' + str(deflection.upper()) + \'_\' + str(surf.upper()) + \'S_NEG_\' + str(i_deg).replace(\'.\', \'p\') + \'_deg.scdoc\'\n')
    f.write('\t\toptions = ExportOptions.Create()\n')
    f.write('\t\tDocumentSave.Execute(saveName, options)\n\n')

  f.close()

  return(SC_FILE)

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

  for file in os.listdir(currDir):
    tempFile = os.path.join(currDir, file)
    if os.path.isfile(tempFile) and os.path.splitext(tempFile)[1] == '.scdoc':
      if 'DUMMY' in file:
        print('\n>> Removing the \'DUMMY\' file ...\n')
        os.remove(tempFile)



if __name__ == '__main__':

  currDir = os.getcwd()

  ## read command line arguments
  case, case_path, CADdir = Read_Inputs(currDir)

  ## Incase of a nested directory, update the currDir
  currDir = CADdir

  ## Create the python file to be passed to spaceclaim in batch
  SC_FILE = Create_Rotation_File(currDir, case_path)

  ## run SpaceClaim in batch mode
  Run_Space_Claim(currDir, SC_FILE)
