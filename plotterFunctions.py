# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 10:54:29 2025

@author: pmills
"""

import os
import re
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator, FuncFormatter
import statistics


####################################################
## CREATE_RESIDUAL_PLOTS
####################################################
def Create_Residual_Plots(rescale_file_dir):

  ##PCM: THINGS TO ADD
  ## 4) IFF back pressuring is being applied, add when the back pressure was applied to convergecne plots (probably only on the user defined plots.))

  ##FUTURE work, if there are multiple process_output files, parse all of them or only the most recent?

  print('\n>>Creating residual plots...')

  ## create the residual plots directory
  residual_plotDir = os.path.join(rescale_file_dir, 'residual_plots')
  if not os.path.exists(residual_plotDir):
    os.mkdir(residual_plotDir)


  ## parse the working directory and grab the process_output file.
  for file in os.listdir(rescale_file_dir):
    tempFile = os.path.join(rescale_file_dir, file)
    if os.path.isfile(tempFile):
      if 'process_output' in file:
        outFile = tempFile

        break

  ## read the whole file into data
  data = []
  with open(outFile, 'r') as fin:
    for line in fin:
      data.append(line.strip().split())

  fin.close()

  ## check for timestamps in the process_output file
  offsetIdx = 0
  if data[0][0].startswith('[') and data[0][0].endswith(']:'):
    offsetIdx = 1

  ## Grab the user defined convergence metrics
  convMetrics = []
  for i_line in range(len(data)):
    if 'solve/monitors/residual/convergence-criteria' in data[i_line]:
      for i_data in range(offsetIdx + 2, len(data[i_line])):
        convMetrics.append(float(data[i_line][i_data]))

  ## If no user defined convergence metics, use the default convergence metrics
  if len(convMetrics) == 0:
    convMetrics = [1e-3, 1e-3, 1e-3, 1e-3, 1e-6, 1e-3, 1e-3] ## default values

  ## get the residual names
  convCrit = []
  for i_line in range(len(data)):
    if 'iter' in data[i_line]:
      if not 'scalar-0' in data[i_line]:
        for convName in range(len(data[i_line])):
          convCrit.append(data[i_line][convName])
        break

  ## get the default vs user defined residual names
  endDefault = 0
  for i_conv in range(len(convCrit)):
    if convCrit[i_conv] == 'omega': ## omega SHOULD always be the last default convergence metric
      endDefault = i_conv + 1
      break

  ## get the number of user defined convergence criteria
  nconvCrit_user = 0
  try:
    nconvCrit_user = len(convCrit[endDefault:len(convCrit)-1])
  except IndexError:
    doNothing=True

  ## get the user defined convergence criteria conditions
  user_convNames = []
  for i_conv in range(endDefault, len(convCrit)-1):
    user_convNames.append(convCrit[i_conv])

  user_convMetrics_full = []
  user_convMetrics = []
  for i_line in range(len(data)):
    if 'stop-criterion' in data[i_line]:
      user_currConv = data[i_line][offsetIdx].split('/')[-1].replace('>','')
      for i_conv in range(len(user_convNames)):
        if user_currConv == user_convNames[i_conv]:
          user_convDict = {}
          user_convDict['user_convMetric'] = user_convNames[i_conv]
          user_convDict['stop_condition'] = float(data[i_line][-1])
          user_convMetrics_full.append(user_convDict)
      if len(user_convMetrics_full) == len(user_convNames):
        break

  condFound = False
  isUserMetrics = False
  if len(user_convMetrics_full) == 0:
    isUserMetrics = False
  else:
    if not len(user_convMetrics_full) == len(user_convNames):
      for check in user_convNames:
        condFound = False
        for i_name in range(len(user_convMetrics_full)):
          if user_convMetrics_full[i_name]['user_convMetric'] == check:
            user_convMetrics.append(user_convMetrics_full[i_name])
            condFound = True
            break
        if condFound == False:
          dummy_dict = {}
          dummy_dict['user_convMetric'] = check
          dummy_dict['stop_condition'] = 1e-3
          user_convMetrics_full.append(dummy_dict)
    else:
      isUserMetrics = True

  ## if all values are the same, condense into a single value

  ## PCM: add logic here to consense all similar values into 'user defined convergence criteria'...
  ## that way, if there is a different value, it wont print all names to plot
  if isUserMetrics == True:
    allSame = True
    for i_name in range(len(user_convMetrics_full)):
      if not user_convMetrics_full[i_name]['stop_condition'] == user_convMetrics_full[0]['stop_condition']:
        allSame = False
        break

    if allSame == True:
      user_convMetrics_full = {'stop_condition': user_convMetrics_full[0]['stop_condition']}
    else:
      stopConds = []
      for i_name in range(len(user_convMetrics_full)):
        stopConds.append(user_convMetrics_full[i_name]['stop_condition'])
      default = statistics.mode(stopConds)

      try:
        dummy_CONV = []
        dummy_conv = {}
        dummy_conv['user_convMetric'] = 'user defined convergence criteria'
        dummy_conv['stop_condition'] = default
        dummy_CONV.append(dummy_conv)
        for i_name in range(len(user_convMetrics_full)):
          if not user_convMetrics_full[i_name]['stop_condition'] == default:
            dummy_CONV.append(user_convMetrics_full[i_name])

        ## clear the original list of convergence metrics
        user_convMetrics_full = []
        for i_conv in range(len(dummy_CONV)):
          user_convMetrics_full.append(dummy_CONV[i_conv])

      except statistics.StatisticsError:
        ## If there is no mode, there isnt much to do, it might just be a messy plot
        doNothing = True


  ## get the case data
  caseData_default = {}
  caseData_user = {}

  nIters_default = 1
  nIters_user = 1

  nCases = 0
  convergeVals_default = {}
  convergeVals_user = {}

  ## initialize the default convergence dict
  convergeVals_default = {}
  for i_conv in range(offsetIdx+1,endDefault):
    convergeVals_default[convCrit[i_conv]] = []
  convergeVals_default['iters'] = []
  convergeVals_default['iters_true'] = []

  ## initailize the user convergence dict
  convergeVals_user = {}
  for i_conv in range(endDefault, len(convCrit)-1):
    convergeVals_user[convCrit[i_conv]] = []
  convergeVals_user['iters'] = []
  convergeVals_user['iters_true'] = []

  writeDataFlag = False ## flag to tell the code to stop appending data (for simplicity, convergence data from the 1 iter run from writing report files will not be counted)
  ## loop over every line in data and parse the convergence lines
  for i_line in range(len(data)):

    try:
      if writeDataFlag == False:
        if int(data[i_line][offsetIdx]):
          ## with exactly 1 user defined convergence criteria, parse the line and update the convergenceVals dict
          if nconvCrit_user == 1:
            if len(data[i_line]) == len(convCrit):
              ## only default convergence criteria
              for i_conv in range(offsetIdx+1, endDefault):
                convergeVals_default[convCrit[i_conv]].append(float(data[i_line][i_conv]))
              convergeVals_default['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_default['iters'].append(nIters_default)
              nIters_default += 1

              ## user defined convergence criteria
              for i_conv in range(endDefault, len(convCrit)-1):
                convergeVals_user[convCrit[i_conv]].append(np.nan)
              convergeVals_user['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_user['iters'].append(nIters_user)
              nIters_user += 1

            elif len(data[i_line]) == len(convCrit)+nconvCrit_user:
              ## default convergence criteria AND user convergence criteria(s) present
              for i_conv in range(offsetIdx+1, endDefault):
                convergeVals_default[convCrit[i_conv]].append(float(data[i_line][i_conv]))
              convergeVals_default['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_default['iters'].append(nIters_default)
              nIters_default += 1

              ## user defined convergence criteria
              for i_conv in range(endDefault, len(convCrit)-1):
                convergeVals_user[convCrit[i_conv]].append(float(data[i_line][i_conv]))
              convergeVals_user['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_user['iters'].append(nIters_user)
              nIters_user += 1


          ## with no user defined convergence criteria, parse the line and update the convergenceVals dict
          elif nconvCrit_user == 0:
            if len(data[i_line]) == len(convCrit)+1:
              for i_conv in range(offsetIdx+1, endDefault):
                convergeVals_default[convCrit[i_conv]].append(float(data[i_line][i_conv]))

              convergeVals_default['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_default['iters'].append(nIters_default)
              nIters_default += 1

          ## with n+1 user defined convergence criteria, parse the line and update the convergenceVals dict
          else:
            if len(data[i_line]) == len(convCrit)+1:
              ## default only convergence criteria
              for i_conv in range(offsetIdx+1, endDefault):
                convergeVals_default[convCrit[i_conv]].append(float(data[i_line][i_conv]))
              convergeVals_default['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_default['iters'].append(nIters_default)
              nIters_default += 1

              ## user defined convergence criteria
              for i_conv in range(endDefault,len(convCrit)-1):
                convergeVals_user[convCrit[i_conv]].append(float(data[i_line][i_conv]))
              convergeVals_user['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_user['iters'].append(nIters_user)
              nIters_user += 1


            elif len(data[i_line]) == endDefault+2:
              ## default convergence criteria
              for i_conv in range(offsetIdx+1, endDefault):
                convergeVals_default[convCrit[i_conv]].append(float(data[i_line][i_conv]))
              convergeVals_default['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_default['iters'].append(nIters_default)
              nIters_default += 1

              ## user defined convergence criteria
              for i_conv in range(endDefault,len(convCrit)-1):
                convergeVals_user[convCrit[i_conv]].append(np.nan)
              convergeVals_user['iters_true'].append(int(data[i_line][offsetIdx]))
              convergeVals_user['iters'].append(nIters_user)
              nIters_user += 1


    except (IndexError, ValueError): ## when encountering an empty line, do nothing
      doNothing = True


    if 'file/write-case-data' in data[i_line]: ## when writing the case data
      if writeDataFlag == False: ## If there are 2 consecutive writeDataFlag == True, it means you are writing a restart file
        ## Update the default converge criteria dictionary
        caseName = str(re.findall(r'"(.*?)"', data[i_line][offsetIdx+2]))
        caseData_default['case_'+str(nCases)] = {}
        caseData_default['case_'+str(nCases)]['caseName'] = str(caseName).replace('[\'','').replace('\']','')

        caseData_user['case_'+str(nCases)] = {}
        caseData_user['case_'+str(nCases)]['caseName'] = str(caseName).replace('[\'','').replace('\']','')

        for i_conv in range(offsetIdx+1, endDefault):
          caseData_default['case_'+str(nCases)][convCrit[i_conv]] = convergeVals_default[convCrit[i_conv]]

        caseData_default['case_'+str(nCases)]['iters'] = convergeVals_default['iters']
        caseData_default['case_'+str(nCases)]['iters_true'] = convergeVals_default['iters_true']

        for i_conv in range(endDefault, len(convCrit)-1):
          caseData_user['case_'+str(nCases)][convCrit[i_conv]] = convergeVals_user[convCrit[i_conv]]

        caseData_user['case_'+str(nCases)]['iters'] = convergeVals_user['iters']
        caseData_user['case_'+str(nCases)]['iters_true'] = convergeVals_user['iters_true']

        ## clear the default convergence values dictionary for the next case
        convergeVals_default = {}
        for i_conv in range(offsetIdx+1,endDefault):
          convergeVals_default[convCrit[i_conv]] = []
        convergeVals_default['iters'] = []
        convergeVals_default['iters_true'] = []

        ## clear the user convergence values dictionary for the next case
        convergeVals_user = {}
        for i_conv in range(endDefault, len(convCrit)-1):
          convergeVals_user[convCrit[i_conv]] = []
        convergeVals_user['iters'] = []
        convergeVals_user['iters_true'] = []

        ## update the write flag
        writeDataFlag = True

        nCases += 1
        nIters_default = 1
        nIters_user = 1

    ## continue reading file to prevent the iterations of the report file being added to the next case
    if writeDataFlag == True:
      try:
        if 'define/boundary-conditions/' in data[i_line][offsetIdx+1]: ## this will catch either setting the farfield or setting the backpressure
          ## reset the writeDataFlag
          writeDataFlag = False

      except IndexError:
        doNothing = True


################################################################
## Create the default convergence criteria plot
################################################################
  default_colors = ['cyan', 'mediumpurple', 'red', 'deepskyblue', 'orange', 'greenyellow', 'pink'] ## these colors are intended to match the default colors that fluent uses
  for i_case in range(nCases):
    titleName = str(caseData_default['case_'+str(i_case)]['caseName']) + ' Residuals'
    iters = caseData_default['case_'+str(i_case)]['iters']

    A = []
    labels = []
    conv_metric_idx = 0

    plt.figure(figsize=(24, 12))
    ## unpack the convergenceVals dict into a nested loop for easier plotting
    for i_conv in range(offsetIdx+1,endDefault):
      A.append(caseData_default['case_'+str(i_case)][convCrit[i_conv]])
      resid_val = caseData_default['case_'+str(i_case)][convCrit[i_conv]][-1]

      if isinstance(convMetrics[1], list):
        if resid_val <= convMetrics[i_case-1][conv_metric_idx]:
          converged = True
        else:
          converged = False
      else:
        if resid_val <= convMetrics[conv_metric_idx]:
          converged = True
        else:
          converged = False
      lab = str(convCrit[i_conv]) + ' (' + str(format(resid_val, '.3e')) + '), Converged = ' + str(converged)
      labels.append(lab)

      conv_metric_idx += 1

    defaultMetricApplied = False
    for j in range(len(A)):
      plt.plot(iters, A[j], color = default_colors[j], label = labels[j])
      if isinstance(convMetrics[1], list):
        if convMetrics[i_case-1][j] == 1e-3:
          if defaultMetricApplied == False:
            plt.plot([iters[0], iters[-1]], [convMetrics[i_case-1][j], convMetrics[i_case-1][j]], 'k--')
            text = 'default convergence metric: ' + str(format(convMetrics[i_case-1][j], '.2e'))
            plt.text(iters[0], convMetrics[i_case-1][j]*1.1, text)
            defaultMetricApplied = True
        else:
          plt.plot([iters[0], iters[-1]], [convMetrics[i_case-1][j], convMetrics[i_case-1][j]], '--', color = default_colors[j])
          text = str(convCrit[j+offsetIdx+1]) + ' convergence metric: ' + str(format(convMetrics[i_case-1][j], '.2e'))
          plt.text(iters[0], convMetrics[i_case-1][j]*1.1, text)

      else:
        if convMetrics[j] == 1e-3:
          if defaultMetricApplied == False:
            plt.plot([iters[0], iters[-1]], [convMetrics[j], convMetrics[j]], 'k--')
            text = 'default convergence metric: ' + str(format(convMetrics[j], '.2e'))
            plt.text(iters[0], convMetrics[j]*1.1, text)
            defaultMetricApplied = True
        else:
          plt.plot([iters[0], iters[-1]], [convMetrics[j], convMetrics[j]], '--', color = default_colors[j])
          text = str(convCrit[j+offsetIdx+1]) + ' convergence metric: ' + str(format(convMetrics[j], '.2e'))
          plt.text(iters[0], convMetrics[j]*1.1, text)

    plt.yscale('log')


    ## NOTE: If the difference between the smallest and largest residual is greater than 8 orders of magnitude, minor tick marks on the y-axis will not appear
    flat_vals = [val for inner in A for val in inner]
    minY = min(flat_vals) / 10
    maxY = max(flat_vals) * 10
    plt.ylim(minY, maxY)

    ax = plt.gca()
    # format the y-axis
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))
    ax.yaxis.set_minor_formatter(NullFormatter())

    ## format the x-axis
    ax.xaxis.set_minor_locator(MultipleLocator(25))

    # Custom formatter: show label only if divisible by 50
    def label_every_50(x, pos):
        return f"{int(x)}" if x % 50 == 0 else ""

    ax.xaxis.set_minor_formatter(FuncFormatter(label_every_50))
    ax.xaxis.set_major_locator(MultipleLocator(100))

    plt.grid( which='major', axis='y')
    ax.tick_params(axis='y', which='minor', length=4, color='gray')
    plt.legend(loc = 'center left', bbox_to_anchor =(1, 0.5))
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)

    plt.xlabel('Iterations')
    plt.ylabel('Residuals')

    if not caseData_default['case_'+str(i_case)]['iters_true'][0] == 1:
      titleName += '\n(Solution started from a restart)'
    plt.title(titleName)
    plt.subplots_adjust(top=0.9, bottom=0.1, left = 0.05)

    ## to save figures instead of showing them, replace plt.show() with plt.savefig(saveName), where saveName is the uniqie path + Name of the picture
    saveName = os.path.join(residual_plotDir, str(caseData_default['case_'+str(i_case)]['caseName']) +'_default_residuals_plot.png')
    plt.savefig(saveName)
    plt.close()


  ##########################################################################
  ## Create the User defined convergence criteria plot (ONLY if they exist)
  if isUserMetrics == True:
    default_colors = plt.get_cmap('hsv', len(convCrit) - endDefault)
    for i_case in range(nCases): ##NOTE: since nCases updates when writing cases data, nCases will be 1 more than the total number of cases, hence the -1
      titleName = str(caseData_user['case_'+str(i_case)]['caseName']) + ' Residuals'
      iters = caseData_user['case_'+str(i_case)]['iters_true']

      A = []
      labels = []
      conv_metric_idx = 0

      ## unpack the convergenceVals dict into a nested loop for easier plotting
      for i_conv in range(endDefault, len(convCrit)-1):
        A.append(caseData_user['case_'+str(i_case)][convCrit[i_conv]])
        resid_val = caseData_user['case_'+str(i_case)][convCrit[i_conv]][-1]


        if isinstance(user_convMetrics_full, list):
          for i_name in range(len(user_convMetrics_full)):
              if user_convMetrics_full[i_name]['user_convMetric'] == convCrit[i_conv]:
                break
          if resid_val <= user_convMetrics_full[i_name]['stop_condition']:
            converged = True
          else:
            converged = False
        else:
          if resid_val <= user_convMetrics_full['stop_condition']:
            converged = True
          else:
            converged = False
        lab = str(convCrit[i_conv]) + ' (' + str(format(resid_val, '.3e')) + '), Converged = ' + str(converged)
        labels.append(lab)

        conv_metric_idx += 1


      if not all(np.isnan(val) for currConv in A for val in currConv): ## make sure the list is not only nans
        plt.figure(figsize=(24, 12))
        defaultMetricApplied = False
        for j in range(len(A)):
          plt.plot(iters, A[j], color = default_colors(j), label = labels[j])
        if isinstance(user_convMetrics_full, list):
          for i_metric in range(len(user_convMetrics_full)):
            if user_convMetrics_full[i_metric]['user_convMetric'] == 'user defined convergence criteria':
              plt.plot([iters[0], iters[-1]], [user_convMetrics_full[i_metric]['stop_condition'], user_convMetrics_full[i_metric]['stop_condition']], 'k--')
            else:
              plt.plot([iters[0], iters[-1]], [user_convMetrics_full[i_metric]['stop_condition'], user_convMetrics_full[i_metric]['stop_condition']], '--', color = default_colors(j))
            text = str(user_convMetrics_full[i_metric]['user_convMetric']) + ' convergence metric: ' + str(format(user_convMetrics_full[i_metric]['stop_condition'], '.2e'))
            plt.text(iters[0], user_convMetrics_full[i_metric]['stop_condition']*1.025, text)

        else:
          if user_convMetrics_full['stop_condition'] == 1e-3:
            if defaultMetricApplied == False:
              plt.plot([iters[0], iters[-1]], [user_convMetrics_full['stop_condition'], user_convMetrics_full['stop_condition']], 'k--')
              text = 'default convergence metric: ' + str(format(convMetrics[j], '.2e'))
              plt.text(iters[0], user_convMetrics_full['stop_condition']*1.025, text)
              defaultMetricApplied = True
          else:
            plt.plot([iters[0], iters[-1]], [user_convMetrics_full['stop_condition'], user_convMetrics_full['stop_condition']], 'k--')
            text = 'User defined convergence metric: ' + str(format(user_convMetrics_full['stop_condition'], '.2e'))
            plt.text(iters[0], user_convMetrics_full['stop_condition']*1.025, text)
        plt.yscale('log')

        ## NOTE: If the difference between the smallest and largest residual is greater than 8 orders of magnitude, minor tick marks on the y-axis will not appear
        flat_vals = [val for inner in A for val in inner if not np.isnan(val)]
        minY = min(flat_vals) / 10
        maxY = max(flat_vals) * 10
        plt.ylim(minY, maxY)

        ax = plt.gca()
        # format the y-axis
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))
        ax.yaxis.set_minor_formatter(NullFormatter())

        ## format the x-axis
        ax.xaxis.set_minor_locator(MultipleLocator(25))
        plt.xlim(min(iters), max(iters))

        # Custom formatter: show label only if divisible by 50
        def label_every_50(x, pos):
            return f"{int(x)}" if x % 50 == 0 else ""
        ax.xaxis.set_minor_formatter(FuncFormatter(label_every_50))
        ax.xaxis.set_major_locator(MultipleLocator(100))

        plt.grid( which='major', axis='y')
        ax.tick_params(axis='y', which='minor', length=4, color='gray')
        plt.legend(loc = 'center left', bbox_to_anchor =(1, 0.5))
        plt.tight_layout()
        plt.subplots_adjust(right=0.75)

        plt.xlabel('Iterations')
        plt.ylabel('Residuals')

        if not caseData_default['case_'+str(i_case)]['iters_true'][0] == 1:
          titleName += '\n(Solution started from a restart)'
        plt.title(titleName)
        plt.subplots_adjust(top=0.9, bottom=0.1, left = 0.05)
        ## to save figures instead of showing them, replace plt.show() with plt.savefig(saveName), where saveName is the unique path + Name of the picture
        saveName = os.path.join(residual_plotDir, str(caseData_user['case_'+str(i_case)]['caseName']) +'_user_defined_residuals_plot.png')
        plt.savefig(saveName)
        plt.close()

  ## Before exiting the function, reset the plot configurations
  plt.rcParams.update(plt.rcParamsDefault)
####################################################
## END CREATE RESIDUAL PLOTS
####################################################



####################################################
## CREATE CANE CURVES
####################################################
def Create_Cane_Curves(plotDir, df, interface_surfs):

  print('\n>> Generating Cane Curves...')

  INTERFACES = []
  for i_int in range(len(interface_surfs)):
    interface = interface_surfs[i_int].split('_')[-1]
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
  vmag_0 = []
  u_0 = list(df['u_0'])
  v_0 = list(df['v_0'])
  w_0 = list(df['w_0'])
  for i_vel in range(len(u_0)):
    vmag_0.append(math.sqrt(u_0[i_vel]**2 + v_0[i_vel]**2 + w_0[i_vel]**2))
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

    ## the critical point will be defined using the mdot and P_static at the AIP
    ## WHAT HAPPENS IF THE AIP is not present ???
    aipFound = False
    for int_surf in INTERFACES:
      if 'aip' in int_surf.lower() or 'api' in int_surf.lower():
        mdot_ST = list(df[df['mach_0_BC'] == i_mach]['mdot_' + str(int_surf)])
        p_ST = list(df[df['mach_0_BC'] == i_mach]['p_' + str(int_surf)])
        UNSTART_VARS=Calculate_critical_value(mdot_ST, p_ST)
        aipFound = True
        break

    if aipFound == False:
      raise Exception('NO AIP PLANE FOUND')

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
      vmag_0 = []
      u_0 = list(df[df['mach_0_BC'] == i_mach]['u_0'])
      v_0 = list(df[df['mach_0_BC'] == i_mach]['v_0'])
      w_0 = list(df[df['mach_0_BC'] == i_mach]['w_0'])
      for i_vel in range(len(u_0)):
        vmag_0.append(math.sqrt(u_0[i_vel]**2 + v_0[i_vel]**2 + w_0[i_vel]**2))
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
      maxBP = backPress_ST[UNSTART_VARS[1]]
      caseVars = [currMach, currAlt, currInterface, maxBP]

      ## plot the different cane curves at each station
      Plot_Cane_Curves(currIntPlotDir, x_var = mdot_ST, y_var = p_ST, x_name = 'Mdot', y_name = 'static_pressure', case_vars = caseVars, unstart_vars=UNSTART_VARS)
      Plot_Cane_Curves(currIntPlotDir, x_var = mdot_ST, y_var = pt_ST, x_name = 'Mdot', y_name = 'total_pressure', case_vars = caseVars, unstart_vars=UNSTART_VARS)

      # calculate MFR
      MFR = []
      PTR = []
      PR = []
      for i in range(len(rho_0)):
        MFR.append(mdot_ST[i]/ (rho_0[0] * vmag_0[0] * 0.00553)) ## PCM: NEED TO GET AREF FROM CODE. THIS CAN NOT BE HARD CODED!!!!!
        PTR.append(pt_ST[i] / pt_0[i])
        PR.append(p_ST[i] / p_0[i])

      Plot_Cane_Curves(currIntPlotDir, x_var = MFR, y_var = PTR, x_name = 'MFR', y_name = 'PTR', case_vars = caseVars, unstart_vars=Calculate_critical_value(MFR, PTR))
      Plot_Cane_Curves(currIntPlotDir, x_var = MFR, y_var = PR, x_name = 'MFR', y_name = 'PR', case_vars = caseVars, unstart_vars=Calculate_critical_value(MFR, PR))

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

    if abs(x_var[i]) >= 1.5*critical:
      break
      #elif x_var[i] >= critical and y_var[i] <= y_var[i-1]:
      #  max_idx = i

  return([critical, max_idx])

def Plot_Cane_Curves(plotDir, x_var = [], y_var = [], x_name = '', y_name = '', case_vars = None, unstart_vars=None):
  ## reset plot parameters to default
  plt.rcParams.update(plt.rcParamsDefault)

  mach = case_vars[0]
  alt = case_vars[1]
  interface = case_vars[2]
  maxBP = case_vars[3]

  critical = unstart_vars[0]
  max_idx = unstart_vars[1]

  plt.figure(figsize=(24, 12))
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

  ## CREATE PLOT ZOOMED IN ON THE CRITICAL POINT
  plt.figure(figsize=(24, 12))
  plt.grid(True, which = 'major')
  plt.plot([critical, critical], [min(y_var), max(y_var)], 'k--', label = 'Unstart limit (2% deviation)')
  plt.plot(x_var, y_var, 'x-')
  plt.plot(x_var[max_idx], y_var[max_idx], 'rx', label = 'Critical Point. (Back Pressure: ' + str(maxBP) + ' Pa)')
  title = 'Mach ' +str(mach) + ' Alt ' +str(alt) + ' ft\nStation ' + str(interface) + ': ' + str(x_name) +' vs. ' +str(y_name)
  plt.title(title)
  plt.xlabel(x_name)
  plt.ylabel(y_name)
  plt.xlim(0.95*critical, 1.05*critical)
  plt.ylim(y_var[0], y_var[max_idx]*1.01)

  plt.legend()
  plt.savefig(os.path.join(plotDir, 'Station_' + str(interface) + '_' +str(x_name) + '_vs_' + str(y_name) +'_zoomed_critical.png'))
  plt.close()




  ####################################################
  ## END CREATE CANE CURVES
  ####################################################