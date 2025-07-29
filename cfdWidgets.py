# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 07:16:56 2025

@author: mfund
"""

import math

def make_rescale_cli_submission(f,rescale_name='test_job',journal_name='journal_name',rescale_cores=8,rescale_analysis='ansys_fluent',rescale_analysis_version='2024r2',rescale_core_type='Emerald',rescale_walltime_hours=10,rescale_program_flags=None,rescale_project_name='Invictus Aero'):

    header_string = """#!/bin/bash
#RESCALE_NAME={rescale_name}
#RESCALE_CORES={rescale_cores}
#RESCALE_ANALYSIS={rescale_analysis}
#RESCALE_ANALYSIS_VERSION={rescale_analysis_version}
#RESCALE_CORE_TYPE={rescale_core_type}
#RESCALE_ENV_ANSYSLMD_LICENSE_FILE=1055@tiberius-ansys
#RESCALE_ENV_ANSYSLI_SERVERS=1055@tiberius-ansys
#{rescale_license_queue}
#RESCALE_WALLTIME={rescale_walltime_hours}
#RESCALE_PROJECT_ID="{rescale_project_name}"
"""


    if rescale_analysis == 'ansys_fluent':
        rescale_base_feature = 'cfd_base'

    if rescale_core_type == 'Emerald':
        rescale_cores = math.ceil(rescale_cores/36)*36
    elif rescale_core_type == 'Luna':
        rescale_cores = math.ceil(rescale_cores/48)*48

    if 'ansys' in rescale_analysis:
        if rescale_cores <= 12:
            rescale_hpc_packs = 1
        elif rescale_cores <= 36:
            rescale_hpc_packs = 2
        elif rescale_cores <= 132:
            rescale_hpc_packs = 3
        elif rescale_cores <= 516:
            rescale_hpc_packs = 4
        elif rescale_hpc_packs <= 2052:
            rescale_hpc_packs = 5

    replacements = {
        'rescale_name':rescale_name,
        'rescale_cores':str(int(rescale_cores)),
        'rescale_analysis':rescale_analysis,
        'rescale_analysis_version':rescale_analysis_version,
        'rescale_core_type':rescale_core_type,
        'rescale_walltime_hours':str(int(rescale_walltime_hours)),
        'rescale_project_name':rescale_project_name,
        'rescale_license_queue':'RESCALE_USER_DEFINED_LICENSE_SETTINGS={"featureSets":[{"name":"USER_SPECIFIED","features":[{"name":"{rescale_base_feature}","count":1},{"name":"anshpc_pack","count":{rescale_hpc_packs}}]}]}}]}'.replace('{rescale_base_feature}',rescale_base_feature).replace('{rescale_hpc_packs}',str(int(rescale_hpc_packs)))
        }

    header_string = header_string.format(**replacements)

    f.writelines(header_string+'\n')
    f.writelines('\n')

    if 'fluent' in rescale_analysis:
        if rescale_program_flags == None:
            rescale_program_flags = '3ddp'
        rescale_program = 'fluent'
        rescale_program_flags = rescale_program_flags + ' -gu -ssh -cnf=$FLUENT_HOSTS -t$RESCALE_CORES_PER_SLOT'

        f.writelines(rescale_program + ' ' + rescale_program_flags + ' -i ' + journal_name + '.jou > process_output.log')

# %%

# f = open('test.txt','w')
# make_rescale_cli_submission(f,rescale_cores=12,rescale_walltime_hours=40,rescale_program_flags='2ddp')
# f.close()
