# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
def standard_init(f):
    f.writelines('solve/init/initialize\n')

def set_standard_init(f,boundary_type='pressure-far-field',boundary_name='pressure_farfield'):
    f.writelines('solve/init/compute-defaults ' + boundary_type + ' ' + boundary_name + '\n')

def set_standard_init_reduced_velocity(f,boundary_type='pressure-far-field',boundary_name='pressure_farfield', velocity=100):
    f.writelines('solve/init/compute-defaults ' + boundary_type + ' ' + boundary_name + '\n')
    f.writelines('solve/init/set-defaults/ x-velocity ' +str(velocity) +'\n')
    f.writelines('solve/init/reference-frame relative\n')

def fmg_init(f):
    f.writelines('solve/init/fmg-init yes\n')

def hyb_init(f):
    f.writelines('solve/init/hyb-init\n')

def set_under_relaxation_temperature(f, temp=1.0):
  f.writelines('solve/set/under-relaxation/temperature ' + str(temp) + ' quit quit\n')

def set_pressure_farfield(f,boundary_name='farfield_fluid',Mach=3,p_Pa=101325,T_K=288.15,TI=0.05,TVR=10,x_component_direction=1,y_component_direction=0,z_component_direction=0,mode='3D',turbulence='SST',species_transport=False,dpm=False):

    if mode=='3D' and turbulence == 'SA':
        string = 'define/boundary-conditions/pressure-far-field ' + boundary_name + ' no ' +\
            str(p_Pa) + ' no ' + str(Mach) + ' no ' + str(T_K) + ' yes no ' + str(x_component_direction) +\
                ' no ' + str(y_component_direction) + ' no ' + str(z_component_direction) +\
                    ' no no yes no ' + str(TVR)
    elif mode=='axisymmetric' and turbulence == 'SST':
        string = 'define/boundary-conditions/pressure-far-field ' + boundary_name + ' no ' +\
            str(p_Pa) + ' no ' + str(Mach) + ' no ' + str(T_K) + ' no ' + str(x_component_direction) +\
                ' no ' + str(y_component_direction) + ' no no yes ' + str(TI) + ' ' + str(TVR)
    elif mode=='3D' and turbulence == 'SST':
        string = 'define/boundary-conditions/pressure-far-field ' + boundary_name + ' no ' +\
            str(p_Pa) + ' no ' + str(Mach) + ' no ' + str(T_K) + ' yes no ' + str(x_component_direction) +\
                ' no ' + str(y_component_direction) + ' no ' + str(z_component_direction) +\
                    ' no no yes ' + str(TI) + ' ' + str(TVR)

    if species_transport:
        string += ' no no 0 no 0.23 no 0 no 0' #weight fractions of diesel, o2, co2, h2o, all else nitrogen (this is an air mix)

    if dpm:
        string += ' no'

    string += '\n'

    f.writelines(string)

def set_mass_flow_inlet(f,boundary_name='fuel_inlet',mdot_kgps=0.45,Tt_K=790,p_Pa=101325,TI=0.05,TVR=10,mode='3D',turbulence='SST',species_transport=False,dpm=False,mixture=False):

    if mode in ['axisymmetric','3D'] and turbulence == 'SST' and not mixture:
         string = 'define/boundary-conditions/mass-flow-inlet ' + boundary_name + ' yes yes no ' +\
             str(mdot_kgps) + ' no ' + str(Tt_K) + ' no ' + str(p_Pa) + ' no yes no no yes ' + str(TI) + ' ' + str(TVR)

    if mode in ['3D'] and turbulence == 'SST' and mixture:
         string = 'define/boundary-conditions/mass-flow-inlet ' + boundary_name + ' phase-1 yes no ' +\
             str(mdot_kgps)

    if species_transport:
        string += ' no no 0 no 0.23 no 0 no 0' #weight fractions of diesel, o2, co2, h2o, all else nitrogen (this is pure fuel)

    if dpm:
        string += ' no'

    string += '\n'

    f.writelines(string)


def set_pressure_inlet(f,boundary_name='pressure_inlet',Tt_K=790,p_Pa=101325,pt_Pa=101325,TI=0.05,TVR=10,mode='3D',turbulence='SST',species_transport=False,dpm=False):

    if mode in ['axisymmetric','3D'] and turbulence == 'SST':
         string = 'define/boundary-conditions/pressure-inlet ' + boundary_name + ' yes no ' +\
             str(pt_Pa) + ' no ' + str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes ' + str(TI) + ' ' + str(TVR)

    if species_transport:
        string += ' no no 0 no 0.23 no 0 no 0' #weight fractions of diesel, o2, co2, h2o, all else nitrogen (this is pure fuel)

    if dpm:
        string += ' no'

    string += '\n'

    f.writelines(string)

def set_pressure_outlet(f,boundary_name='pressure_outlet_fluid',p_Pa=101325,Tt_K=288.15,TI=0.05,TVR=10,turbulence='SA',species_transport=False,dpm=False,target_mass_flow=False,mdot_kgps=1):

    if target_mass_flow:
        p_lower_Pa = 1
        p_upper_Pa = 5e6

    if not target_mass_flow and turbulence == 'SA':
        string = 'define/boundary-conditions/pressure-outlet ' + boundary_name + ' yes no ' +\
            str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes no ' + str(TVR) + \
                ' yes no yes no\n'

    elif not target_mass_flow and turbulence == 'SST':
        string = 'define/boundary-conditions/pressure-outlet ' + boundary_name + ' yes no ' +\
          str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes ' + str(TI) + ' ' + str(TVR) +\
              ' no yes no yes no\n'
            #str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes ' + str(TI) + ' ' + str(TVR) +\
            #    ' yes yes no\n'

    elif target_mass_flow and turbulence == 'SST':
        string = 'define/boundary-conditions/pressure-outlet ' + boundary_name + ' yes no ' +\
            str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes ' + str(TI) + ' ' + str(TVR) +\
                ' yes yes yes no ' + str(mdot_kgps) + ' no ' + str(p_upper_Pa) + ' no ' + str(p_lower_Pa) + '\n'

    elif target_mass_flow and turbulence == 'SA':
        string = 'define/boundary-conditions/pressure-outlet ' + boundary_name + ' yes no ' +\
            str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes no ' + str(TVR) + \
                ' yes no yes yes no ' + str(mdot_kgps) + ' no ' + str(p_upper_Pa) + ' no ' + str(p_lower_Pa) + '\n'

    elif not target_mass_flow and turbulence == 'SST' and species_transport == True and dpm == True:
        string = 'define/boundary-conditions/pressure-outlet ' + boundary_name + ' yes no ' +\
            str(p_Pa) + ' no ' + str(Tt_K) + ' no yes no no yes ' + str(TI) + ' ' + str(TVR) + \
                ' no no 0 no 0.23 no 0 no 0 no yes no no no\n'

    f.writelines(string)

def set_pressure_outlet_2D(f, pressure, Tt_K):
  string = 'define/boundary-conditions/pressure-outlet outlet yes no ' +\
    str(pressure) + ' no ' + str(Tt_K) + ' no yes no no yes ' + str(0.001) + ' ' + str(10) +\
        ' no yes yes no\n'
  f.writelines(string)

def solve(f,iterations=10):
    iterations = int(iterations)
    f.writelines('solve/iterate ' + str(iterations) + '\n')

def read_case(f,case_name='Inlet_Aero_155_Initial_v015_flowthrough'):
    f.writelines('file/read-case "' + str(case_name) + '"\n')

def read_data(f,data_name='Inlet_Aero_155_Initial_v015_flowthrough'):
    f.writelines('file/read-data "' + str(data_name) + '"\n')

def read_case_data(f,case_name='Inlet_Aero_155_Initial_v011_Full_Body'):
    f.writelines('file/read-case-data "' + str(case_name) + '"\n')

def set_CFL(f,CFL=5): #must be in implicit timestepping mode
    momentum_URF = 0.5
    pressure_URF = 0.5

    f.writelines('solve/set/p-v-controls ' + str(CFL) + ' ' + str(momentum_URF)\
                 + ' ' + str(pressure_URF) + '\n')

def set_tui_version(f):
    tui_version = '"24.2"'
    f.writelines('file/set-tui-version ' + str(tui_version) + '\n')

def set_batch_options(f):
    f.writelines('file/set-batch-options no no yes no \n')

def write_case_data(f,case_name='Inlet_Aero_155_Initial_v015_flowthrough'):
    f.writelines('file/write-case-data "' + str(case_name) + '"\n')

def write_output_parameters(f):
    output_name = 'Inlet_Aero_155_Initial_v011_Full_Body'
    f.writelines('define/parameters/output-parameters/write-all-to-file "' + str(output_name) + '.out"\n')

def set_reference_values(f,boundary_type='pressure-far-field',boundary_name='farfield_fluid'):
    location_name = boundary_type + ' ' + boundary_name
    f.writelines('report/reference-values/compute/' + str(location_name) + '\n')

def set_CFL_ramp(f,CFLs=[40],iterations=[10]):
    iterable=list(zip(CFLs,iterations))

    for i_step in range(len(iterable)):
        set_CFL(f,iterable[i_step][0])
        solve(f,iterable[i_step][1])

def exit_job(f):
    f.writelines('exit\n')

def check_convergence(f,boolean=True,mode='3D',turbulence='SST',species_transport=False,species_off=False):
    if boolean and mode=='axisymmetric' and turbulence=='SST' and species_transport and not species_off:
        f.writelines('solve/monitors/residual/check-convergence yes yes yes yes yes yes yes yes yes yes\n')

    elif boolean and mode=='axisymmetric' and turbulence=='SST' and not species_transport:
        f.writelines('solve/monitors/residual/check-convergence yes yes yes yes yes yes\n')
    elif not boolean and mode=='axisymmetric' and turbulence=='SST' and not species_transport:
        f.writelines('solve/monitors/residual/check-convergence no no no no no no\n')

    elif boolean and mode=='axisymmetric' and turbulence=='SST' and species_transport and species_off:
        f.writelines('solve/monitors/residual/check-convergence yes yes yes yes yes yes no no no no\n')
    elif not boolean and mode=='axisymmetric' and turbulence=='SST' and species_transport:
        f.writelines('solve/monitors/residual/check-convergence no no no no no no no no no no\n')

    elif boolean and mode == '3D'and turbulence=='SST' and species_transport and not species_off:
        f.writelines('solve/monitors/residual/check-convergence yes yes yes yes yes yes yes yes yes yes yes\n')
    elif boolean and mode == '3D'and turbulence=='SST' and species_transport and species_off:
        f.writelines('solve/monitors/residual/check-convergence yes yes yes yes yes yes yes no no no no\n')
    elif not boolean and mode == '3D'and turbulence=='SST' and species_transport:
        f.writelines('solve/monitors/residual/check-convergence no no no no no no no no no no no\n')

def change_boundary_type(f,boundary_name='fuel_inlet',boundary_type='mass-flow-inlet'):
    string = 'define/boundary-conditions/modify-zones/zone-type ' + boundary_name + ' ' + boundary_type + '\n'
    f.writelines(string)

def set_AMR(f,iterations=500): #for shock refinement only
    f.writelines('mesh/adapt/predefined-criteria/aerodynamics/shock-indicator density-based\n')
    f.writelines('mesh/adapt/manage-criteria/edit/density_0 frequency ' + str(iterations) + '\n')

def patch(f,cell_register='ignitor',quantity='temperature',value=2500,custom_field_functions = False,no_attempts=1):
    if not custom_field_functions:
        for i in range(no_attempts):
            f.writelines('solve/patch/() ' + cell_register + ' () ' + quantity + ' ' + str(value) + '\n')
            solve(f)
    elif custom_field_functions:
        for i in range(no_attempts):
            f.writelines('solve/patch/() ' + cell_register + ' () ' + quantity + ' no ' + str(value) + '\n')
            solve(f)

def make_force_report(f,report_name,force_vector,surfaces):
    f.writelines('solve/report-definitions add ' + report_name + ' force\n')

    string = 'thread-names '
    for i in range(len(surfaces)):
        string+=surfaces[i] + ' '
    string+='()\n'

    f.writelines(string)

    string = 'force-vector '
    for i in range(len(force_vector)):
        string+=str(force_vector[i]) + ' '
    string+='() quit\n'

    f.writelines(string)

def make_moment_report(f,report_name,moment_vector,MRC,surfaces):
    f.writelines('solve/report-definitions add ' + report_name + ' moment\n')

    string = 'thread-names '
    for i in range(len(surfaces)):
        string+=surfaces[i] + ' '
    string+='()\n'

    f.writelines(string)

    string = 'mom-axis '
    for i in range(len(moment_vector)):
        string+=str(moment_vector[i]) + ' '
    string+='()\n'

    f.writelines(string)

    string = 'solve/report-definition edit ' + str(report_name) + ' mom-center '
    for i in range(len(MRC)):
        string+=str(MRC[i]) + ' '
    string+='quit\n'

    f.writelines(string)

    string = 'solve/report-definition edit ' + str(report_name) + ' scaled no quit\n'

    f.writelines(string)

def change_report_file_name(f,report_file_name):
    f.writelines('solve/report-files/edit reports file-name ' + report_file_name + '.out quit\n')

def change_all_report_file_names(f, report_file_name, report_file_num):
  f.writelines('solve/report-files/edit ' + str(report_file_num) + ' file-name ' + str(report_file_name) + '_' + str(report_file_num) + '.out quit\n')

def change_residual_file_name(f,residual_file_name):
    f.writelines('plot/residuals-set/plot-to-file ' + residual_file_name + '.out\n')

def change_turb_chem_int(f,eddy_dissipation=False):
    if eddy_dissipation:
        f.writelines('define/models/species/set-turb-chem-interaction no yes\n')
    else:
        f.writelines('define/models/species/set-turb-chem-interaction yes\n')

def make_fuel_injector_dpm(f,injector_name,boundary_name,mdot_fuel_kgps,V_fuel_mps,T_fuel_K,D_fuel_m):
    f.writelines('define/models/dpm/interaction coupled-calculations yes \n') #turns on interaction with continuous phase (needed for particle evap to vapor fuel)
    f.writelines('define/models/dpm/injections create-injection ' + str(injector_name) + ' yes droplet yes surface diesel-liquid ' + str(boundary_name) + ' () no yes no no no no c10h22 no '\
                 + str (D_fuel_m) + ' ' + str(T_fuel_K) + ' ' + str(V_fuel_mps) + ' ' + str(mdot_fuel_kgps) + '\n')

def modify_fuel_injector_dpm(f,injector_name,boundary_name,mdot_fuel_kgps,V_fuel_mps,T_fuel_K,D_fuel_m):
    f.writelines('define/models/dpm/injections set-injection-properties ' + str(injector_name) + ' ' + str(injector_name) + ' no no no ' + str(boundary_name) + ' () yes yes no no no no no no '\
                 + str (D_fuel_m) + ' ' + str(T_fuel_K) + ' ' + str(V_fuel_mps) + ' ' + str(mdot_fuel_kgps) + '\n')

def make_report(f,report_name,report_type,report_surfaces,report_field):

    surface_str = ''
    for i_surface in range(len(report_surfaces)):
        surface_str += report_surfaces[i_surface] + ' '

    if report_type != 'surface-area' and report_type != 'surface-massflowrate':
        f.writelines('solve/report-definitions add ' + str(report_name) + ' ' + report_type + ' surface-names ' + surface_str + \
                     '() field ' + str(report_field) + ' quit\n')

    else:
        f.writelines('solve/report-definitions add ' + str(report_name) + ' ' + report_type + ' surface-names ' + surface_str + \
                     '()' + ' quit\n')

    #surface-massavg
    #surface-massflowrate

def make_report_file(f, freq, report_name, fname):
  if len(report_name) > 0:
    f.writelines('solve/report-files add ' + str(fname) + ' quit\n')
    f.writelines('solve/report-files edit ' + str(fname) + ' frequency ' + str(freq) + ' quit\n')
    f.writelines('solve/report-files edit ' + str(fname) +' report-defs ')
    for i_surface in range(len(report_name)):
      f.writelines(str(report_name[i_surface])+' ')
    f.writelines('\n\nquit\n')

def make_report_file_to_monitor(f, surfs):
  f.writelines('solve/report-files add reports_for_monitoring_only quit\n')
  f.writelines('solve/report-files edit reports_for_monitoring_only frequency 1 quit\n')
  f.writelines('solve/report-files edit reports_for_monitoring_only report-defs ')
  for i_surf in surfs:
    f.writelines(i_surf + ' ')
  f.writelines('\n\nactive? yes\n')
  f.writelines('write-instantaneous-values? yes\n')
  f.writelines('print? yes quit quit quit quit\n')

def write_force_files(f, fName):
  f.writelines('report/forces/wall-forces/ yes 1 0 0 yes \"' + str(fName) + '\" quit quit\n')

def make_field_function(f,function_name,function):
    f.writelines('define/custom-field-functions define "' + str(function_name) + '" "' + str(function) + '" quit\n')

def add_reax_to_mixture(f,eddy_dissipation_A=4,eddy_dissipation_B=0.5):
    f.writelines(f'define/materials change-create diesel-air diesel-air no yes 1 reaction-1 2 c10h22 1 0.25 o2 15.5 1.5 2 co2 10 0 h2o 11 0 {eddy_dissipation_A} {eddy_dissipation_B} no no no no no no no\n')

def add_reax_to_mixture_stull(f,eddy_dissipation_A=4,eddy_dissipation_B=0.5):
    f.writelines(f'define/materials change-create kerosene-air kerosene-air no yes 1 reaction-1 2 c12h23 1 0.25 o2 17.75 1.5 2 co2 12 0 h2o 11.5 0 {eddy_dissipation_A} {eddy_dissipation_B} no no no no no no no\n')

def write_converged_reports(f,report_file_name):
    if isinstance(report_file_name, str):
      f.writelines('solve/report-files/edit ' + str(report_file_name) + ' frequency 1 quit\n' )
    elif isinstance(report_file_name, list):
      for i_file in report_file_name:
        f.writelines('solve/report-files/edit ' + str(i_file) + ' frequency 1 quit\n' )
    else:
      raise Exception('!! ERROR !! Invalid report_file_name type passed to write_converged_reports(). The only valid types are str or list.')

    solve(f,iterations=1)

    if isinstance(report_file_name, str):
      f.writelines('solve/report-files/edit ' + str(report_file_name) + ' frequency 100 quit\n' )
    elif isinstance(report_file_name, list):
      for i_file in report_file_name:
        f.writelines('solve/report-files/edit ' + str(i_file) + ' frequency 100 quit\n' )

def make_iso_surface(f,surface_of='pressure-coefficient',surface_name='Cp_eta_0p1_wall_wing_upper',surface='wall_wing_upper',zone='volume_zone',value=-0.7):
    f.writelines('surface/iso-surface ' + surface_of + ' ' + surface_name + ' ' + surface + ' () ' + zone + '() ' + str(value) + ' () quit\n')

def write_strip_quantity(f,filename,quantity='pressure-coefficient',vector=[1,0,0],surface='wall_wing_upper_eta_0p0'):
    f.writelines('plot/plot yes ' + filename + '.xy yes no no ' + quantity + ' yes ' + str(vector[0]) + ' ' + str(vector[1]) + ' ' + str(vector[2]) + ' ' + surface + ' () quit\n')

def activate_report_file(f,report_name,active='yes'):
    f.writelines('solve/report-files/edit/' + report_name + ' active ' + active + ' quit quit quit \n')

def make_isothermal_walls(f,wall_names,T_K): #walls must initially be adiabatic
    string = ''
    for i_wall in range(len(wall_names)):
        string += wall_names[i_wall] + ' '
    f.writelines('define/boundary-conditions/wall ' + string + '() no 0 no 0 no yes temperature no ' + str(T_K) + ' no no no no no 0 no 0.5 no polynomial 1 1 polynomial 1 1 no yes yes yes yes no 0 no 0 no 0 no 0 no 1\n')

def export_solution_data_ascii(f,file_name,wall_names,var_names,cell_centered=False):
    if cell_centered:
        cell_centered_str = 'yes'
    else:
        cell_centered_str = 'no'

    wall_str = ''
    for i_wall in range(len(wall_names)):
        wall_str+=wall_names[i_wall] + ' '

    var_str = ''
    for i_var in range(len(var_names)):
        var_str+=var_names[i_var] + ' '

    f.writelines('file/export/ascii ' + file_name + '.csv ' + wall_str + '() yes ' + var_str + '() ' + cell_centered_str + '\n')

def make_convergence_criterion(f,criterion_name,report_name,previous_iterations,stop_value,active=True,console=True,ignore_iterations_before=None):
    f.writelines('solve/convergence-conditions/conv-reports/add ' + criterion_name + '\n')
    f.writelines('report-defs ' + report_name + '\n')

    if not ignore_iterations_before == None:
      f.writelines('initial-values-to-ignore ' + str(ignore_iterations_before) + '\n')
    f.writelines('previous-values-to-consider ' + str(previous_iterations) + '\n')
    f.writelines('stop-criterion ' + str(round(stop_value,3)) + '\n')

    if console:
        f.writelines('print yes\n')
    else:
        f.writelines('print no\n')

    if active:
        f.writelines('active? yes quit quit quit quit\n')
    else:
        f.writelines('active? no quit quit quit quit\n')

def modify_convergence_criterion(f,criterion_name,stop_value,active=True,console=True):
    f.writelines('solve/convergence-conditions/conv-reports/edit ' + criterion_name + '\n')
    f.writelines('stop-criterion ' + str(round(stop_value,3)) + '\n')
    if console:
        f.writelines('print yes\n')
    else:
        f.writelines('print no\n')

    if active:
        f.writelines('active? yes quit quit quit quit\n')
    else:
        f.writelines('active? no quit quit quit quit\n')

def modify_convergence_initial_iterations(f, criterion_name, iters_to_ignore):
  f.writelines('solve/convergence-conditions/conv-reports/edit ' + str(criterion_name) + '\n')
  f.writelines('initial-values-to-ignore ' + str(iters_to_ignore) + ' quit quit quit quit\n')

def toggle_convergence_criterions(f, orig_criterion_name, new_criterion_name):
  f.writelines('solve/convergence-conditions/conv-reports/edit ' + str(orig_criterion_name) + '\n')
  f.writelines('active? no quit quit quit quit\n')
  f.writelines('solve/convergence-conditions/conv-reports/edit ' + str(new_criterion_name) + '\n')
  f.writelines('active? yes\n')
  f.writelines('print? yes quit quit quit quit\n')

def modify_residual_criterions(f, cont, x_vel, y_vel, z_vel, energy, k, omega):
  f.writelines('solve/monitors/residual/convergence-criteria ' + str(cont) + ' ' + str(x_vel) + ' ' + str(y_vel) + ' ' + str(z_vel) + ' ' + str(energy) + ' ' + str(k) + ' ' + str(omega) + '\n')

def modify_residual_criterions_2D(f, cont, x_vel, y_vel, energy, k, omega):
  f.writelines('solve/monitors/residual/convergence-criteria ' + str(cont) + ' ' + str(x_vel) + ' ' + str(y_vel) + ' ' + str(energy) + ' ' + str(k) + ' ' + str(omega) + '\n')

def modify_convergence_criterion_iterations(f, criterion_name, previous_iterations):
  f.writelines('solve/convergence-conditions/conv-reports/edit ' + str(criterion_name) + '\n')
  f.writelines('previous-values-to-consider ' + str(previous_iterations) + ' quit quit quit quit\n')

def make_vis_plane(f,plane_name,normal=[1,0,0],distance_from_origin_m = 0):
    f.writelines('surface/plane-slice ' + plane_name + ' ' + str(normal[0]) + ' ' + str(normal[1]) + ' ' + str(normal[2]) + ' ' + str(distance_from_origin_m) + '\n')

def active_graphics_object(f,object_name):
    f.writelines('display/objects/display ' + object_name + '\n')

def save_picture(f,file_name):
    f.writelines('display/save-picture ' + file_name + '\n')

def make_point(f,point_name,point_location=[0,0,0]):
    f.writelines('surface/point-surface ' + point_name + ' ' + str(point_location[0]) + ' ' + str(point_location[1]) + ' ' + str(point_location[2]) + '\n')

def make_graphics_object(f,object_type,object_name,display_state_name,field,surface_list):
    surface_str = ''
    for i_surface in range(len(surface_list)):
        surface_str+=surface_list[i_surface] + ' '

    f.writelines('display/objects/create ' + object_type + ' ' + object_name + ' field ' + field + ' ' + ' display-state-name ' + display_state_name + ' surfaces-list ' + surface_str + '() range-option auto-range-on global-range no ' + 'quit quit\n')

def edit_graphics_object(f,object_type,object_name,display_state_name=None,field=None,surface_list=None):

    edit_str = 'display/objects/edit ' + object_name

    if field != None:
        edit_str+=' field ' + field + ' '
    if display_state_name != None:
        edit_str+=' display-state-name ' + display_state_name + ' '
    if surface_list != None:
        surface_str = ''
        for i_surface in range(len(surface_list)):
            surface_str+=surface_list[i_surface] + ' () '

        edit_str+= ' surfaces-list ' + surface_str + ' '

    f.writelines(edit_str+'quit quit\n')


def create_centerline_graphics(f, save_name):
  ## create the contours
  f.writelines('/preferences/graphics/colormap-settings/number-format-precision 2\n')
  f.writelines('/preferences/graphics/colormap-settings/number-format-type float\n')
  f.writelines('/preferences/graphics/colormap-settings/text-font-automatic-vertical-size 0.125\n')
  f.writelines('/preferences/graphics/colormap-settings/text-font-name \"Helvetica-Bold\"\n')
  f.writelines('/preferences/graphics/lighting/headlight \"On\"\n')
  f.writelines('display/views/read-views \"cfd_view\"\n')
  f.writelines('restore-view \"cfd_view\" quit quit\n')

  contours = ['mach_number', 'static_pressure', 'total_pressure', 'static_temp', 'total_temp', 'vel_mag', 'x_vel_mag', 'y_vel_mag', 'z_vel_mag', 'schlieren']
  for i_contour in contours:
    match i_contour:
      case 'mach_number':
        object_name = 'mach_contour'
        field = 'mach-number'
      case 'static_pressure':
        object_name = 'static_pressure_contour'
        field = 'pressure'
      case 'total_pressure':
        object_name = 'total_pressure_contour'
        field = 'total-pressure'
      case 'static_temp':
        object_name = 'static_temperature_contour'
        field = 'temperature'
      case 'total_temp':
        object_name = 'total_temperature_contour'
        field = 'total-temperature'
      case 'vel_mag':
        object_name = 'velocity_magnitude_contour'
        field = 'velocity-magnitude'
      case 'x_vel_mag':
        object_name = 'x_velocity_magnitude_contour'
        field = 'x-velocity'
      case 'y_vel_mag':
        object_name = 'y_velocity_magnitude_contour'
        field = 'y-velocity'
      case 'z_vel_mag':
        object_name = 'z_velocity_magnitude_contour'
        field = 'z-velocity'
      case 'schlieren':
        object_name = 'schlieren_contour'
        field = 'schlieren' ##PCM: when doing schlieren contours, we will have to update the colorbar to grey-reverse and set the limits 0-50 and turn off the colorbar

    f.writelines('display/objects/create contour ' + object_name + ' surfaces-list symmetry* () field ' + field + ' ' + 'range-option auto-range-on global-range no ' + 'quit quit\n')
    f.writelines('display/objects/edit/ ' + object_name + ' color-map/title-elements/ \"Variable Only\" quit quit\n')
    f.writelines('display/views/read-views \"cfd_view\"\n')
    f.writelines('display/views/restore-view cfd_view\n')
    f.writelines('display/objects/display ' + object_name + '\n')
    file_name = '\"' + save_name + '_' + object_name + '\"'
    save_picture(f, file_name)


#hooks for axi
#hooks for SST
#change wall to mass flow inlet
#set mass flow inlet including species
#change wall to heat source
#set heat flux
#AMR


def create_residuals_plot_2D(f, case_name, continuity, x_vel, y_vel, energy, k, omega, user_defined = None):
  f.writelines('plot/residuals ' + continuity + ' ' + x_vel + ' ' + y_vel + ' ' + energy + ' ' + k + ' ' + omega + '\n')
  saveName = case_name + '_residuals.png'
  f.writelines('display/save-picture ' + saveName + 'quit quit\n')


def create_residuals_plot(f, case_name, continuity, x_vel, y_vel, z_vel, energy, k, omega):

  f.writelines('plot/residuals ' + continuity + ' ' + x_vel + ' ' + y_vel + ' ' + z_vel + ' ' + energy + ' ' + k + ' ' + omega + '\n')
  saveName = case_name + '_residuals.png'
  f.writelines('display/save-picture ' + saveName + 'quit quit\n')


def create_angled_plane(f, degree, pt1, pt2, pt3):
  f.writelines('surface/plane-surface/ plane_' + str(degree) + '_deg three-points ' + str(pt1[0]) + ' ' + str(pt1[1]) + ' ' + str(pt1[2]) + ' ' + str(pt2[0]) + ' ' + str(pt2[1]) + ' ' + str(pt2[2]) + ' ' + str(pt3[0]) + ' ' + str(pt3[1]) + ' ' + str(pt3[2]) + ' no quit\n')

def create_symmetry_contours(f, field):
  f.writelines('display/objects/create/contour ' + str(field) +'\n')
  f.writelines('surfaces-list symmetry*\n\n')
  f.writelines('field ' + str(field) + '\n')
  f.writelines('coloring smooth\n')
  f.writelines('color-map/title-elements \"Variable Only" quit\n')
  f.writelines('range-option auto-range-on\n')
  f.writelines('global-range? yes quit quit quit quit\n')

def create_station_contours(f, station, field):
  f.writelines('display/objects/create/contour ' + str(station) + '_' + str(field) +'\n')
  f.writelines('surfaces-list ' + str(station) + '\n\n')
  f.writelines('field ' + str(field) + '\n')
  f.writelines('coloring smooth\n')
  f.writelines('color-map/title-elements \"Variable Only" quit\n')
  f.writelines('range-option auto-range-on\n')
  f.writelines('global-range? yes quit quit quit quit\n')

def create_plane_contours(f, plane, field):
  f.writelines('display/objects/create/contour ' + str(plane) + '_' + str(field) +'\n')
  f.writelines('surfaces-list ' + str(plane) + '\n\n')
  f.writelines('field ' + str(field) + '\n')
  f.writelines('coloring smooth\n')
  f.writelines('color-map/title-elements \"Variable Only" quit\n')
  f.writelines('range-option auto-range-on\n')
  f.writelines('global-range? yes quit quit quit quit\n')
#functions needed


def create_scene(f, planes, field):
  f.write('display/objects/create scene no_struts_' + str(field) + '\n')
  for i_plane in range(len(planes)):
    if i_plane == 0:
      f.writelines('graphics-objects add \"' + str(planes[i_plane]) + '\"\n')
      f.writelines('quit\n')
    else:
      f.writelines('add\n')
      f.writelines('\"' + str(planes[i_plane]) + '\"\n')
      f.writelines('quit\n')
  f.writelines('quit quit\n')
  f.writelines('display/objects/display no_struts_' + str(field) + '\n')


#write residuals /plot/residuals-set plot to file then plot/residuals yes yes yes whatever you want output

#turn residual monitor off/on for back pressuring cases

#start by loading case, defining bcs
#run a solution all the way out with fmg initialization fuel/ignitors off
#set fuel/ignitor bcs
#start with low fuel flow rate, 2000 iterations, write case/data
#increase fuel flow up to some amount, writing cases
#then sweep back through with AMR, writing out AMR after 2k iterations