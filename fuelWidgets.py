# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 07:24:57 2024

@author: mfund
"""

import numpy as np
from scipy.interpolate import interp1d,interpn

# datathief from Lee handbook 

I1_viscosityCompensationFactor = np.concatenate((np.arange(0.01,0.1,0.01),np.arange(0.1,1.1,0.1)),axis=0)
I2_pressureDifferential_kPa = [0.7,6.9,41,172,690,2760,6900,20700]
D_kinematicViscosity_cSt = np.zeros((len(I1_viscosityCompensationFactor),len(I2_pressureDifferential_kPa)))

D_kinematicViscosity_cSt[:,0] = [101.16028, 50.28087, 33.65535, 25.16527, 20.30557, 16.72811, 14.26621, 12.50819, 11.19691, 10.0231, 5.05134, 3.33462, 2.49341, 1.97056, 1.66896, 1.42333, 1.23933, 1.10175, 0.9931]
D_kinematicViscosity_cSt[:,1] = [457.42578, 202.1209, 122.79434, 88.69432, 67.24427, 54.25864, 45.00975, 38.65221, 33.42321, 29.71296, 13.03858, 7.97634, 5.19315, 3.69944, 2.74712, 2.14122, 1.64601, 1.2566, 1]
D_kinematicViscosity_cSt[:,2]= [1000, 870.72248, 688.1377, 562.99042, 473.53349, 415.18005, 364.01748, 330.39851, 306.17656, 274.07899, 142.99163, 90.55529, 59.36727, 37.33742, 23.48235, 13.97303, 7.97634, 3.99211, 1]
D_kinematicViscosity_cSt[:,3] = [993.1023, 1000, 993.1023, 1013.93946, 1000, 907.64461, 823.81873, 742.57698, 688.1377, 628.92259, 351.63507, 233.74239, 157.54143, 106.91979, 69.13204, 40.85286, 22.06416, 9.22422, 1]
D_kinematicViscosity_cSt[:,4] = [1000, 1006.94561, 993.1023, 1006.94561, 1000, 1006.94561, 993.1023, 1000, 993.1023, 993.1023, 907.64461, 599.17697, 415.18005, 285.70104, 189.91406, 113.79213, 66.78044, 32.06358, 1.00695]
D_kinematicViscosity_cSt[:,5] = [1000, 1013.93946, 1000, 993.1023, 1006.94561, 1006.94561, 993.1023, 993.1023, 986.25218, 1006.94561, 1013.93946, 1013.93946, 829.54066, 603.33862, 398.29092, 255.74999, 147.00587, 71.0728, 1]
D_kinematicViscosity_cSt[:,6] = [1000, 1000, 1000, 993.1023, 1000, 1013.93946, 993.1023, 1000, 1000, 993.1023, 1000, 993.1023, 1000, 993.1023, 993.1023, 628.92259, 361.5066, 174.77728, 1]
D_kinematicViscosity_cSt[:,7] = [1000, 1006.94561, 986.25218, 1013.93946, 993.1023, 1013.93946, 1006.94561, 993.1023, 1000, 1000, 993.1023, 993.1023, 1006.94561, 993.1023, 1013.93946, 1000, 993.1023, 467.02344, 1.00695]
    
points = (I1_viscosityCompensationFactor,I2_pressureDifferential_kPa)

def compute_Lohms(I_m3ps,H_Pa,nu_m2ps=4e-6,density_kgpm3=875): 
    #I is volumetric fuel flow rate
    #H is differential pressure pt_meter_US-pt_meter_DS in Pa (will know p2,p0 = f(t), target p1)
    
    #L = KV/I*sqrt(H/S) 
    # nu, S, V are actually lookups as function of T_fuel (and H_Pa for V)
    K = 28.8 #depends on pressure units and flow units. Here is for kPa,L/min
    nu_cSt = nu_m2ps*1e6 #kinematic viscosity from page O14, diesel at 60F
    S = density_kgpm3/1000 #specific gravity, diesel at 60F
        
    for i_NR in range(100):
        if i_NR==0:
            factorGuess = 0.9
        
        point = (factorGuess,H_Pa/1000)
        f = interpn(points,D_kinematicViscosity_cSt,point)-nu_cSt
        
        if abs(f) < 1e-2:
            V = factorGuess
            break
        
        factorGuess += 0.01
        point = (factorGuess,H_Pa/1000)
        df = (interpn(points,D_kinematicViscosity_cSt,point)-nu_cSt)-f
        
        fprime = df/0.01
        
        factorGuess+=-0.01
        factorGuess+=-f/fprime
        
    L_Lohms = K*V/(I_m3ps*1000*60)*np.sqrt((H_Pa/1000)/S)
    
    return L_Lohms

def compute_flow_rate(L_lohms,H_Pa,nu_m2ps=4e-6,density_kgpm3=875):
    
    K = 28.8 #depends on pressure units and flow units. Here is for kPa,L/min
    nu_cSt = nu_m2ps*1e6 #kinematic viscosity from page O14, diesel at 60F
    S = density_kgpm3/1000 #specific gravity, diesel at 60F
        
    for i_NR in range(100):
        if i_NR==0:
            factorGuess = 0.9
        
        point = (factorGuess,H_Pa/1000)
        f = interpn(points,D_kinematicViscosity_cSt,point)-nu_cSt
        
        if abs(f) < 1e-2:
            V = factorGuess
            break
        
        factorGuess += 0.01
        point = (factorGuess,H_Pa/1000)
        df = (interpn(points,D_kinematicViscosity_cSt,point)-nu_cSt)-f
        
        fprime = df/0.01
        
        factorGuess+=-0.01
        factorGuess+=-f/fprime
                
    I_m3ps = K*V/(L_lohms*1000*60)*np.sqrt((H_Pa/1000)/S)
        
    return I_m3ps

def compute_pt_drop_Pa(L_lohms,I_m3ps,nu_m2ps=4e-6,density_kgpm3=875):
    K = 28.8 #depends on pressure units and flow units. Here is for kPa,L/min
    nu_cSt = nu_m2ps*1e6 #kinematic viscosity from page O14, diesel at 60F
    S = density_kgpm3/1000 #specific gravity, diesel at 60F
    
    try:
        
        for i_NR_outer in range(100):
            
            V = np.nan
            
            if i_NR_outer==0:
                pressureDifferentialGuess_kPa = 50
                
            #here need to find the V that gives right centiStokes for this dP
            
            for i_NR_inner in range(100):
                
                if i_NR_inner == 0:
                    factorGuess = 0.9
            
                point = (factorGuess,pressureDifferentialGuess_kPa)
                f_inner = interpn(points,D_kinematicViscosity_cSt,point,bounds_error=False,fill_value=None)-nu_cSt
                #print('inner loop', point, f_inner)
    
                if abs(f_inner) < 1e-2:
                    V = factorGuess
                    #print('Inner loop converged')
                    break
                
                factorGuess += -0.01
                point = (factorGuess,pressureDifferentialGuess_kPa)
                df_inner = (interpn(points,D_kinematicViscosity_cSt,point,bounds_error=False,fill_value=None)-nu_cSt)-f_inner
                
                fprime_inner = df_inner/-0.01
                
                factorGuess+=0.01
                factorGuess+=-f_inner/fprime_inner*0.25
                
                if factorGuess < 0.02:
                    factorGuess = 0.02
                if factorGuess > 1.0:
                    factorGuess = 1.0
                
            f_outer = pressureDifferentialGuess_kPa*1000 - ((I_m3ps*1000*60*L_lohms)/(K*V))**2*S*1000
            #print('outer loop',pressureDifferentialGuess_kPa,f_outer)
            
            if abs(f_outer) < 1e-2:
                H_Pa = pressureDifferentialGuess_kPa*1000
                break
            
            pressureDifferentialGuess_kPa += 50
            V = np.nan
    
            for i_NR_inner in range(100):
                
                if i_NR_inner == 0:
                    factorGuess = 0.9
            
                point = (factorGuess,pressureDifferentialGuess_kPa)
                f_inner = interpn(points,D_kinematicViscosity_cSt,point,bounds_error=False,fill_value=None)-nu_cSt
                #print('inner loop', point, f_inner)
    
                if abs(f_inner) < 1e-2:
                    V = factorGuess
                    #print('Inner loop converged')
                    break
                
                factorGuess += -0.01
                point = (factorGuess,pressureDifferentialGuess_kPa)
                df_inner = (interpn(points,D_kinematicViscosity_cSt,point,bounds_error=False,fill_value=None)-nu_cSt)-f_inner
                
                fprime_inner = df_inner/-0.01
                
                factorGuess+=0.01
                factorGuess+=-f_inner/fprime_inner
                
                if factorGuess < 0.02:
                    factorGuess = 0.02
                if factorGuess > 1.0:
                    factorGuess = 1.0
                
            df_outer = pressureDifferentialGuess_kPa*1000 - ((I_m3ps*1000*60*L_lohms)/(K*V))**2*S*1000 - f_outer
              
            fprime_outer = df_outer/50
            
            pressureDifferentialGuess_kPa+=-50
            pressureDifferentialGuess_kPa+=-f_outer/fprime_outer*0.25
            
            # if pressureDifferentialGuess_kPa > 20650:
            #     pressureDifferentialGuess_kPa = 20650
            
        return H_Pa, V
    except Exception as error:
        print(error,pressureDifferentialGuess_kPa,factorGuess)
    
def plain_orifice_atomizer(D_m,L_m,r_m,rho_l_kgpm3,mu_kgpms,sigma_kgps2,p_vapor_Pa,mdot_fuel_kgps,pinput_Pa,force_single_phase=False,input_mode='p2',Cd_factor=1):
    #D_m is orifice physical diameter
    #L_m is orifice length
    #r_m is orifice inlet radius of cruvature
    #rho is liquid fuel density
    #mu is "" dynamic viscosity
    #p_vapor is "" vapor pressure
    #p2 is orifice exit gas pressure
    
    Cct = 0.611
    CA = 3+L_m/3.6*D_m #suggested by Reitz from Fluent page 
    phi_start_deg = 0
    phi_stop_deg = 360

    
    dp_Pa_guess = 2.5e4
    A_m2 = np.pi*D_m**2/4

    for i_NR in range(10000): #determine the nozzle state
    
        
        if input_mode=='p2':
            p2_Pa = pinput_Pa
            p1_Pa = pinput_Pa+dp_Pa_guess
        elif input_mode=='p1':
            p1_Pa = pinput_Pa
            p2_Pa = p1_Pa-dp_Pa_guess
        
        #these calcs should always be the same
        
        K = (p1_Pa-p_vapor_Pa)/(p1_Pa-p2_Pa)
        Re_h = D_m*rho_l_kgpm3/mu_kgpms*np.sqrt(2*(dp_Pa_guess)/rho_l_kgpm3)
        
        K_incep = 1.9*(1-r_m/D_m)**2-1000/Re_h
        
        if force_single_phase:
            K_incep = 1 #you'll never get here
        
        K_crit = 1+1/((1+L_m/(4*D_m))*(1+2000/Re_h)*np.exp(70*r_m/D_m))
        
        #single phase calcs
        
        nozzle_state = ''
        
        if K > K_incep: 
            if K >= K_crit: #single phase
                Cdu = 0.827-0.0085*L_m/D_m
                Cd = 1/(1/Cdu+20*(1+2.25*L_m/D_m)/Re_h)
                nozzle_state = 'single phase'
            elif K < K_crit: #flipped
                Cd = Cct
                nozzle_state = 'flipped'
        
        elif K <= K_incep: 
            if K >= K_crit: #cavitating
                Cc = 1/np.sqrt(1/(Cct**2)-11.4*r_m/D_m)
                if Cc > 1: #jet area can't increase due to blunt lip
                    Cc = 1
                Cd = Cc*np.sqrt(K)
                if Cd > 1:
                    Cd = 1 #jet area can't increase                
                nozzle_state = 'cavitating'
                
            elif K < K_crit: #flipped
                Cd = Cct
                nozzle_state = 'flipped'
        
        #else:
        #    print('What have you done?')
        
        Cd=Cd*Cd_factor
            
        mdot_eff_kgps = 2*np.pi*mdot_fuel_kgps/np.radians(phi_stop_deg-phi_start_deg)
        dp_actual_Pa = (mdot_eff_kgps/(Cd*A_m2))**2/(2*rho_l_kgpm3) #now iterate in NR until closed 
        
        f = dp_actual_Pa-dp_Pa_guess
        #print(i_NR,f,dp_Pa_guess,nozzle_state,K,K_incep,Cd)
        if abs(f) < 100: #was 100
            dp_converged_Pa = dp_Pa_guess
            nozzle_state_converged = nozzle_state
            #print(K,K_incep,K_crit)
            break
        
        dp_Pa_guess += 2.5e4
        p1_Pa = p2_Pa+dp_Pa_guess
        
        #these calcs should always be the same
        
        K = (p1_Pa-p_vapor_Pa)/(p1_Pa-p2_Pa)
        Re_h = D_m*rho_l_kgpm3/mu_kgpms*np.sqrt(2*(dp_Pa_guess)/rho_l_kgpm3)
        
        K_incep = 1.9*(1-r_m/D_m)**2-1000/Re_h
        K_crit = 1+1/((1+L_m/(4*D_m))*(1+2000/Re_h)*np.exp(70*r_m/D_m))
        
        #single phase calcs
        
        nozzle_state = ''
        
        if K > K_incep: 
            if K >= K_crit: #single phase
                Cdu = 0.827-0.0085*L_m/D_m
                Cd = 1/(1/Cdu+20*(1+2.25*L_m/D_m)/Re_h)
                nozzle_state = 'single phase'
            elif K < K_crit: #flipped
                Cd = Cct
                nozzle_state = 'flipped'
        
        elif K <= K_incep: #flipped
            if K >= K_crit: #cavitating
                Cc = 1/np.sqrt(1/(Cct**2)-11.4*r_m/D_m)
                Cd = Cc*np.sqrt(K)
                if Cd > 1:
                    Cd = 1 #jet area can't increase
                nozzle_state = 'cavitating'
                #print(K,nozzle_state,Cc,Cd,Cct)
                
            elif K < K_crit: #flipped
                Cd = Cct
                nozzle_state = 'flipped'
        
        # else:
        #     print('What have you done?')
        
        Cd=Cd*Cd_factor
            
        mdot_eff_kgps = 2*np.pi*mdot_fuel_kgps/np.radians(phi_stop_deg-phi_start_deg)
        dp_actual_Pa = (mdot_eff_kgps/(Cd*A_m2))**2/(2*rho_l_kgpm3) #now iterate in NR until closed 
        
        df = (dp_actual_Pa-dp_Pa_guess)-f
        fprime = df/2.5e4
        
        dp_Pa_guess+=-2.5e4
        dp_Pa_guess+=-f/fprime*0.25
        #print('fuelWidget',K,K_incep,K_crit)
        
    #could compute spray angle here also
    
    if nozzle_state_converged == 'single phase':
        V_fuel_mps = mdot_eff_kgps/(rho_l_kgpm3*A_m2)#*Cd)
        lambdaa = D_m/8
        We = rho_l_kgpm3*V_fuel_mps**2*lambdaa/sigma_kgps2
        d_SMD_m = 133*lambdaa*We**(-0.74)
        s = 3.5
        d0_m = 1.2726*d_SMD_m*(1-1/s)**(1/s)
        
    elif nozzle_state_converged == 'cavitating':
        V_fuel_mps = (2*Cc*p1_Pa-p2_Pa+(1-2*Cc)*p_vapor_Pa)/(Cc*np.sqrt(2*rho_l_kgpm3*(p1_Pa-p_vapor_Pa)))
        lambdaa = np.sqrt(4*Cd*A_m2/np.pi)/8
        We = rho_l_kgpm3*V_fuel_mps**2*lambdaa/sigma_kgps2
        d_SMD_m = 133*lambdaa*We**(-0.74)
        s = 1.5
        d0_m = 1.2726*d_SMD_m*(1-1/s)**(1/s)
        
    elif nozzle_state_converged == 'flipped':
        V_fuel_mps = mdot_eff_kgps/(rho_l_kgpm3*Cct*A_m2)
        d0_m = D_m*np.sqrt(Cct)
    # else:
    #     print('What have you done?')
        
    return p1_Pa,p2_Pa,dp_converged_Pa,nozzle_state_converged,mdot_eff_kgps,Cd,V_fuel_mps,d0_m,K,K_incep,K_crit

class fuelPropertyLibrary:
    def __init__(self):
        self.fuelTypes = ['diesel','water','jet_A1']
        
    def getProperties(self,fuelType,temperature):
        
        if fuelType in self.fuelTypes:
            #have density, viscosity, surface tension as functions of temp
            
            fuelProperties = {}
            
            if fuelType == 'diesel':

                #kinematic viscosity from Lee isn't datathiefable, using the equation from Razzaq instead: ln(eta) = -f*T+y where for diesel f = 0.018 and Y = 2.092
                fuelProperties['kinematicViscosity_m2ps'] = np.exp((temperature-273.15)*-0.018+2.092)*1e-6
                
                temperature_F = [-22.52707581227437, 0.21660649819494182, 20.361010830324915, 40.07220216606498, 60.64981949458483, 80.57761732851986, 100.28880866425993, 119.78339350180505,140.04689331770223, 160.09378663540446, 180.14067995310668, 199.97655334114887, 219.60140679953105, 239.6483001172333, 260.1172332942556, 279.74208675263776, 300.21101992966, 319.4138335287222]
                temperature_K = np.multiply(np.array(temperature_F)+459.67,5/9) #from Lee handbook
                specificGravity = np.array([0.9196165191740413, 0.909882005899705, 0.8985250737463126, 0.8887905604719764, 0.8790560471976401, 0.8676991150442477, 0.8579646017699114, 0.8476892822025565,0.8374519230769231, 0.8284615384615385, 0.8173557692307691, 0.8067788461538461, 0.7962019230769231, 0.7850961538461538, 0.775576923076923, 0.7655288461538461, 0.754951923076923, 0.7438461538461538])
                specificGravityInterpolant = interp1d(temperature_K,specificGravity,fill_value='extrapolate')
                fuelProperties['density_kgpm3'] = 1000*specificGravityInterpolant(temperature)
                fuelProperties['dynamicViscosity_kgpms'] = fuelProperties['density_kgpm3']*fuelProperties['kinematicViscosity_m2ps']
                
                temperature_K = np.array([293.89379, 323.92861, 373.40627, 423.46431, 448.13059]) #from https://doi.org/10.1016/j.fuel.2012.05.006
                surfaceTension_Npm = np.multiply([25.74826, 23.10905, 20.17981, 18.03364, 15.68445],1/1000)
                surfaceTensionInterpolant = interp1d(temperature_K,surfaceTension_Npm,fill_value='extrapolate')
                fuelProperties['surfaceTension_Npm'] = surfaceTensionInterpolant(temperature)   
                fuelProperties['vaporPressure_Pa'] = 300
                
                return fuelProperties
            
            if fuelType == 'water':
                
                #all quantities via interpolant (eng Toolbox)
                
                temperature_K_density = [272.57573, 297.26945, 322.53744, 348.09257, 373.07343, 398.05429, 422.74801, 447.72887, 472.99686, 498.26485, 522.95858, 547.93943, 572.34602, 597.32688, 612.25796, 623.16914, 633.21891]
                density_kgpm3 = [1000.55494, 998.11321, 987.12542, 974.91676, 958.43507, 938.90122, 917.53607, 892.50832, 864.42841, 832.6859, 798.50166, 757.60266, 713.6515, 655.04994, 612.93008, 569.58935, 526.85905]
                fuelProperties['density_kgpm3'] = interp1d(temperature_K_density,density_kgpm3,fill_value='extrapolate')(temperature)
                
                temperature_K_mu = [273.39776, 277.36198, 282.06948, 286.28145, 293.46659, 303.12935, 318.24291, 332.36542, 345.49687, 359.37161, 372.2553, 385.88228, 402.48242, 422.05571, 447.0798, 472.84718, 500.3489, 522.89535, 550.14931, 572.20024, 598.95867, 622.99171, 632.90224]
                dynamic_viscosity_kgpm3 = [0.00179, 0.0016, 0.00141, 0.0012, 0.001, 0.0008, 0.0006, 0.00047, 0.0004, 0.00033, 0.00028, 0.00025, 0.00021, 0.00018, 0.00015, 0.00013, 0.00012, 0.0001, 0.00009, 0.00008, 0.00008, 0.00007, 0.00006]
                fuelProperties['dynamicViscosity_kgpms'] = interp1d(temperature_K_mu,dynamic_viscosity_kgpm3,fill_value='extrapolate')(temperature)
                fuelProperties['kinematicViscosity_m2ps'] = fuelProperties['dynamicViscosity_kgpms']/fuelProperties['density_kgpm3']

                temperature_K_sigma = [273.15, 323.81667, 373.27121, 423.21061, 473.15, 522.60455, 572.54394, 622.72576, 647.69545]
                surfaceTension_Npm = [0.07577, 0.06781, 0.05885, 0.04802, 0.03782, 0.0265, 0.01481, 0.00386, 0]
                fuelProperties['surfaceTension_Npm'] = interp1d(temperature_K_sigma,surfaceTension_Npm,fill_value='extrapolate')(temperature)

                temperature_K_vapor_pressure = [272.42985, 289.95354, 308.19738, 319.47975, 337.72359, 358.60799, 373.01102, 397.97628, 423.18159, 446.46649, 473.1121, 497.59725, 523.28266, 548.72802, 573.45322, 600.33888, 623.38373, 643.54798]
                vapor_pressure_Pa = np.array([0.60607, 2.02253, 5.79711, 10.28478, 25.02499, 59.48264, 99.53309, 239.36742, 482.99551, 866.97811, 1593.07335, 2484.99319, 4109.80426, 6117.66041, 8792.36464, 12636.47715, 16733.13195, 21393.63357])*1000
                fuelProperties['vaporPressure_Pa'] =  interp1d(temperature_K_vapor_pressure,vapor_pressure_Pa,fill_value='extrapolate')(temperature)
                
                return fuelProperties
            
            if fuelType == 'jet_A1':
                
                #density per CRC report 530 pg 25
                temperature_K_density = np.array([-40.3709, -35.17831, -30.54208, -25.16405, -20.15692, -15.14978, -10.32811, -5.50641, -0.31384, 4.69329, 9.70042, 14.89301, 19.90014, 24.90728, 28.80171, 34.36519, 39.92867, 44.93581, 49.75748, 54.76462, 59.9572, 64.77888, 69.78601, 74.60771, 79.80028, 84.99287, 90])+273.15
                density_kgpm3 = [849.32629, 845.2576, 842.29855, 837.85998, 834.90093, 831.572, 827.87319, 824.17438, 820.47557, 817.14664, 813.44783, 810.1189, 806.42008, 802.35139, 799.39235, 795.69353, 791.99472, 788.66579, 784.59709, 781.26816, 777.56936, 774.24043, 770.54162, 767.58257, 763.88376, 759.81507, 756.85601]
                density_interpolant = interp1d(temperature_K_density,density_kgpm3,fill_value='extrapolate')
                fuelProperties['density_kgpm3'] = density_interpolant(temperature)
                
                #viscosity per "" pg 34
                temperature_K_nu = [223.154943440255, 232.904058430212, 243.057175381875, 252.385425989265, 262.584469024998, 272.922000180299, 282.824976557966, 293.088837302231, 302.82659382285, 313.511247753011, 323.291794514894, 333.049327718686, 342.762929794964, 354.162336905824, 363.409984852314, 373.268002990598, 383.391465364396, 392.617808645073, 403.658611510332, 412.955275261538, 417.475202889444]
                kinematic_viscosity_m2ps = [1.52245393496157E-05, 8.88124043941142E-06, 5.75188196691694E-06, 4.13187947466026E-06, 3.02924332143278E-06, 2.32359552433833E-06, 1.88887466166847E-06, 1.5309716663995E-06, 1.31185593829369E-06, 1.11648234934541E-06, 9.85814435461779E-07, 8.73256182215237E-07, 7.87957843915543E-07, 7.07259105456088E-07, 6.46221556525333E-07, 5.85417384044032E-07, 5.33813141844797E-07, 4.83193148887555E-07, 4.35341296461287E-07, 3.99545366751347E-07, 3.82146966410891E-07]
                nu_interpolant = interp1d(temperature_K_nu,kinematic_viscosity_m2ps,fill_value='extrapolate')
                
                fuelProperties['kinematicViscosity_m2ps'] = nu_interpolant(temperature)
                fuelProperties['dynamicViscosity_kgpms'] = fuelProperties['kinematicViscosity_m2ps']*fuelProperties['density_kgpm3'] 
                
                #surface tension per "" 38
                temperature_K_sigma = np.array([-38.52112, -28.80281, -18.66196, -8.94366, 0.77465, 10.9155, 20.84507, 32.88732, 45.77465, 60.56339, 74.92958, 80])+273.15
                surfaceTension_Npm = np.array([28.04255, 27.25532, 26.46808, 25.6383, 24.87234, 24.10638, 23.2766, 22.31915, 21.31915, 20.14894, 19, 18.59575])/1000
                fuelProperties['surfaceTension_Npm'] = interp1d(temperature_K_sigma,surfaceTension_Npm,fill_value='extrapolate')(temperature)
                
                #vapr pressures per "" pg 46
                temperature_K_vapor_pressure = 1/(np.array([3.09727, 3.00096, 2.91814, 2.83917, 2.76019, 2.68122, 2.61766, 2.54831, 2.48282, 2.40963])/1000)
                vapor_pressure_Pa = np.array([1.31063, 2.01245, 2.93122, 4.13091, 5.89894, 8.76382, 11.56212, 16.08065, 21.07571, 29.70157])*1000
                fuelProperties['vaporPressure_Pa'] =  interp1d(temperature_K_vapor_pressure,vapor_pressure_Pa,fill_value='extrapolate',kind='linear')(temperature)
                
                if fuelProperties['vaporPressure_Pa'] < 0:
                    fuelProperties['vaporPressure_Pa'] = 0 
                
                return fuelProperties

        else:
            print('Fuel type not in library')
        
