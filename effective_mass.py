# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import random
import matplotlib.pyplot as plt
import numpy as np
import pylab
from matplotlib import rc
from numpy import diff



z_reion=10

with open('output_xe_new.txt') as f:
    #Read the given file
    data=[[float(x) for x in line.split()] for line in f]

    redshift=['']
    ioniz=['']
for num in range(0,len(data)):
    z,f=data[num]
    if z > z_reion:
            redshift.append(z)
            ioniz.append(f)            
    elif z_reion > z and z > z_reion-4:
        redshift.append(z)
        ioniz.append(1.07- 0.02*z)
    else:
        redshift.append(z)
        ioniz.append(1.07)
        
        
        
 
redshift.pop(0)
ioniz.pop(0)



def proton_density(z): 
    baryon_to_photon_ratio=6.7e-10
    T_CMB=2.7
    Y_p=0.25
    cm3_to_ev3_conversion=5.32e11
    pd=cm3_to_ev3_conversion*0.25*(1-Y_p/2)*baryon_to_photon_ratio*(T_CMB*(1+z))**3
    return pd

def plasma_mass_squared(ion_frac,proton_den):
    plasma_mass_squared= 1.79e-07*ion_frac*proton_den
    return plasma_mass_squared


def plasma_mass_squared_2(ion_frac,proton_den):
    plasma_mass_squared= 1.4e-21*ion_frac*proton_den*(5.24e11)
    return plasma_mass_squared


def neutral_atoms_mass_squared(ion_frac,proton_den,photon_freq,z):
    neutral_atom_mass_squared= -5.36e-12*(photon_freq*(1+z))**2*(1-ion_frac)*proton_den
    return neutral_atom_mass_squared

def mass_total(ion_frac,photon_freq,proton_den,z):
    mass_total=1.4e-21* (ion_frac - (1-ion_frac)*7.3e-3 * (photon_freq*(1+z))**2)*proton_den 
    return mass_total

#Compute with omega/Tcmb=1#
omega_over_TCMB=1
gamma_freq=omega_over_TCMB*2.73*8.61e-05

cosmo_mass=['']
cosmo_mass_squared_list=['']
for num in range(0, len(redshift)):
    cosmo_mass_squared=mass_total(ioniz[num], gamma_freq, proton_density(redshift[num]), redshift[num] )
    cosmo_mass_squared_list.append(cosmo_mass_squared)
    cosmo_mass.append(np.sqrt(cosmo_mass_squared))

cosmo_mass.pop(0)
cosmo_mass_squared_list.pop(0)


#now changing omega/Tcmb=3#

omega_over_TCMB=3
gamma_freq=omega_over_TCMB*2.73*8.61e-05

cosmo_mass_2=['']
cosmo_mass_squared_list_2=['']
for num in range(0, len(redshift)):
    cosmo_mass_squared=mass_total(ioniz[num], gamma_freq, proton_density(redshift[num]), redshift[num] )
    cosmo_mass_squared_list_2.append(cosmo_mass_squared)
    cosmo_mass_2.append(np.sqrt(cosmo_mass_squared))

cosmo_mass_2.pop(0)
cosmo_mass_squared_list_2.pop(0)


#now changing omega/Tcmb=4#

omega_over_TCMB=4
gamma_freq=omega_over_TCMB*2.73*8.61e-05

cosmo_mass_3=['']
cosmo_mass_squared_list_3=['']
for num in range(0, len(redshift)):
    cosmo_mass_squared=mass_total(ioniz[num], gamma_freq, proton_density(redshift[num]), redshift[num] )
    cosmo_mass_squared_list_3.append(cosmo_mass_squared)
    cosmo_mass_3.append(np.sqrt(cosmo_mass_squared))
    
cosmo_mass_3.pop(0)
cosmo_mass_squared_list_3.pop(0)


#now changing omega/Tcmb=10#

omega_over_TCMB=10
gamma_freq=omega_over_TCMB*2.73*8.61e-05
redshift_final=['']
cosmo_mass_4=['']
cosmo_mass_squared_list_4=['']
for num in range(0, len(redshift)):
    cosmo_mass_squared=mass_total(ioniz[num], gamma_freq, proton_density(redshift[num]), redshift[num] )
    if cosmo_mass_squared < 0:
        redshift_new= redshift[num] - 0.0078*redshift[num]
        cosmo_mass_squared=mass_total(ioniz[num], gamma_freq, proton_density(redshift[num]), redshift_new)
        redshift_final.append(redshift_new)
        if cosmo_mass_squared <0:
                cosmo_mass_squared=mass_total(ioniz[num], gamma_freq, proton_density(redshift[num]), redshift[num] )
                redshift_final.pop()
                redshift_final.append(redshift[num])
    else:
        redshift_final.append(redshift[num])    
    cosmo_mass_squared_list_4.append(cosmo_mass_squared)
    cosmo_mass_4.append(np.sqrt(cosmo_mass_squared))
    
    
cosmo_mass_4.pop(0)
cosmo_mass_squared_list_4.pop(0)
redshift_final.pop(0)



    


plt.plot(redshift,ioniz,'blue', label='Ionization fraction')
plt.yscale('log')
plt.xscale('log')


dxdy=diff(np.log(cosmo_mass_squared_list))/diff(redshift)


def log_derivative(redshift, dlndz):
    h=0.7
    h_to_eV=1.33e-32
    omega_lambda=0.7
    omega_rad= 0.0005
    omega_matter=0.3
    log_derivative=-dlndz* h*h_to_eV*(1+redshift)*np.sqrt(omega_lambda+ omega_matter*(1+redshift)**3 + omega_rad*(1+z)**4)
    
    return log_derivative

def probability(log_der, photon_frequency,chi_0,DP_mass):
    probability= np.pi*(DP_mass**2)*(chi_0**2)*(log_der)**(-1)*photon_frequency**(-1)
    return probability


        



effective_mass=plt.figure(2)
plt.plot(redshift,cosmo_mass,'black', label='Photon effective mass')
plt.plot(redshift,cosmo_mass_2,'red', label='Photon effective mass')
plt.plot(redshift,cosmo_mass_3,'blue', label='Photon effective mass')
plt.plot(redshift_final,cosmo_mass_4,'green', label='Photon effective mass')
plt.xlabel('Redshift ', fontsize=18)
plt.ylabel('$\omega/T$', fontsize=18)
plt.yscale('log')
plt.xscale('log')




redshift.pop(0)

mass_squared_derivative=plt.figure(3) #photon mass squared derivative wrt redshift#
plt.plot(redshift,dxdy,'green', label='Photon mass squared derivative')
plt.yscale('log')
plt.xscale('log')






