# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 12:33:20 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib.ticker import FormatStrFormatter
import datetime
t_Start = datetime.datetime.now()

# Constants
R = 8.314           # J/(mol*K)
cm_per_inch = 2.54  # cm per inch

# Temperature range
T_start = 300       # °C
T_end = 1200        # °C
temperatures_C = np.linspace(T_start, T_end, int((T_end-T_start)+1))
temperatures_K = temperatures_C + 273.15 # conversion to Kelvin

# Initial pressures
# unit shoult not matter since the values are relative, 
# so if p_standard is equal to p_total the result should always be the same...
# ...but it is not, I am not sure why. So use bar.
# pressure effect according to formula (32): https://link.springer.com/article/10.1007/s10973-012-2334-2
p_standard = 1  # bar
p_total = 1   # bar


# different Equilibrium constant functions
# Erste Ulich’sche Näherung, https://av.tib.eu/media/15706
d_R_H = 172530  # J/mol
d_R_S = 176.69  # J/(mol*K)
def K_eq_formula1(T):
    d_R_G = d_R_H - d_R_S * T
    return np.exp(-d_R_G / (R * T))

# Andrzej Mianowski et al., formula (34): https://link.springer.com/article/10.1007/s10973-012-2334-2
def K_eq_formula2(T):
    return np.exp(-20780.9 / T + 20.32)


# Function to calculate and plot equilibrium for a given K_eq function
def calculate_mole_fractions(K_eq_func, p_total, label, color):
    mole_fractions_CO = []
    mole_fractions_CO2 = []
    p_CO = []
    p_CO2 = []
    p_CO2_initial = p_total*0.5
    p_CO_initial = 0.0
    for T in temperatures_K:
        # Use fsolve to adjust p_CO2_initial until p_total_eq = p_total
        p_CO2_init_adjusted = fsolve(equilibrium_function, p_CO2_initial, args=(K_eq_func, T, p_total, p_CO_initial))[0]
        # Calculate equilibrium pressures with the adjusted p_CO2_initial
        K_eq = K_eq_func(T)
        K_eq_p_corrected = K_eq * (p_standard/p_total)
        x = 1/8 * (K_eq_p_corrected**0.5 * (8 * p_CO_initial + 16 * p_CO2_init_adjusted + K_eq_p_corrected)**0.5 - 4 * p_CO_initial - K_eq_p_corrected)
        p_CO2_eq = p_CO2_init_adjusted - x
        p_CO_eq = p_CO_initial + 2 * x
        #print(p_CO2_eq+p_CO_eq, p_CO_eq, p_total-p_CO2_eq)
        p_CO.append(p_CO_eq)
        p_CO2.append(p_CO2_eq)
        # Mole fractions
        mole_fraction_CO = p_CO_eq / p_total
        mole_fractions_CO.append(mole_fraction_CO)
        mole_fraction_CO2 = p_CO2_eq / p_total
        mole_fractions_CO2.append(mole_fraction_CO2)
    
    if color == "none":
        ax1.plot(temperatures_C, mole_fractions_CO, label=label)
    else:
        ax1.plot(temperatures_C, mole_fractions_CO, label=label, color=color)



# isobaric calculation: We need to iterate to find the correct initial pressure 
# to end up with the sum of both partial pressures equal to p_total
# Function to first calculate equilibrium pressures and then the error to the target pressure to correct for changes in volume
def equilibrium_function(p_CO2_initial, K_eq_func, T, p_total, p_CO_initial):
    K_eq = K_eq_func(T)
    K_eq_p_corrected = K_eq * (p_standard/p_total)
    x = 1/8 * (K_eq_p_corrected**0.5 * (8 * p_CO_initial + 16 * p_CO2_initial + K_eq_p_corrected)**0.5 - 4 * p_CO_initial - K_eq_p_corrected)
    p_CO2_eq = p_CO2_initial - x
    p_CO_eq = p_CO_initial + 2 * x
    p_total_eq = p_CO2_eq + p_CO_eq
    return p_total_eq - p_total


# prepare Plot for the results
fig, ax1 = plt.subplots(dpi=200, figsize=(12/cm_per_inch, 9/cm_per_inch))

# actual calculation
calculate_mole_fractions(K_eq_formula1, p_total, "Erste Ulich’sche Näherung", "r")
calculate_mole_fractions(K_eq_formula2, p_total, "Andrzej Mianowski et al.", "b")


# Plot settings and annotation
ax1.set_ylabel('Mole fraction of CO', color='k')
ax1.tick_params(axis='y', colors='k')
ax1.tick_params(axis='x', colors='k')
ax1.xaxis.set_major_formatter(FormatStrFormatter('%g °C'))
ax1.set_ylim(ymin=0, ymax=1)
ax1.set_xlim(xmin=T_start, xmax=T_end)

# Secondary y-axis for CO2 mole fraction
ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim()[::-1])  # Reverse the y-axis
ax2.set_ylabel('Mole fraction of CO$_2$', color='k')
ax2.tick_params(axis='y', colors='k')

ax1.legend(loc='best');
ax1.grid(True)
plt.tight_layout()
plt.title(f'Boudouard Equilibrium, p = {p_total:.0f} bar', color='k')





# calculation and plot for different pressures
p_start = 0.01
p_end = 100
p_array = np.logspace(np.log10(p_start), np.log10(p_end), int((np.log10(p_end)-np.log10(p_start))+1), base=10)
# prepare Plot for the results
fig, ax1 = plt.subplots(dpi=300, figsize=(12/cm_per_inch, 9/cm_per_inch))

for p in p_array:
    if p >= 1:
        label = f"{p:.0f} bar"
    else:
        n_digits = int(-np.log10(p))
        label = f"{p:.{n_digits}f} bar"
    calculate_mole_fractions(K_eq_formula2, p, label, "none")
    
# Plot settings and annotation
ax1.set_ylabel('Mole fraction of CO', color='k')
ax1.tick_params(axis='y', colors='k')
ax1.tick_params(axis='x', colors='k')
ax1.xaxis.set_major_formatter(FormatStrFormatter('%g °C'))
ax1.set_ylim(ymin=0, ymax=1)
ax1.set_xlim(xmin=T_start, xmax=T_end)

# Secondary y-axis for CO2 mole fraction
ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim()[::-1])  # Reverse the y-axis
ax2.set_ylabel('Mole fraction of CO$_2$', color='k')
ax2.tick_params(axis='y', colors='k')

ax1.legend(loc='best');
ax1.grid(True)
plt.tight_layout()
plt.title('Boudouard Equilibrium at different pressures', color='k')


plt.show()

t_Ende = datetime.datetime.now()
print("Runtime total: ", t_Ende - t_Start) # about 1 second for me
