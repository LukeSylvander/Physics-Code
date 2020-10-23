#    Import Numpy
import numpy as np
import math

"""
     MeasureSimulation.py
---------------------------------------------------------------------------------------------------
     - Kinetic_Energy
     - Localised_Temperature
     - Localised_Pressure
"""

"""
     Function Kinetic_Energy
------------------------------------------------------------------------------------------------------
     - Calculates the kinetic energy of the simulation
          Inputs:
               - grid data
               - reduced mass
               - inertia
          Outputs:
               - Kinetic_Energy = T
"""

def potential_energy(r, epsilon, r_0):
        U = epsilon * ((r_0/np.float_(r))**12 - 2 * (r_0 / np.float_(r)) ** 6)
        return U

def sum_potential_energy(list_of_distant, epsilon, r_0):
        i = 0
        U = []
        #       Calculate mag distace
        mag_dis_1 = 0
        mag_dis_2 = 0
        #Set up while loop
        while len(list_of_distant[1][0])> i:
                #       print(list_of_distant[2][0][i][0])
                mag_dis_1 =  np.float_(math.sqrt(list_of_distant[1][0][i]**2 +
                                                 list_of_distant[1][1][i]**2))
                mag_dis_2 =  np.float_(math.sqrt(list_of_distant[2][0][i]**2 +
                                                 list_of_distant[2][1][i]**2))
                
                U.append(potential_energy(np.float_(mag_dis_1), epsilon, r_0) +
                         potential_energy(np.float_(mag_dis_2), epsilon, r_0))

                i += 1
        U_int = sum(U)
        return U_int


def unit_diatomic_kinetic_energy(velocity, omega, mass, a):
     n = 0
     T = 0
     m = mass
     mu = (np.float_(m) * np.float_(m))/(np.float_(m) + np.float_(m))
     I = mu * a ** 2
     while n != len(omega):
##          velocity_mag = math.sqrt(velocity[0][n]**2 + velocity[1][n]**2)
##          velocity_dir = math.atan2(velocity[1][n], velocity[0][n])
          angular_velocity = omega[n]

          v_x = velocity[0][n] 
          v_y = velocity[1][n]
         # print('x = ' + str(v_x) + ', y = ' + str(v_y) + ' omega = ' + str(omega)) 
          
          T = T + (1/2 * mass * (v_x ** 2 + v_y ** 2) + 1/2 * I * angular_velocity ** 2)
          n = n + 1
     return T

def diatomic_kinetic_energy(grid_data, mass, a):
     n = 0
     T = 0
     m = mass
     mu = (np.float_(m) * np.float_(m))/(np.float_(m) + np.float_(m))
     I = mu * a ** 2
     while n != len(grid_data):
          velocity_mag = grid_data[n][2][0]
          velocity_dir = grid_data[n][2][1]
          angular_velocity = grid_data[n][4]

          v_x = velocity_mag * math.cos(velocity_dir)
          v_y = velocity_mag * math.sin(velocity_dir)
          
          T = T + (1/2 * mass * (v_x ** 2 + v_y ** 2) + 1/2 * I * angular_velocity ** 2)
          n = n + 1
     return T
##def diatomic_kinetic_energy(linear_v, angular_velocity, mass, a):
##     i = 0
##     linear_velocity = []
##     while i != len(angular_velocity):
##          linear_velocity.appen((np.float_(linear_v[i][0]) ** 2 + np.float_(linear_v[i][1]) ** 2 ) ** (1/2))
##          i = i + 1
##     
##     n = 0
##     T = 0
##     interia = ((mass/2) * ((a/2) ** 2))
##     while n != len(angular_velocity):
##          T = T + (1/2 * mass * linear_velocity[n] ** 2 + 1/2 * interia * angular_velocity[n] ** 2)
##          n = n + 1
##     return T

"""
     Function localised temperature
---------------------------------------------------------------------------------------------------------
     - Calculates the localised temp
          Inputs:
               - kinetic energy
          Outputs:
               - Localised Temperature
"""

def localised_temperature(kinetic_energy):
     boltz_con = 1.380649e-23
     T = (2 / (3 * boltz_con)) * kinetic_energy
     return T

##def localised_temperature(grid_data, mass, a):
##     n = 0
##     E = 0
##     boltz_con = 1.380649e-23
##     interia = ((mass/2) * ((a/2) ** 2))  * 2 
##     while n != len(grid_data):
##          velocity_mag = grid_data[n][2][0]
##          velocity_dir = grid_data[n][2][1]
##          angular_velocity = grid_data[n][4]
##
##          v_x = velocity_mag * math.cos(velocity_dir)
##          v_y = velocity_mag * math.sin(velocity_dir)
##          
##          E = E + (1/2 * mass * (v_x ** 2 + v_y ** 2) + 1/2 * interia * angular_velocity ** 2)
##          n = n + 1
##     T = 2/(3 * boltz_con) * E 
##     return T

"""
     Function localised pressure
-------------------------------------------------------------------------------------------------------
     - calculates the pressure of the simulation
          Inputs:
               - temperature
               - Volume
               - number of particles
               - change in position
               - force
          Outputs:
               - Pressure
"""

def localised_pressure(temperature, area, work, no_particles):
        k_B = 1.380649e-23
        P = (no_particles * np.float_(temperature) * k_B)/(np.float_(area)
                                                            ) + (1/(3*np.float_(area)) * np.float_(work))
        #       Not remove abs from pressure
        return P #   N/m



