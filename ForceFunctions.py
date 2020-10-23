"""
     ForceFunctions.py
"""
import math
import random
import numpy as np


from MeasureSimulation import *
"""
     LJ_Force
--------------------------------------------------------------------------------------------------
     - Put the Lennard-Jones potential to calculate the force of the
     interaction.
     - Inputs
          r - distance between two atoms
          epsilon - Potential energy at the equilibrium bond length
          r_0 - Distance at which the potential energy is zeor

     - Outputs
          f - force between two atoms, this is a van der Waals interaction and
          electron orbits of the atom
"""

def LJ_Force(r, epsilon, r_0):
     sigma = r_0/((2 ** (1/6)))
     f =  (12 * epsilon) * (((r_0**12)/(np.float_(r)**13)) - (r_0**6)/(np.float_(r)**7))
##     if abs(r) < 1.25 * sigma:
##          print('High Force = ' + str(f))
     return np.float_(f)


"""
     Torque
---------------------------------------------------------------------------------------------------
     - Calculate the torque for the equation.
     - Inputs:
          a - distance between the atoms in the molecule
          epsilon - potential energy at the equilibrium bond length
          r_0 - distance at which the potential energy is zero
          r_a - magnitude of atom a
          r_b - magnitude of atom b
          alpha_a - angle of atom a
          alpha_b - angle of atom b

     - Outputs:
          t - torque of the molecule
"""

def Torque (a, epsilon, r_0, r_a, r_b, alpha_a, alpha_b):
     t = (a / 2) * (LJ_Force(r_a, epsilon, r_0) * np.float_(math.sin(alpha_a)
                                                            ) + LJ_Force(r_b, epsilon, r_0) * np.float_(math.sin(alpha_b)))
     return t


              
"""
     Inertia
------------------------------------------------------------------------------------------------------------
     - Calculate the inertia for the equation.
     - Inputs:
          a - distance between two atoms in a molecule.
          m - reduced mass of the molecule.

     - Outputs:
          I - inertia of the diatomic molecule
"""

def Inertia (m, a):
     mu = (np.float_(m) * np.float_(m))/(np.float_(m) + np.float_(m))
     I = mu * a ** 2
     return I

"""
     sum_of_force
-------------------------------------------------------------------------------------------------------------
     Sums the forces between the atoms.  
     Inputs:
          - distance between atoms (list inputs)
          - epsilon
          - r_0
     Outputs:
          - Force = [x_force, y_force]
"""

def sum_of_force(list_of_distant, epsilon, r_0):
     #    Set up counters and set it to zero
     i = 0
     x_force = []
     y_force = []
     #    Calculate the work
     work = []
     mol_pos = []

     #    While look with the i counter
     while len(list_of_distant[1][0]) > i:   #while len(list_of_distant[1][0]) > i:
          atom_1_mag = np.float_(list_of_distant[1][0][i])
          atom_2_mag = np.float_(list_of_distant[2][0][i])
          atom_1_dir = np.float_(list_of_distant[1][1][i])
          atom_2_dir = np.float_(list_of_distant[2][1][i])
          
          x_force.append((LJ_Force(atom_1_mag, epsilon, r_0)* math.cos(atom_1_dir)) +
                         LJ_Force(atom_2_mag, epsilon, r_0) * math.cos(atom_2_dir))
          y_force.append((LJ_Force(atom_1_mag, epsilon, r_0) * math.sin(atom_1_dir))+
                          LJ_Force(atom_2_mag, epsilon, r_0) * math.sin(atom_2_dir))

          mol_pos.append([(np.float_(atom_1_mag) * math.cos(atom_1_dir) +
                                     (np.float_(atom_2_mag) * math.cos(atom_2_dir))),
                          (np.float_(atom_1_mag) * math.sin(atom_1_dir) +
                                     (np.float_(atom_2_mag) * math.sin(atom_2_dir)))])
          #print(mol_pos[-1])
          
          work.append(abs(np.dot(mol_pos[i], [x_force[i], y_force[i]])))
          #print(work[-1])
          
          i = i + 1

##     print(atom_1_mag)
##     print(atom_1_dir)

     work_total = sum(work)
     x_force_int = sum(x_force)
     y_force_int = sum(y_force)
     #    Return with force
     return [[x_force_int, y_force_int], work_total]


def sum_potential_energy(list_of_distant, epsilon, r_0):
     #    Set up counters and set it to zero
     i = 0
     x_energy = []
     y_energy = []
     mag_energy = []
     #    While look with the i counter
     while len(list_of_distant[1][0]) > i:
          x_energy.append((potential_energy(np.float_(list_of_distant[1][0][i]),
                                                epsilon, r_0) + potential_energy(np.float_(list_of_distant[2][0][i]),
                                                                         epsilon, r_0)))
          y_energy.append((potential_energy(list_of_distant[1][0][i],
                                                epsilon, r_0) + potential_energy(list_of_distant[2][0][i],
                                                                         epsilon, r_0)))

          mag_energy.append(math.sqrt(x_energy ** 2 + y_energy ** 2))
          
          i = i + 1

     potential_energy_sum = sum(mag_energy)
     #    Return with force
     return potential_energy_sum


"""
     sum_of_torque
------------------------------------------------------------------------------------------------------------
     Sums the torques acting on a molecule
     Inputs:
          - distance between atoms (list inputs)
          - epsilon
          - r_0
     Outputs:
          - torque_int
"""

def sum_of_torque(list_of_distant, epsilon, r_0, a):
     #    Counter_max and counter
     max_counter = len(list_of_distant[1][0])
##     print(len(list_of_distant[1][0]))
     i = 1
     torque_int = 0
##     polar_dist_a = []
##     polar_dist_b = [] 
     while max_counter > i:
          atom_1_mag = np.float_(list_of_distant[1][0][i])
          atom_2_mag = np.float_(list_of_distant[2][0][i])
          atom_1_dir = np.float_(list_of_distant[1][1][i])
          atom_2_dir = np.float_(list_of_distant[2][1][i])
##          print(str(atom_1_mag))
          
##          polar_dist_a.append([atom_1_mag, atom_1_dir])
##          polar_dist_b.append([atom_2_mag, atom_2_dir])
##          print(polar_dist_a)
          torque_int = torque_int + Torque(a, epsilon, r_0,
                                           atom_1_mag, atom_2_mag, atom_1_dir, atom_2_dir)
##          print(i)
          i = i + 1
     #print('Torque = ' + str(torque_int))
     return torque_int

"""
     check_force
--------------------------------------------------------------------------------------------------------------
     Checks if the force is above a threshold if so then then dt is shortened
     Inputs:
          - dt
          - list of force
     Outputs:
          - dt_new
"""

def check_force(dt, initial_dt, list_of_force, check_value, maximum_force):
     force_max = max(max(list_of_force))


     #print("The Max Force is = " + str(force_max))
     if check_value <= 10:
          dt = abs((initial_dt - dt) * (check_value/10))
          check_value = check_value + 1
          
          print("dt is now " + str(dt))
          print("The Check is " + str(check_value))
     
     elif abs(force_max) > (maximum_force):
          ratio_force = (maximum_force)/abs(force_max)
          dt = abs(ratio_force * dt)
          check_value = 1
          print("dt is lowered to " + str(abs(dt)))
          
     else:
          dt = initial_dt

     return [dt, check_value]


