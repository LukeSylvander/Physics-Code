###             Impoer Functions
import random
import math
import numpy as np

from BoltzmannMaxwell import *
"""
PlaceMolecules

--------------------------------------------------------------------------------------------------------------

Four Different Functions
     - Find_Atoms
     - Place_Positions_Random
     - Positions_Grid
     - Place_Positions_Grid
     
"""


"""
Function Find_Atoms
-----------------------------------------------------------------------------------------------------------------
        - The position is at the centre of mass of the molecule

        - find_atoms(molecule_position, phi, a)
                Gives -    [[Atom 1 Position], [Atom 2 Position]]
"""

def find_atoms(molecule_position, molecule_phi, a):
        #       Define the constants used in the calculation
        half_length = a / 2
        relative_degree_1 = np.float_(molecule_phi - math.pi / 2)
        relative_degree_2 = np.float_(molecule_phi + math.pi / 2)

        #       Caclulate the Relative Position of the atoms
        relative_atom_1 = [np.float_(half_length * math.cos(relative_degree_1)),
                           np.float_(half_length * math.sin(relative_degree_1))]
        relative_atom_2 = [np.float_(half_length * math.cos(relative_degree_2)),
                           np.float_(half_length * math.sin(relative_degree_2))]

        #       Caclulate the atoms position
        atom_1 = [np.float_(molecule_position[0] + relative_atom_1[0]),
                  np.float_(molecule_position[1] + relative_atom_1[1])]
        atom_2 = [np.float_(molecule_position[0] + relative_atom_2[0]),
                  np.float_(molecule_position[1] + relative_atom_2[1])]

        #       Return Functon
        return [atom_1, atom_2]


"""
Function Place_Positions_Random
-------------------------------------------------------------------------------------------------------------
          - Functions is not used
"""

def place_positions_random(total_no, length_of_molecule):
        molecule_list = give_random_molecule(total_no)
        atom_list = give_atom_position(molecule_list, length_of_molecule)
        i = 0
        mol_str_list = []
        phi_str_list = []
        Atom_a_list = []
        Atom_b_list = []
        position_list = []
        while i != total_no:
                position_list.append([[mol_str_list[i], molecule_list[i][0]],
                                        [phi_str_list[i], molecule_list[i][1]],
                                        [Atom_a_list[i], atom_list[i][0]],
                                        [Atom_b_list[i],  atom_list[i][1]]])
                i = i + 1
        return position_list


"""
     Positions_Grid
----------------------------------------------------------------------------------------------------------
     - Inputs the total number of molecules, the length of the simulation
     - Outputs the positions of the molecules
"""

def positions_grid(total_no, length_of_x_axis, length_of_y_axis):
        length_of_side = length_of_x_axis
        total_no_square = total_no ** 0.5
        #       Note this measures the area/number of particles
        axis_grid_x = np.float_(length_of_x_axis / total_no_square)
        axis_grid_y = np.float_(length_of_y_axis / total_no_square)
        #       Counter & List position
        i = 0
        x = 0
        y = 0
        position_list = []
        y_axis = axis_grid_y
        while i != (total_no):
                x_axis = np.float_(x * axis_grid_x)
                y_axis = np.float_(y * axis_grid_y)
                x = x + 1
                if total_no_square <= x:
                        x = 0
                        #x_axis = x * axis_grid_x
                        y = y + 1
                
                position_list.append([x_axis, y_axis])
                i = i +1

        return position_list

"""
Place_Positions_Grid
-------------------------------------------------------------------------------------------------------------
          Inputs:
                    - Total Number of Molecules
                    - Length of x axis
                    - Length of y axis
                    - Starting Velocity
                    - Length of Molecules
          Outputs:
                    - [[Molecule Position], [Atoms Position], [velocity], [Phi], [Omega]]

"""

def place_positions_grid(total_no, length_of_x_axis, length_of_y_axis, a, m, T):
     position_grid = positions_grid(total_no, length_of_x_axis, length_of_y_axis)
     i = 0

     position_lists = []
     atom_list = []
     random_phi = []
     n = 0
     omega_list = []
     n_list = []
     velocity_list = []
     while i < total_no:
             x_mol_pos = position_grid[i][0]
             y_mol_pos = position_grid[i][1]
             #random_phi.append((random.randint(0, 2 * 10 ** 3) * 10 ** - 3) * math.pi)
             random_phi.append(0)
             atom_list.append(find_atoms(position_grid[i], random_phi[i], a))
             #print(atom_list[i])
             i = i + 1


     while n < total_no:
             omega_list.append(normalised_start_omega(T, m, a, total_no))
             n_list.append(n+1)
             velocity_list.append(normalised_start_velocity(T, m, total_no))
             n = n + 1

     i = 0
     final_position = []
     while i < total_no:
             final_position.append([position_grid[i], atom_list[i], velocity_list[i],
                                    random_phi[i], omega_list[i], n_list[i], "Unit Cell"])
             i = i +1
     
     return final_position



"""
change_velocity_keep_positions
-------------------------------------------------------------------------------------------------------------
          Inputs:
                    - Unit_Cell_Unit 
                    - Temperature
          Outputs:
                    - [position_grid[i], atom_list[i], velocity_list[i],
                            random_phi[i], omega_list[i], n_list[i], "Unit Cell"]
"""

def change_velocity_keep_positions(unit_cell, T, m, a):
     i = 0

     position_lists = []
     atom_list = []
     phi_list = []
     omega_list = []
     n_list = []
     velocity_list = []
     while i < len(unit_cell):
             position_lists.append(unit_cell[i][0])
             atom_list.append(unit_cell[i][1])

             phi_list.append(unit_cell[i][3])

             n_list.append(unit_cell[i][5])
             i = i + 1
             
     n = 0
     total_no = len(unit_cell)
     
     while n < total_no:
             omega_list.append(normalised_start_omega(T, m, a, total_no))
             
             velocity_list.append(normalised_start_velocity(T, m, total_no))
             n = n + 1

     i = 0
     final_position = []
     while i < total_no:
             final_position.append([position_lists[i], atom_list[i], velocity_list[i],
                                    phi_list[i], omega_list[i], n_list[i], "Unit Cell"])
             i = i +1
     
     return final_position










##"""
##Function Find_Triatomic_Atoms
##-----------------------------------------------------------------------------------------------------------------
##        - The position is at the centre of mass of the molecule
##
##        - find_atoms(molecule_position, phi, a)
##                Gives -    [[Atom 1 Position], [Atom 2 Position], [Atom 3 Position]]
##"""
##def find_tratomic_atoms(position, phi, a):
##        #       Define the constants
##        tri_r = math.acos(math.atan(((3) ** (1/2))/2)) * a/2
##        relative_angle_1 = phi + math.pi / 2
##        relative_angle_2 = phi + (math.pi * 7) / 6
##        relative_angle_3 = phi + (math.pi * 11) / 6
##
##        #       Calculate the relative position
##        relative_atom_1 = [tri_r * math.cos(relative_angle_1),
##                           tri_r * math.sin(relative_angle_1)]
##        relative_atom_2 = [tri_r * math.cos(relative_angle_2),
##                           tri_r * math.sin(relative_angle_2)]
##        relative_atom_3 = [tri_r * math.cos(relative_angle_3),
##                           tri_r * math.sin(relative_angle_3)]
##
##        #       Calculate the atoms Position
##        atom_1 = [position[0] + relative_atom_1[0],
##                  position[1] + relative_atom_1[1]]
##        atom_2 = [position[0] + relative_atom_2[0],
##                  position[1] + relative_atom_2[1]]
##        atom_3 = [position[0] + relative_atom_3[0],
##                  position[1] + relative_atom_3[1]]
##
##        return [atom_1, atom_2, atom_3]
##
##
##"""
##Place_Positions_Grid_Triatomic
##-------------------------------------------------------------------------------------------------------------
##          Inputs:
##                    - Total Number of Molecules
##                    - Length of x axis
##                    - Length of y axis
##                    - Starting Velocity
##                    - Length of Molecules
##          Outputs:
##                    - [[Molecule Position], [Atoms Position], [velocity], [Phi], [Omega]]
##
##"""
##
##def place_positions_grid_triatomic(total_no, length_of_x_axis, length_of_y_axis, starting_velocity, a):
##        position_grid = positions_grid(total_no, length_of_x_axis, length_of_y_axis)
##        velocity_direction = []
##
##        i = 0
##        position_lists = []
##        atom_list = []
##        phi = []
##
##        omega_list = 0
##        n_list = []
##        velocity_list = []
##
##        while i < total_no:
##                x_mol_pos = position_grid[i][0]
##                y_mol_pos = position_grid[i][1]
##                velocity_direction.append((random.randint(0, 2 * 10 ** 3) * 10 ** - 3) * math.pi)
##                #phi.append((random.randint(0, 2 * 10 ** 3) * 10 ** - 3) * math.pi)
##                phi.append(0)
##                atom_list.append(find_atoms(position_grid[i], phi[i], a))
##                omega_list.append(0)
##                n_list.append(i)
##                velocity_list.append([starting_velocity, velocity_direction[i]])
##                i = i + 1
##
##        i = 0
##        final_position = []
##
##        while i < total_no:
##                final_position.append([position_grid[i], atom_list[i], velocity_list[i],
##                                       phi[i], omega_list[i], n_list[i], "Unit Cell"])
##                i = i + 1
##        return final_position

