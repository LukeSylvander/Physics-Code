import math
import numpy as np
"""
     Distance Between Functions
------------------------------------------------------------------------------------------------------------

     - Find Specific Counter
     - Distance between all Neighbours
     - Distance between closer neighbours
     - Closest neighbours

"""

"""
Find Specific Counter
------------------------------------------------------------------------------------------------------------
     Function:   Finds the specific grid within the grid_data

     Inputs:       grid_data, specific_string
     Outputs:    List positions of specific grid
"""

def find_specific_counter(grid_data, grid_specific):
     counter = len(grid_data)
     string = grid_specific
     counter_list = []
     i = 0
     while i != counter:
          if grid_data[i][6] == string:
               counter_list.append(i)
          i = i + 1
     return counter_list


"""
Distance between All Neighbours
----------------------------------------------------------------------------------------------------------
     Function:     Finds the distance between the atoms within the molecule
                              Interaction between all neighbours

     Inputs:         grid_ data, specific_int
     Outputs:      [[molecule_r, molecule_phi], [Atom_1_r, Atom_1_phi], [Atom_2_r, Atom_2_phi]]

     Note -
          When editing the code this is not true the output is in cart coordinates
"""

def dis_between_all_unit(unit_data, counter):
     """
     distance between molecules in untis

     """
     i = 0
     mol_pos = []
     atom_1_pos = []
     atom_2_pos = []
     while i != len(unit_data):
          mol_pos.append([np.float_(unit_data[i][0][0]),
                          np.float_(unit_data[i][0][1])])
          atom_1_pos.append([np.float_(unit_data[i][1][0][0]),
                             np.float_(unit_data[i][1][0][1])])
          atom_2_pos.append([np.float_(unit_data[i][1][1][0]),
                             np.float_(unit_data[i][1][1][1])])
          i += 1
     #    Calculate the Distance between units
     ii = 0

     dist_b_mol =[]
     dist_b_a1 = []
     dist_b_a2 = []
     while ii != len(unit_data):
          if ii != counter:
               dist_b_mol.append([np.float_(mol_pos[counter][0]) - np.float_(mol_pos[ii][0]),
                                  np.float_(mol_pos[counter][1]) - np.float_(mol_pos[ii][1])])
               dist_b_a1.append([np.float_(atom_1_pos[counter][0]) - np.float_(atom_1_pos[ii][0]),
                                 np.float_(atom_1_pos[counter][1]) - np.float_(atom_1_pos[ii][1])])
               dist_b_a2.append([np.float_(atom_2_pos[counter][0]) - np.float_(atom_2_pos[ii][0]),
                                 np.float_(atom_2_pos[counter][1]) - np.float_(atom_2_pos[ii][1])])
          ii += 1

     iii = 0
     while iii != len(unit_data):
          if iii != counter:
               dist_b_mol.append([np.float_(mol_pos[counter][0]) - np.float_(mol_pos[iii][0]),
                                  np.float_(mol_pos[counter][1]) - np.float_(mol_pos[iii][1])])
               dist_b_a1.append([np.float_(atom_1_pos[counter][0]) -
                                 np.float_(atom_2_pos[iii][0]),
                                 np.float_(atom_1_pos[counter][1]) -
                                 np.float_(atom_2_pos[iii][1])])
               dist_b_a1.append([np.float_(atom_2_pos[counter][0]) -
                                 np.float_(atom_1_pos[iii][0]),
                                 np.float_(atom_2_pos[counter][1]) -
                                 np.float_(atom_1_pos[iii][1])])
          iii += 1

     #    Return Function

     return [[dist_b_mol], [dist_b_a1], [dist_b_a2]]
     
                                       



def distance_between_all_neighbours(grid_data, specific_int):
     i = 0
     mol_positions = []
     atom_1_positions = []
     atom_2_positions = []
     while i != len(grid_data):
          mol_positions.append([np.float_(grid_data[i][0][0]),np.float_(grid_data[i][0][1])])
          atom_1_positions.append([np.float_(grid_data[i][1][0][0]),np.float_(grid_data[i][1][0][1])])
          atom_2_positions.append([np.float_(grid_data[i][1][1][0]),np.float_(grid_data[i][1][1][1])])
          i = i + 1
     i = 0
     dist_between_mol = []
     dist_between_atom_1 = []
     dist_between_atom_2 = []
     while i != len(grid_data):       #       change in position
          if i != specific_int:
               dist_between_mol.append([np.float_(mol_positions[specific_int][0]) -
                                        np.float_(mol_positions[i][0]),
                                        np.float_(mol_positions[specific_int][1]) -
                                        np.float_(mol_positions[i][1])])
               dist_between_atom_1.append([np.float_(atom_1_positions[specific_int][0]) -
                                           np.float_(atom_1_positions[i][0]),
                                           np.float_(atom_1_positions[specific_int][1]) -
                                           np.float_(atom_1_positions[i][1])])
               dist_between_atom_2.append([np.float_(atom_2_positions[specific_int][0]) -
                                           np.float_(atom_2_positions[i][0]),
                                           np.float_(atom_2_positions[specific_int][1]) -
                                           np.float_(atom_2_positions[i][1])])
          i = i + 1
     i = 0
     while i != len(grid_data):
          if i != specific_int:
               dist_between_atom_1.append([np.float_(atom_1_positions[specific_int][0]) -
                                           np.float_(atom_2_positions[i][0]),
                                           np.float_(atom_1_positions[specific_int][1]) -
                                           np.float_(atom_2_positions[i][1])])
               dist_between_atom_2.append([np.float_(atom_2_positions[specific_int][0]) -
                                           np.float_(atom_1_positions[i][0]),
                                           np.float_(atom_2_positions[specific_int][1]) -
                                           np.float_(atom_1_positions[i][1])])
          i = i + 1
     mol_r = []
     mol_phi = []
     atom_1_r = []
     atom_1_phi = []
     atom_2_r = []
     atom_2_phi= []
     i = 0
     while i < (len(grid_data) - 1):
          mol_r.append((np.float_(dist_between_mol[i][0]) ** 2 +
                        np.float_(dist_between_mol[i][1]) ** 2) ** 0.5)
          mol_phi.append(math.atan2(np.float_(dist_between_mol[i][1]),
                                    np.float_(dist_between_mol[i][0])))

          atom_1_r.append((np.float_(dist_between_atom_1[i][0]) ** 2 +
                           np.float_(dist_between_atom_1[i][1]) ** 2) ** 0.5)
          atom_1_phi.append(math.atan2(np.float_(dist_between_atom_1[i][1]),
                                       np.float_(dist_between_atom_1[i][0])))

          atom_2_r.append((np.float_(dist_between_atom_2[i][0]) ** 2 +
                           np.float_(dist_between_atom_2[i][1]) ** 2) ** 0.5)
          atom_2_phi.append(math.atan2(np.float_(dist_between_atom_2[i][1]),
                                       np.float_(dist_between_atom_2[i][0])))

          i = i + 1

     return [[mol_r, mol_phi], [atom_1_r, atom_1_phi], [atom_2_r, atom_2_phi]]


"""
Distance between Nearest Neighbours
----------------------------------------------------------------------------------------------------------
     Function:     Finds the distance between the nearest atoms within the molecule
                              Interaction between the nearest neighbours

     Inputs:         grid_ data, specific_int, max_range
     Outputs:      [[molecule_r, molecule_phi], [Atom_1_r, Atom_1_phi], [Atom_2_r, Atom_2_phi]]
"""

##def distance_between_nearest_neighbours(grid_data, specific_int, max_range):
##     i = 0
##     mol_positions = []
##     atom_1_positions = []
##     atom_2_positions = []
##        
##     while i != len(grid_data):
##          mol_positions.append([np.float_(grid_data[i][0][0]),np.float_(grid_data[i][0][1])])
##          atom_1_positions.append([np.float_(grid_data[i][1][0][0]),np.float_(grid_data[i][1][0][1])])
##          atom_2_positions.append([np.float_(grid_data[i][1][1][0]),np.float_(grid_data[i][1][1][1])])
##          i = i + 1
##     
##     i = 0
##     dist_between_mol = []
##     dist_between_atom_1 = []
##     dist_between_atom_2 = []
##     while i < len(grid_data):
##          if i != specific_int:
##               if abs(mol_positions[specific_int][0] - mol_positions[i][0]) <= max_range:
##                    dist_between_mol.append([mol_positions[specific_int][0] - mol_positions[i][0],
##                                             mol_positions[specific_int][1] - mol_positions[i][1]])
##                    dist_between_atom_1.append([atom_1_position[specific_int][0] - atom_1_position[i][0],
##                                                atom_1_position[specific_int][1] - atom_1_position[i][1]])
##                    dist_between_atom_2.append([atom_2_position[specific_int][0] - atom_2_position[i][0],
##                                                atom_2_position[specific_int][1] - atom_2_position[i][1]])
##          i = i + 1
##     i = 0
##     while i < len(grid_data):
##          if i != specific_int:
##               if abs(mol_positions[specific_int][0] - mol_positions[i][0]) <= max_range:
##                    dist_between_atom_1.append([atom_1_positions[specific_int][0] - atom_2_positions[i][0],
##                                                atom_1_positions[specific_int][1] - atom_2_positions[i][1]])
##                    dist_between_atom_1.append([atom_2_positions[specific_int][0] - atom_1_positions[i][0],
##                                                atom_2_positions[specific_int][1] - atom_1_positions[i][1]])
##          i = i + 1
##     i = 0
##     mol_r = []
##     mol_alpha = []
##     atom_1_r = []
##     atom_1_alpha = []
##     atom_2_r = []
##     atom_2_alpha = []
##     while i <= len(grid_data):
##          mol_r.append(dist_between_mol[i][0])
##          mol_alpha.append(dist_between_mol[i][1])
##          atom_1_r.append(dist_between_atom_1[i][0])
##          atom_1_alpha.append(dist_between_atom_1[i][1])
##          atom_2_r.append(dist_between_atom_2[i][0])
##          atom_2_alpha.append(dist_between_atom_2[i][1])
##          i = i + 1
##     
##     sum_mol = [sum(mol_r), sum(mol_alpha)]
##
##     sum_atom_1 = [sum(atom_1_r), sum(atom_1_alpha)]
##
##     sum_atom_2 = [sum(atom_2_r), sum(atom_2_alpha)]
##
##     return [sum_mol, sum_atom_1, sum_atom_2]
##
##"""
##Closest Neighbours Int
##----------------------------------------------------------------------------------------------------------
##     Function:     Function returns with the number of near neighbours near the molecule
##
##     Inputs:         grid_ data, specific_int, max_range
##     Outputs:      Number of neighbours near
##"""

def near_neighbours_int(grid_data, specific_int, max_range):
     i = 0
     counter_neighbours = 0
     mol_positions = grid_data
     while i < len(grid_data):
          if 1 != specific_int:
               if abs(np.float_(mol_positions[specific_int][0][0]) -
                      np.float_(mol_positions[i][0][0])) <= max_range and abs(
                           np.float_(mol_positions[specific_int][0][1]) -
                           np.float_(mol_positions[i][0][1])) <= np.float_(max_range):
                    counter_neighbours = counter_neighbours + 1
          i = i + 1
     return counter_neighbours



"""

"""

def distance_between_unit_molecules(unit_data):
     i = 0
     mol_positions = []

     while i != len(unit_data):
          mol_positions.append([np.float_(unit_data[i][0][0]), np.float_(unit_data[i][0][1])])

          i = i + 1
     i = 0
     dist_between_mol = []

     while i != len(unit_data):       #       change in position
          ii = 0
          while ii != len(unit_data):
               
               if i != ii:
                    dist_between_mol.append([mol_positions[ii][0] - mol_positions[i][0],
                                             mol_positions[ii][1] - mol_positions[i][1]])
               ii = ii + 1
          i = i + 1

     return dist_between_mol
