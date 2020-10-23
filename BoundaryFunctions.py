import math
import numpy as np


"""
BoundaryFunctions.py

--------------------------------------------------------------------------------------------------------------

Three different functions within the doc
     - periodic_boundary_conditions
     - find_specific_grid
     - periodic_boundary_outside

"""




"""
     This file has the periodic boundary condtions, instead it being one massive file
     This function is broken up into three different parts
     ---------------------------------------------------------------------------------------------------------
     1.    Convert the input format into lists and specify the axis.
     2.    Check the data in the list to see if the positions are outside the conditions.
     3.    Convert the checked data back into the input data given.
     
"""

def periodic_boundary_conditions(specific_grid_data, length_of_axis, specific_axis, length_of_molecule):
     #    Initial Counter & Start list

     molecular_position = []
     atom_1_position = []
     atom_2_position = []

     #    Specific Axis String Shortened
     short_specific_axis = specific_axis.replace(" ", "")

     #    Specify axis
     if short_specific_axis == 'x' or short_specific_axis == 'xaxis' or short_specific_axis == 'x-axis':
          specific_axis_int = 0
##          print("x axis")
     elif short_specific_axis == 'y' or short_specific_axis == 'yaxis' or short_specific_axis == 'y-axis':
          specific_axis_int = 1
##          print("y axis")
     elif short_specific_axis == 'z' or short_specific_axis == 'zaxis' or short_specific_axis == 'z-axis':
          specific_axis_int = 2
     else:
          return None

     for i in specific_grid_data:
          atom_1_position = np.float_(i[1][0][specific_axis_int])
          atom_2_position = np.float_(i[1][1][specific_axis_int])
          molecule_position = np.float_(i[0][specific_axis_int])
          ratio_position = np.float_(atom_1_position / (length_of_axis))
          #print(np.float_(atom_1_position))

          if -1 * length_of_molecule - 1 * length_of_axis < atom_1_position < -1 * length_of_molecule:
               i[1][0][specific_axis_int] = atom_1_position + np.float_(length_of_axis)  #+ length_of_molecule#   atom_1_position + np.float_((length_of_axis + length_of_molecule))
               i[1][1][specific_axis_int] = atom_2_position + np.float_(length_of_axis) #+ length_of_molecule
               i[0][specific_axis_int] = molecule_position + np.float_(length_of_axis) #+ length_of_molecule
               print("Atom 1 Position Lesser Than " + str(atom_1_position) + ".  Atom 2 Position " + str(atom_2_position))
               print("Moved to: " + str(i[1][0][specific_axis_int]) + ".  " + str(i[1][1][specific_axis_int]))

          elif length_of_molecule + length_of_axis < atom_1_position < length_of_molecule + 2 * length_of_axis:
               i[1][0][specific_axis_int] = atom_1_position - np.float_(length_of_axis)  #- length_of_molecule   # atom_1_position - np.float_((length_of_axis + length_of_molecule))
               i[1][1][specific_axis_int] = atom_2_position - np.float_(length_of_axis)#- length_of_molecule
               i[0][specific_axis_int] = molecule_position - np.float_(length_of_axis)# - length_of_molecule
               print("Atom 1 Position  Greater Than " + str(atom_1_position) + ".  Atom 2 Position " + str(atom_2_position))
               print("Moved to: " + str(i[1][0][specific_axis_int]) + ".  " + str(i[1][1][specific_axis_int]))

          elif ratio_position <= -1 and -1 * length_of_molecule - 1 * length_of_axis > atom_1_position:
               i[1][0][specific_axis_int] = atom_1_position + np.float_(length_of_axis) * abs(math.floor(ratio_position)) #+ length_of_molecule
               i[1][1][specific_axis_int] = atom_2_position + np.float_(length_of_axis) * abs(math.floor(ratio_position)) #+ length_of_molecule
               i[0][specific_axis_int] = molecule_position + np.float_(length_of_axis) * abs(math.floor(ratio_position)) #+ length_of_molecule
               print("Atom 1 Position Too Small" + str(atom_1_position) + ".  Atom 2 Position " + str(atom_2_position))
               print("Moved to: " + str(i[1][0][specific_axis_int]) + ".  " + str(i[1][1][specific_axis_int]))

          elif ratio_position >= 2 and atom_1_position > length_of_molecule + 2 * length_of_axis:
               i[1][0][specific_axis_int] = atom_1_position + np.float_(length_of_axis)  * -1 * abs(math.floor(ratio_position))# - length_of_molecule
               i[1][1][specific_axis_int] = atom_2_position + np.float_(length_of_axis)  * -1 * abs(math.floor(ratio_position)) #- length_of_molecule
               i[0][specific_axis_int] = molecule_position + np.float_(length_of_axis) * -1 * abs(math.floor(ratio_position)) #- length_of_molecule
               print("Atom 1 Position Too Large " + str(atom_1_position) + ".  Atom 2 Position " + str(atom_2_position))
               print("Moved to: " + str(i[1][0][specific_axis_int]) + ".  " + str(i[1][1][specific_axis_int]))

          else:
               i[1][0][specific_axis_int] = atom_1_position 
               i[1][1][specific_axis_int] = atom_2_position 
               i[0][specific_axis_int] = molecule_position
               #print("Atom 1 Position " + str(atom_1_position) + ".  Atom 2 Position " + str(atom_2_position))
               #print("Good")


     return specific_grid_data
               


##
##     while i != len(specific_grid_data):
##          atom_1_position = specific_grid_data[i][1][0][specific_axis_int]
##          atom_2_position = specific_grid_data[i][1][1][specific_axis_int]
##          molecule_position = specific_grid_data[i][0][specific_axis_int]
##
##          if 0 > atom_1_position and 0 > atom_2_position: #and 0 > molecule_position:
##               if -1 <= atom_1_position / length_of_axis and -1 <= atom_2_position / length_of_axis:
##                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position + length_of_axis
##                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position + length_of_axis
##                    
##               else:
##                    atom_1_multiply = abs(math.ceil(atom_1_position / length_of_axis))
##                    atom_2_multiply = abs(math.ceil(atom_2_position / length_of_axis))
##                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position + length_of_axis * atom_1_multiply
##                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position + length_of_axis * atom_2_multiply
##
##               specific_grid_data[i][0][specific_axis_int] = (specific_grid_data[i][1][0][specific_axis_int] +
##                                                              specific_grid_data[i][1][1][specific_axis_int]) / 2
##               
##               print("Lesser than: Initial Position " + str(molecule_position) + ".  Moved to Position " + str(specific_grid_data[i][0][specific_axis_int]))
##               print("Lesser than: Initial Position Atom 1 " + str(atom_1_position) + ".  Moved to Position " + str(specific_grid_data[i][1][0][specific_axis_int]))
##               print("Lesser than: Initial Position Atom 2 " + str(atom_2_position) + ".  Moved to Position " + str(specific_grid_data[i][1][1][specific_axis_int]))
##
##
##          elif length_of_axis < atom_1_position and length_of_axis < atom_2_position: #and length_of_axis < molecule_position:
##               if 2 > atom_1_position / length_of_axis and 2 > atom_2_position / length_of_axis:
##                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position - length_of_axis
##                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position - length_of_axis
##                    
##               else:
##                    atom_1_multiply = abs(math.ceil(atom_1_position / length_of_axis))
##                    atom_2_multiply = abs(math.ceil(atom_2_position / length_of_axis))
##                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position - length_of_axis * atom_1_multiply
##                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position - length_of_axis * atom_2_multiply
##                                   #print("Greater than: Initial Position " + str(molecule_position) + ".  Moved to Position " + str(specific_grid_data[i][0][specific_axis_int]))
##
##               specific_grid_data[i][0][specific_axis_int] = (specific_grid_data[i][1][0][specific_axis_int] +
##                                                         specific_grid_data[i][1][1][specific_axis_int]) / 2
##
##               #print("Too Big")
##
##          elif 0 <= atom_1_position <= length_of_axis or 0 <= atom_2_position <= length_of_axis:
##               specific_grid_data[i][1][0][specific_axis_int] = atom_1_position 
##               specific_grid_data[i][1][1][specific_axis_int] = atom_2_position
##               specific_grid_data[i][0][specific_axis_int] = molecule_position
##               #print("Stays the same. " + str(molecule_position))
##
##          
##          
##          i = i + 1
##
##     #    Return with the specific grid
##     return specific_grid_data


"""
     Find Specific Grid
     ------------------------------------------------------------------------------------------------------------------
     Function finds the specific grid within the grid_data

     Inputs:        grid_data, specific_string
     Outputs:     Final_specific_grid

"""

def find_specific_grid(grid_data, grid_specific):
     #    Max Counter
     counter = len(grid_data)
     i = 0
     final_grid = []
     
     while i != counter:
          if grid_data[i][6] == grid_specific:
               final_grid.append(grid_data[i])
          i = i + 1
     return final_grid



"""
     Outside Boundary Periodic Conditions
     ------------------------------------------------------------------------------------------------------------------
     Copy the positions, velocity, angular displacement, angular velocity
     1.    Copy the Unit Cell
     2.    Return grid data
"""

def periodic_boundary_outside(cell_data, length_of_x_axis, length_of_y_axis, length_of_molecule):
        i = 0

        x_axis = length_of_x_axis# + length_of_molecule
        y_axis = length_of_y_axis# + length_of_molecule
        position_lists = []
        atom_lists = []
        phi = []
        velocity = []
        omega = []
        n = []
        while i < len(cell_data):
                position_lists.append(cell_data[i][0])
                atom_lists.append(cell_data[i][1])
                velocity.append(cell_data[i][2])
                phi.append(cell_data[i][3])
                omega.append(cell_data[i][4])
                n.append(cell_data[i][5])
                i = i + 1
                
        middle_positions = position_lists
        i = 0
        grid_molecule_list = []
        grid_atom_list = []
        grid_phi_list = []
        grid_velocity_list = []
        grid_omega_list = []
        grid_n = []
        grid_specific = []
        grid_format = []
        
        while i != len(position_lists):
                #Top-left Positions (x-xaxis, y+yaxis)
                grid_molecule_list.append([middle_positions[i][0] - x_axis, middle_positions[i][1] + y_axis])
                grid_atom_list.append([[atom_lists[i][0][0] - x_axis, atom_lists[i][0][1] + y_axis],
                                       [atom_lists[i][1][0] - x_axis, atom_lists[i][1][1] + y_axis]])

                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])                
                grid_specific.append("TL Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Top-mid Positions (x, y+yaxis)
                grid_molecule_list.append([middle_positions[i][0], middle_positions[i][1] + y_axis])
                grid_atom_list.append([[atom_lists[i][0][0], atom_lists[i][0][1] + y_axis],
                                       [atom_lists[i][1][0], atom_lists[i][1][1] + y_axis]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("TM Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Top-right Positions (x+xaxis, y+yaxis)
                grid_molecule_list.append([middle_positions[i][0] + x_axis, middle_positions[i][1] + y_axis])
                grid_atom_list.append([[atom_lists[i][0][0] + x_axis, atom_lists[i][0][1] + y_axis],
                                      [atom_lists[i][1][0] + x_axis, atom_lists[i][1][1] + y_axis]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("TR Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Mid-left Positions (x-xaxis, y)
                grid_molecule_list.append([middle_positions[i][0] - x_axis, middle_positions[i][1]])
                grid_atom_list.append([[atom_lists[i][0][0] - x_axis, atom_lists[i][0][1]],
                                       [atom_lists[i][1][0] - x_axis, atom_lists[i][1][1]]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("ML Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Unit Cell
        
                grid_molecule_list.append(middle_positions[i])
                grid_atom_list.append(atom_lists[i])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("Unit Cell")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])
                
                #Mid-right Positions (x+xaxis, y)
                grid_molecule_list.append([middle_positions[i][0] + x_axis, middle_positions[i][1]])
                grid_atom_list.append([[atom_lists[i][0][0] + x_axis, atom_lists[i][0][1]],
                                       [atom_lists[i][1][0] + x_axis, atom_lists[i][1][1]]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("MR Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Bottom-left Positions (x-xaxis, y-yaxis)
                grid_molecule_list.append([middle_positions[i][0] - x_axis, middle_positions[i][1] - y_axis])
                grid_atom_list.append([[atom_lists[i][0][0] - x_axis, atom_lists[i][0][1] - y_axis],
                                       [atom_lists[i][1][0] - x_axis, atom_lists[i][1][1] - y_axis]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("BL Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Bottom-mid Positions (x, y-yaxis)
                grid_molecule_list.append([middle_positions[i][0], middle_positions[i][1] - y_axis])
                grid_atom_list.append([[atom_lists[i][0][0], atom_lists[i][0][1] - y_axis],
                                       [atom_lists[i][1][0], atom_lists[i][1][1] - y_axis]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("BM Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])

                #Bottom-right Positions (x+xaxis, y-yaxis)
                grid_molecule_list.append([middle_positions[i][0] + x_axis, middle_positions[i][1] - y_axis])
                grid_atom_list.append([[atom_lists[i][0][0] + x_axis, atom_lists[i][0][1] - y_axis],
                                       [atom_lists[i][1][0] + x_axis, atom_lists[i][1][1] - y_axis]])
                grid_velocity_list.append([velocity[i][0], velocity[i][1]])
                grid_phi_list.append(phi[i])
                grid_omega_list.append(omega[i])
                grid_n.append(n[i])       
                grid_specific.append("BR Grid")

                grid_format.append([grid_molecule_list[-1], grid_atom_list[-1],
                                    grid_velocity_list[-1], grid_phi_list[-1],
                                    grid_omega_list[-1],
                                    grid_n[-1], grid_specific[-1]])
                
                i=i+1

        return grid_format
     
"""
     Triatomic Molecule
"""
     
def periodic_boundary_conditions_triatomic(specific_grid_data, length_of_axis, specific_axis):
     #    Initial Counter & Start List
     i = 0
     molecular_position = []
     atom_1_position = []
     atom_2_position = []

     #    Specific Axis String shortened
     short_specific_axis = specific_axis.replace(" ", "")

     #    Specify the axis
     if short_specific_axis == 'x' or short_specific_axis == 'xaxis' or short_specific_axis == 'x-axis':
          specific_axis_int = 0

     elif short_specific_axis == 'y' or short_specific_axis == 'yaxis' or short_specific_axis == 'y-axis':
          specific_axis_int = 1

     elif short_specific_axis == 'z' or short_specific_axis == 'zaxis' or short_specific_axis == 'z-axis':
          specific_axis_int = 2

     else:
          return None

     while i != len(specific_grid_data):
          atom_1_position = specific_grid_data[i][1][0][specific_axis_int]
          atom_2_position = specific_grid_data[i][1][1][specific_axis_int]
          atom_3_position = specific_grid_data[i][1][2][specific_axis_int]
          molecule_position = specific_grid_data[i][0][specific_axis_int]

          if 0 > atom_1_position and 0 > atom_2_position and 0 > atom_3_position:
               if -1 <= atom_1_position / length_of_axis and -1 <= atom_2_position / length_of_axis:
                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position + length_of_axis
                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position + length_of_axis
                    specific_grid_data[i][1][2][specific_axis_int] = atom_3_position + length_of_axis
                    specific_grid_data[i][0][specific_axis_int] = molecule_position + length_of_axis
               else:
                    atom_1_multiply = abs(math.floor(atom_1_position / length_of_axis))
                    atom_2_multiply = abs(math.floor(atom_2_position / length_of_axis))
                    atom_3_multiply = abs(math.floor(atom_3_position / length_of_axis))
                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position + length_of_axis * atom_1_multiply
                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position + length_of_axis * atom_2_multiply
                    specific_grid_data[i][1][2][specific_axis_int] = atom_3_position + length_of_axis * atom_3_multiply
                    specific_grid_data[i][0][specific_axis_int] = (specific_grid_data[i][1][0][specific_axis_int] + specific_grid_data[i][1][1][specific_axis_int] + specific_grid_data[i][1][2][specific_axis_int])/3
               print("Initial Position " + str(molecule_position) + ".  Moved to Position " + str(specific_grid_data[i][0][specific_axis_int]))

          elif length_of_axis < atom_1_position and length_of_axis < atom_2_position and length_of_axis < atom_3_position:
               if 2 > atom_1_position / length_of_axis and 2 > atom_2_position / length_of_axis:
                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position - length_of_axis
                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position - length_of_axis
                    specific_grid_data[i][1][2][specific_axis_int] = atom_3_position - length_of_axis
                    specific_grid_data[i][0][specific_axis_int] = molecule_position - length_of_axis
               else:
                    atom_1_multiply = math.ceil(atom_1_position / length_of_axis)
                    atom_2_multiply = math.ceil(atom_2_position / length_of_axis)
                    atom_3_multiply = math.ceil(atom_3_position / length_of_axis)
                    specific_grid_data[i][1][0][specific_axis_int] = atom_1_position - length_of_axis * atom_1_multiply
                    specific_grid_data[i][1][1][specific_axis_int] = atom_2_position - length_of_axis * atom_2_multiply
                    specific_grid_data[i][1][2][specific_axis_int] = atom_3_position - length_of_axis * atom_3_multiply
                    specific_grid_data[i][0][specific_axis_int] = (specific_grid_data[i][1][0][specific_axis_int] + specific_grid_data[i][1][1][specific_axis_int] + specific_grid_data[i][1][2][specific_axis_int])/3
               print("Initial Position " + str(molecule_position) + ".  Moved to Position " + str(specific_grid_data[i][0][specific_axis_int]))

               #print("Too Big")

          elif 0 <= atom_1_position <= length_of_axis or 0 <= atom_2_position <= length_of_axis and 0 <= atom_3_position <= length_of_axis:
               specific_grid_data[i][1][0][specific_axis_int] = atom_1_position 
               specific_grid_data[i][1][1][specific_axis_int] = atom_2_position
               specific_grid_data[i][1][2][specific_axis_int] = atom_3_position
               specific_grid_data[i][0][specific_axis_int] = molecule_position

          i = i + 1

     #    Return with the specific grid
     return specific_grid_data


##def periodic_boundary_conditions(specific_grid_data, length_of_axis, specific_axis):
##     #    Initial Counter & Start List
##     i = 0
##     molecular_position = []
##     atom_1_position = []
##     atom_2_position = []
##
##     #    Specific Axis String shortened
##     short_specific_axis = specific_axis.replace(" ", "")
##
##     #    Specify the axis
##     if short_specific_axis == 'x' or short_specific_axis == 'xaxis' or short_specific_axis == 'x-axis':
##          specific_axis_int = 0
##
##     elif short_specific_axis == 'y' or short_specific_axis == 'yaxis' or short_specific_axis == 'y-axis':
##          specific_axis_int = 1
##
##     elif short_specific_axis == 'z' or short_specific_axis == 'zaxis' or short_specific_axis == 'z-axis':
##          specific_axis_int = 2
##
##     else:
##          return None
