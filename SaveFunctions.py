#    SaveFunctions.py
"""
     Write functions to a csv function

     - append_new_line_csv
     - open_new_csv
     - save_csv
"""

#    Import csv functions
import csv
import numpy as np


def open_new_csv(file_name):
     row_diatomic = [['dt', 'molecule_pos', 'atom_1_pos', 'atom_2_pos',
                      'molecule_vel', 'molecule_phi', 'molecule_omega', 'temperature',
                     'pressure']]#, 'Total Energy']]

     with open(file_name, 'w+', newline = '') as write_head:
          write_head = csv.writer(write_head)
          write_head.writerow(row_diatomic)



def append_new_line_csv(file_name, appended_data):
     """
     Function append new row with data.
     -----------------------------------------------------------------
     """
     #    Open the file in append & read mode ('a+')
     with open(file_name, "a+", newline = '') as write_obj:
          #    Create a writer object from csv module
          csv_writer = csv.writer(write_obj)
          #    add contents of list as last row in the csv file
          csv_writer.writerow(appended_data)
     return


def save_data_csv(file_name, appended_data, P, temp, dt, counter):#, total_energy):
     if counter == True:
          open_new_csv(file_name)
          counter = False

     else:
          i = 0
          mol_pos = []
          a_1_pos = []
          a_2_pos = []
          mol_velo = []
          mol_phi = []
          mol_ohm = []
          while len(appended_data) > i:
               mol_pos.append(appended_data[i][0])
               a_1_pos.append(appended_data[i][1][0])
               a_2_pos.append(appended_data[i][1][1])
               mol_velo.append([np.float_(appended_data[i][2][0]), appended_data[i][2][1]])
               mol_phi.append(np.float_(appended_data[i][3]))
               mol_ohm.append(np.float_(appended_data[i][4]))

               i = i + 1
          data_structure = [dt, mol_pos, a_1_pos, a_2_pos,
                            mol_velo, mol_phi, mol_ohm, temp, P]#, total_energy]

          #    Save data as csv
          append_new_line_csv(file_name, data_structure)
          
     return counter


def read_line_csv(file_name):
     """
     Read the line csv
     -----------------------------------------------------------------
     """
     #    Open the file in append & read mode ('a+')
     with open(file_name, "r+", newline = '') as read_obj:
          appended_data = []
          #    Create a writer object from csv module
          csv_reader = csv.reader(read_obj)
          print(str(csv_reader))
          #    add contents of list as last row in the csv file
          appended_data.append(csv_reader)
          #csv_writer.writerow(appended_data)
     return appended_data




#x = load_csv_formatted('484ParticleDiatomic2D_atempt_3.csv')
               
               
               
               
                    
                    
                    
          
     
               
##x = load_csv('64ParticleDiatomic2D.csv')  
##          
##          
##format_x = format_loaded_data(x)
##
##
##float_format_x = string_data_to_list(format_x)



          
          




     

##def read_data_csv(file_name):
##     """
##          reads the csv data
##          -----------------------------------------
##          Imputs
##               - file_name
##          Outputs
##               - Unit.data
##     """
##     with open(file_name) as csv_file:
##          csv_reader = csv.reader(csv_file)
##          line_count = 0
##          data_structure = []
##          for row in csv_reader:
##               if line_counter == 0:
##                    line_counter += 1
##               else:
##                    data_structure.append(row)
##                    line_counter += 1
##     return data_structure


##def read_data_csv(file_name):
##     """
##          reads the csv data
##          -----------------------------------------
##          Imputs
##               - file_name
##          Outputs
##               - Unit.data
##     """
##     with open(file_name, newline= '') as File:
##          reader = csv.reader(FIle)
##          data_structure = []
##          line_counter = 0
##          for row in reader:
##               if line_counter == 0:
##                    line_counter += 1
##               else:
##                    data_structure.append(row)
##                    line_counter += 1
##     return data_structure


          
          
               







                    
                    


