#    LoadFunction
import numpy as np
import csv

def load_csv(file_name):
     f = open(file_name, 'r')

     listed_csv = []

     with f:
          reader = csv.reader(f, delimiter = ',')

          for row in reader:
               listed_csv.append(row)

     #    Note remove the first row
##     listed_csv.remove([['dt', 'molecule_pos', 'atom_1_pos', 'atom_2_pos',
##                      'molecule_vel', 'molecule_phi', 'molecule_omega', 'temperature',
##                     'pressure']])
     return listed_csv

def format_loaded_data(csv_data):
     ###  Note that the csv data is formatted like
     #      ["['dt', 'molecule_pos', 'atom_1_pos',
     #      'atom_2_pos', 'molecule_vel', 'molecule_phi',
     #      'molecule_omega', 'temperature', 'pressure']"]
     #    csv_unformatted
     #csv_data.remove(["['dt', 'molecule_pos', 'atom_1_pos', 'atom_2_pos', 'molecule_vel', 'molecule_phi', 'molecule_omega', 'temperature', 'pressure']"])
     max_counter = len(csv_data)
     format_csv_data = []
     i = 1
     #n = 0
     dt = []
     mol_pos = []
     a_1_pos = []
     a_2_pos = []
     mol_vel = []
     mol_phi = []
     mol_ome = []
     temp  = []
     pressure = []
     total_energy = []
     while i < max_counter:
          n = 0
          dt.append(np.float_(csv_data[i][0]))
          temp.append(np.float_(csv_data[i][7]))
          pressure.append(np.float_(csv_data[i][8]))
          total_energy.append(np.float_(csv_data[i][9]))
          #    Temperary list
          m_pos_t = []
          a_1_pos_t = []
          a_2_pos_t = []
          m_vel_t = []
          m_phi_t = []
          m_ohm_t = []
          while n < len(csv_data[i][1]):
               m_pos_t.append(csv_data[i][1].split('], ['))
               a_1_pos_t.append(csv_data[i][2].split('], ['))
               a_2_pos_t.append(csv_data[i][3].split('], ['))
               m_vel_t.append(csv_data[i][4].split('], ['))
               #print(m_vel_t)
               m_phi_t.append(csv_data[i][5].split(', '))
               #print(m_phi_t)
               m_ohm_t.append(csv_data[i][6].split(', '))
               #print(m_ohm_t)
               nn = 0
               while nn < len(m_pos_t[n]):
                    #print(m_pos_t[n][nn])
                    m_pos_t[n][nn] = m_pos_t[n][nn].split(', ')
                    a_1_pos_t[n][nn] = a_1_pos_t[n][nn].split(', ')
                    a_2_pos_t[n][nn] = a_2_pos_t[n][nn].split(', ')
                    m_vel_t[n][nn] = m_vel_t[n][nn].split(', ')
                    #print(m_vel_t[n][nn])
                    #m_phi_t[n][nn] = m_phi_t[n][nn].split(', ')
                    #m_ohm_t[n][nn] = m_ohm_t[n][nn].split(', ')
                    #print(m_pos_t[n][nn])

                    nn += 1
               n += 1
##          mol_pos.append(m_pos_t)
##          a_1_pos.append(a_1_pos_t)
##          a_2_pos.append(a_2_pos_t)
##          mol_vel.append(m_vel_t)
##          mol_phi.append(m_phi_t)
##          mol_ome.append(m_ohm_t)

          format_csv_data.append([dt[-1],
                                  m_pos_t, a_1_pos_t,
                                  a_2_pos_t, m_vel_t, m_phi_t, m_ohm_t,
                                  temp[-1], pressure[-1],
                                  total_energy[-1]])
          print(i)
          i += 1
     #    Return the format csv data
     return format_csv_data


def string_data_to_list(formatted_data_csv):
     max_counter = len(formatted_data_csv)
     i = 0

     dt = []
     m_pos = []
     a_1_pos = []
     a_2_pos = []
     m_v = []
     m_phi = []
     m_ohm = []
     T = []
     P = []
     total_energy = []

     float_data_csv = []
     
     while i < max_counter:
          dt.append(formatted_data_csv[i][0])
          T.append(formatted_data_csv[i][7])
          P.append(formatted_data_csv[i][8])
          total_energy.append(formatted_data_csv[i][9])

          n = 0
          while n < len(formatted_data_csv[i][1][0]):
##               print(formatted_data_csv[i][1][0][n][0][0]) 
               #print(formatted_data_csv[i][1][0][n][0])
               if formatted_data_csv[i][1][0][n][0][0] == '[':
                    #print(formatted_data_csv[i][1][0][n])
                    formatted_data_csv[i][1][0][n][0] = formatted_data_csv[i][1][0][n][0].replace('[[', '')

##               print(formatted_data_csv[i][1][0][n][-1][-1])
               ##print(formatted_data_csv[i][1][0][n][-1][-1])
               if formatted_data_csv[i][1][0][n][-1][-1] == ']':
                    formatted_data_csv[i][1][0][n][-1] = formatted_data_csv[i][1][0][n][-1].replace(']]', '')

               if formatted_data_csv[i][2][0][n][0][0] == '[':
                    #print(formatted_data_csv[i][2][0][n])
                    formatted_data_csv[i][2][0][n][0] = formatted_data_csv[i][2][0][n][0].replace('[[', '')

##               print(formatted_data_csv[i][1][0][n][-1][-1])
               if formatted_data_csv[i][2][0][n][-1][-1] == ']':
                    formatted_data_csv[i][2][0][n][-1] = formatted_data_csv[i][2][0][n][-1].replace(']]', '')


               if formatted_data_csv[i][3][0][n][0][0] == '[':
                    #print(formatted_data_csv[i][3][0][n])
                    formatted_data_csv[i][3][0][n][0] = formatted_data_csv[i][3][0][n][0].replace('[[', '')

##               print(formatted_data_csv[i][1][0][n][-1][-1])
               if formatted_data_csv[i][3][0][n][-1][-1] == ']':
                    formatted_data_csv[i][3][0][n][-1] = formatted_data_csv[i][3][0][n][-1].replace(']]', '')

               if formatted_data_csv[i][4][0][n][0][0] == '[':
                    #print(formatted_data_csv[i][4][0][n])
                    formatted_data_csv[i][4][0][n][0] = formatted_data_csv[i][4][0][n][0].replace('[[', '')

               if formatted_data_csv[i][4][0][n][-1][-1] == ']':
                    formatted_data_csv[i][4][0][n][-1] = formatted_data_csv[i][4][0][n][-1].replace(']]', '')

####               print(formatted_data_csv[i][2][0][n][0][0])
##               if formatted_data_csv[i][2][0][n][0][0] == '[':
##                    formatted_data_csv[i][2][0][n][0] = formatted_data_csv[i][2][0][n][0].replace('[[', '')
##                    
##               if formatted_data_csv[i][2][0][n][-1][-1] == ']':
##                    formatted_data_csv[i][2][0][n][-1] = formatted_data_csv[i][2][0][n][-1].replace(']]', '')
##
##               if formatted_data_csv[i][3][0][n][0][0] == '[':
##                    formatted_data_csv[i][3][n][0][0] = formatted_data_csv[i][3][n][0][0].replace('[[', '')
##                    
##               if formatted_data_csv[i][3][0][n][-1][-1] == ']':
##                    formatted_data_csv[i][3][n][0][-1] = formatted_data_csv[i][3][n][0][-1].replace(']]', '')
##
##               if formatted_data_csv[i][4][0][n][0][0] == '[':
##                    print(formatted_data_csv[i][4][0][n])
##                    formatted_data_csv[i][4][n][0][0] = formatted_data_csv[i][4][n][0][-1].replace('[', '')
##                    
##               if formatted_data_csv[i][4][n][0][-1][-1] == ']':
##                    formatted_data_csv[i][4][n][0][-1] = formatted_data_csv[i][4][n][0][-1].replace(']', '')
                    
               #print(formatted_data_csv[i][5][0][n][0])
               if formatted_data_csv[i][5][0][n][0] == '[':
                    formatted_data_csv[i][5][0][n] = formatted_data_csv[i][5][0][n].replace('[', '')
                    #print(formatted_data_csv[i][5][0][n])
                    
               if formatted_data_csv[i][5][0][n][-1] == ']':
                    formatted_data_csv[i][5][0][n] = formatted_data_csv[i][5][0][n].replace(']', '')

               if formatted_data_csv[i][6][0][n][0] == '[':
                    formatted_data_csv[i][6][0][n] = formatted_data_csv[i][6][0][n].replace('[', '')
                    
               if formatted_data_csv[i][6][0][n][-1] == ']':
                    formatted_data_csv[i][6][0][n] = formatted_data_csv[i][6][0][n].replace(']', '')

                    
               formatted_data_csv[i][1][0][n][0] = np.float_(formatted_data_csv[i][1][0][n][0])
##               print(formatted_data_csv[i][1][0][n])
##               print(formatted_data_csv[i][1][0][n][-1])
               formatted_data_csv[i][1][0][n][-1] = np.float_(formatted_data_csv[i][1][0][n][-1])
               
               formatted_data_csv[i][2][0][n][0] = np.float_(formatted_data_csv[i][2][0][n][0])
               formatted_data_csv[i][2][0][n][-1] = np.float_(formatted_data_csv[i][2][0][n][-1])
               
               formatted_data_csv[i][3][0][n][0] = np.float_(formatted_data_csv[i][3][0][n][0])
               formatted_data_csv[i][3][0][n][-1] = np.float_(formatted_data_csv[i][3][0][n][-1])

               formatted_data_csv[i][4][0][n][0] = np.float_(formatted_data_csv[i][4][0][n][0])
               formatted_data_csv[i][4][0][n][-1] = np.float_(formatted_data_csv[i][4][0][n][-1])

               #formatted_data_csv[i][4][n][0] = np.float_(formatted_data_csv[i][4][n][0])
               formatted_data_csv[i][5][0][n] = np.float_(formatted_data_csv[i][5][0][n])
               formatted_data_csv[i][6][0][n] = np.float_(formatted_data_csv[i][6][0][n])

               n += 1

          m_pos.append(formatted_data_csv[i][1][0])
          a_1_pos.append(formatted_data_csv[i][2][0])
          a_2_pos.append(formatted_data_csv[i][3][0])

          m_v.append(formatted_data_csv[i][4][0])
          m_phi.append(formatted_data_csv[i][5][0])
          m_ohm.append(formatted_data_csv[i][6][0])

          i += 1


     float_data_csv.append([dt, m_pos, a_1_pos, a_2_pos,
                            m_v, m_phi, m_ohm, T, P, total_energy])
     return float_data_csv



#    Let's Combine all of the Loading Functions
def load_csv_formatted(csv_filename):
     #    First Load CSV as a string
     string_csv = load_csv(csv_filename)
     #    Seperate the string into 9 different lists
     string_csv_list = format_loaded_data(string_csv)
     #    Conver the sting variables into floats
     loaded_output = string_data_to_list(string_csv_list)
     return loaded_output[0]
     
     
#    x = load_csv_formatted('100 ParticleDiatomic2D.csv')
