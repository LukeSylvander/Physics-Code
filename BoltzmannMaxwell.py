#    Import
import random
import math
import numpy as np
import scipy.stats 

"""
     BoltzmannMaxwell.py
----------------------------------------------------------------------------------------


"""

####def MaxwellBoltzmannFunction(temp_aim, m, a):
####     k = 1.380649e-23
####     pi = math.pi
##
##def normalised_start_velocity(T, m):
##     k = 1.380649e-23 #  Boltzmann Constant
##     half_a = a/2
##     mu = 2 * math.sqrt((k * T)/m) * math.sqrt(2/(math.pi)) * math.sqrt((3 * k)/2)
##     ##mu = 2 * math.sqrt((k * T)/m) * math.sqrt(2/(math.pi)) 
##     #print('Mu = ' + str(mu))
##     sigma = math.sqrt((k * T)/m * (3 * math.pi - 8)/math.pi) * math.sqrt((3 * k)/2)
##     ##sigma = math.sqrt((k * T)/m * (3 * math.pi - 8)/math.pi)
##     #print('Sigma' + str(sigma))
##     pos_neg_dir = [-1, 1]
##     
##     v_x  = random.normalvariate(mu, sigma) * pos_neg_dir[random.randint(0, 1)]
##     v_y  = random.normalvariate(mu, sigma) * pos_neg_dir[random.randint(0, 1)]
##     v_mag = math.sqrt(v_x ** 2 + v_y ** 2)
##     v_dir  = math.atan2(v_y, v_x)
##
##     #    Calculate omega
##     #omega = (random.normalvariate(mu, sigma) / half_a) * pos_neg_dir[random.randint(0, 1)]
##     
##     #print('Start Velocity = [' + str(v_mag) + ', ' + str(v_dir) + ']')
##     #print('Start Omega = ' + str(omega))
##     return [v_mag, v_dir]


def normalised_start_velocity(T, m, total_no):
     k = 1.380649e-23 #  Boltzmann Constant
     #half_a = a/2
     mu =math.sqrt((k*T)/(m * total_no))
     #mu = 2 * math.sqrt((k * T)/m) * math.sqrt(2/(math.pi)) * math.sqrt((3 * k)/2)
     ##mu = 2 * math.sqrt((k * T)/m) * math.sqrt(2/(math.pi)) 
     #print('Mu = ' + str(mu))
     sigma = mu * 0.0001 #1/total_no * 
     #sigma = math.sqrt((k * T)/m * (3 * math.pi - 8)/math.pi) * math.sqrt((3 * k)/2)
     ##sigma = math.sqrt((k * T)/m * (3 * math.pi - 8)/math.pi)
     #print('Sigma' + str(sigma))
     pos_neg_dir = [-1, 1]

     v_x = (random.normalvariate(mu, sigma)) * pos_neg_dir[random.randint(0, 1)]
     v_y = (random.normalvariate(mu, sigma)) * pos_neg_dir[random.randint(0, 1)]
     
     #v_x  = math.sqrt((k*T)/(m * total_no)) * pos_neg_dir[random.randint(0, 1)]
     #v_y  = math.sqrt((k*T)/(m * total_no)) * pos_neg_dir[random.randint(0, 1)]
     v_mag = math.sqrt(v_x ** 2 + v_y ** 2)
     v_dir  = math.atan2(v_y, v_x)

     #    Calculate omega
     #omega = (random.normalvariate(mu, sigma) / half_a) * pos_neg_dir[random.randint(0, 1)]
     
     #print('Start Velocity = [' + str(v_mag) + ', ' + str(v_dir) + ']')
     #print('Start Omega = ' + str(omega))
     return [v_mag, v_dir]
##diatomic_molar_mass = 39.948 * 2           #       39.948 g/mol Argon
##diatomic_mass = diatomic_molar_mass * 6.022*(10**-23)
##m = diatomic_mass * 1/1000
##k = 1.380649e-23
##T = 200
##total_no = 12 ** 2
##mu =math.sqrt((2*k*T)/(m))
##sigma = mu * math.sqrt(1/2)
##test = scipy.stats.boltzmann.stats(mu, sigma)
##print(test)


def normalised_start_omega(T, m, a, total_no):
     k = 1.380649e-23 #  Boltzmann Constant       238
     #half_a = a/2
     reduced_mass = (np.float_(m) * np.float_(m))/(np.float_(m) + np.float_(m))
     I = reduced_mass * a ** 2
     mu = math.sqrt((k*T)/(I * total_no))
     sigma = mu * 0.0001 # 1/total_no * 
     #sigma = math.sqrt((k * T)/m * (3 * math.pi - 8)/math.pi)* math.sqrt((3 * k)/2)
     ##sigma = math.sqrt((k * T)/m * (3 * math.pi - 8)/math.pi)

     #    Calculate omega
     omega = (random.normalvariate(mu, sigma))
     #print('Omega = ' +str(omega))
     #omega = math.sqrt((k*T)/(I * total_no)) * pos_neg_dir[random.randint(0, 1)]
     
     return omega
     



##def average_start_velocity(temp, m, a, total_no):
##     #    calculate the velocity for the temperate
##     k = 1.380649e-23
##     omg_dir = [-1, 1]
##     v_x = 2 * math.sqrt((3*k*T)/(m * total_no)) * (omg_dir[random.randint(0, 1)])
##     v_y =math.sqrt((3*k*T)/(m * total_no)) * (omg_dir[random.randint(0, 1)])
##     #print(v_x)
##     #print(v_y)
##     ohm = (math.sqrt((3*k*T)/(m * total_no)) * (omg_dir[random.randint(0, 1)]))
##     #    calculate the mag of the linear velocity
##     v_mag = math.sqrt(v_x ** 2 + v_y ** 2)
##
##     v_dir  = (random.randint(0, 2 * 10 ** 3) * 10 ** - 3) * math.pi
##
##     return [v_mag, v_dir, ohm]


##a = 1e-10
##
##m = (2 * 6.63e-23) * 1/1000

##
##print(normalised_start_omega(300, m, a))
##print(normalised_start_omega(200, m, a))

##v = normalised_start_velocity(300, m, total_no)
##ohm = normalised_start_omega(300, m, a, total_no)
##
##I =  ((m/2) * ((a/2) ** 2)) * 2 
##
##EK = 0
##EK = EK + (1/2 * m * v[0] ** 2) + (1/2 * I * ohm ** 2)
##
##temp = (2 / (3 * k)) * EK




##
##print(average_start_velocity(300, m, a, total_no))
##print(average_start_velocity(200, m, a, total_no))
##print(average_start_velocity(100, m, a, total_no))
##print(average_start_velocity(50, m, a, total_no))
##print(average_start_velocity(25, m, a, total_no))
##print(average_start_velocity(10, m, a, total_no))
##print(average_start_velocity(1, m, a, total_no))
##print(average_start_velocity(0, m, a, total_no))
######
##print(normalised_start_velocity(300, m, total_no))
##print(normalised_start_velocity(200, m, total_no))
##print(normalised_start_velocity(100, m, total_no))
##print(normalised_start_velocity(25, m, total_no))
##print(normalised_start_velocity(10, m, total_no))
##print(normalised_start_velocity(1, m, total_no))
##print(normalised_start_velocity(0, m, total_no))
##
##print(normalised_start_omega(300, m, a, total_no))
##print(normalised_start_omega(T, m, a, total_no))
##print(normalised_start_omega(50, m, a, total_no))
