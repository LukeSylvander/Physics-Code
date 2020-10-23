###       Correlation Functions
import math
import numpy as np

#    average function
def average(number_list):
     total_sum = sum(number_list)
     total_len = len(number_list)
     average_of_list = np.float_(total_sum/total_len)
     return average_of_list


#    Velocity autocorrelation function

def velocity_autocorrelation(initial_velocity, current_velocity):
     #    v_intial = dot product of self
     i = 0
     max_counter = len(initial_velocity)
     v_initial_dot = []
     v_change_dot = []

     while i < max_counter:
          v_initial_dot.append(np.dot(initial_velocity[i], initial_velocity[i]))
          v_change_dot.append(np.dot(current_velocity[i], initial_velocity[i]))
          i += 1

     velo_auto = (average(v_change_dot))/(average(v_initial_dot))
     return velo_auto

#    Gaussian Rotational Correlation Function

def gaussian_rotational_correlation(inital_phi, current_phi):

     i = 0
     max_list_phi = len(inital_phi)
     phi_change_dot = []

     while i < max_list_phi:
          phi_change_dot.append(current_phi[i] -inital_phi[i])
          i += 1

     gaussian_rotation_corr = np.cos(np.float_(average(phi_change_dot)))
     #         exp(-2Dt)

     return gaussian_rotation_corr


#    Legendre Polynomials

def legendre_polynomials(inital_phi, current_phi, polynomial_no):

     i = 0
     max_list_phi = len(inital_phi)
     #print(max_list_phi)
     phi_change_dot = []

     while i < max_list_phi:
          phi_change_dot.append(current_phi[i] - inital_phi[i])
          i += 1
     #    theta
     theta = average(phi_change_dot)

     x = np.cos(theta)

     if polynomial_no == 2:
          #    2nd legendre_polynomial
          legendre_polynomial = 1/2 * (3 * x ** 2 - 1)
          #         == exp(-6 * D * t)

     if polynomial_no == 3:
          #    3rd legendre_polynomial
          legendre_polynomial = 1/2 * (5 * x ** 3 - 3 * x)
          #         == exp(-12 * D *t)

     if polynomial_no == 4:
          #    4th legendre_polynomial
          legendre_polynomial = 1/8 * (35 * x ** 4 - 30 * x ** 2 + 3)
          #         == exp(-20 * D *t)
     if polynomial_no == 5:
          #    5th legendre_polynomial
          legendre_polynomial = 1/8 * (63 * x ** 5 - 70 * x ** 3 + 15 * x)
          #         == exp(-30 * D * t)

     return legendre_polynomial






     

     
          
          
