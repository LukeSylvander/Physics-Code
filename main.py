# Import maths library
import numpy as np
import math
import random





###     Import function from other files
from BoundaryFunctions import *
from BoltzmannMaxwell import *
from PlaceMolecules import *
from DistanceBetweenFunctions import *
from ForceFunctions import *

from MeasureSimulation import *

from SaveFunctions import *





### 2D simulations


### Setting Up the Box

#Set up the dimensions
DIM = 2
#Number of Particles
N = 64
#Box Size
box_size = 0.001 # meters
#Area of simulation
area_vol = box_size ** DIM
#density = N / area_vol

## Set up the box
box_x = [0, 10 ** -6]
box_y = [0, 10 ** -6]

diatomic_molar_mass = 39.948 * 2           #       39.948 g/mol Argon
diatomic_mass = diatomic_molar_mass/(6.022*(10**23))

### Constants
a = 1e-10 # m  Length of the molecule #         1e-10
#epsilon = 997   #      ϵ=0.997 kJ/mol 
#       r_0 = 1.2 * 10 ** (-10) # potential function has the value
m_1 = 6.63e-23 # kg  Mass of atom 1 #  Smaller mass smaller timestep
m_2 = 6.63e-23  #kg Mass of atom 2
#       Do not need reduced mass                        #m_total = (m_1 * m_2)/(m_1 + m_2) # Molecule mass
m_total = m_1 + m_2
### Make another function that adds up all the forces and the degree it came from
### make a Lennard-Jones Force equation




##def LJ_Force(r, epsilon, r_0):
##        """
##        Put the Lennard_Jones potential to calculate
##        the force of the interaction.
##
##        Inputs:
##        r - distance between two atoms
##        epsilon - Potential energy at the equilibrium bond length (eV)
##        r_0 - Distance at which the potential energy is zero
##
##        Outputs:
##        f - force between two atoms, this is a van der Waals interaction.  
##        """
##        f = (12 * epsilon) * (((r_0**12)/(np.float_(r)**13)) - (r_0**6)/(np.float_(r)**7))
##        #f = np.float_((12 * epsilon * (np.float_(r_0) ** 6) * ((np.float_(r_0) ** 6) - (np.float_(r) ** 6)))/(np.float_(r) ** 13))
##        return f
        
 

### Define torque of the molecule
def Torque (a, epsilon, r_0, r_a, r_b, alpha_a, alpha_b):
        
        """
        Calculate the torque for the equation.  

        Inputs:
        a - distance between the atoms in the molecule
        epsilon - Potential energy at the equilibrium bond length (eV)
        r_0 - Distance at which the potential energy is zero
        r_a - Position of atom a
        r_b - Position of atom b

        Outputs:
        t - torque of the molecule.  
        """
##        print('mol length ' + str(a))
##        print('epsilon ' + str(epsilon))
##        print('r_0 = ' + str(r_0))
##        print('r_a = ' + str(r_a))
##        print('r_b = ' + str(r_b))
##        print('alpha a = ' + str(alpha_a))
##        print('alpha a = ' + str(alpha_b))
        t = (a / 2) * (LJ_Force(r_a, epsilon, r_0) * math.sin(alpha_a) + LJ_Force(r_b, epsilon, r_0) * math.sin(alpha_b))
##        print('torque = ' + str(t))
        return t


def Inertia(m, a):
        
        """
        Calculate the torque for the equation.
        
        Inputs:
        a - distance between the atoms in the molecule
        epsilon - Potential energy at the equilibrium bond length (eV)
        m - mass of the molecule

        Outputs:
        I - Inertia of the diatomic molceule 
        """
        #I = ((m/2) * ((a/2) ** 2))  * 2
        mu = (np.float_(m) * np.float_(m))/(np.float_(m) + np.float_(m))
        I = mu * a ** 2
        return I


def find_r_and_phi_for_f(position_1, position_2):
        """
        Converts Cartisian coordinates of distance between atoms to polar coordinates 
        
        Inputs:
        Position 1 - cartisian coordinates of atom 1 [x, y]
        Position 2 - cartisian coordinates of atom 2 [x, y]

        Outputs:
        to_polar - is the polar coordinates between two atoms
        """
        dist_between = [(position_2[0] - position_1[0]), (position_2[1] - position_1[1])]
        to_polar = [(((dist_between[0] ** 2) + (dist_between[1] ** 2)) ** 0.5),
                    math.atan2(dist_between[1], dist_between[0])]
        return to_polar


def r_cm_new(r_old, v_old, dt, m_total, force):
        """
        Gives the new position of the molecule.

        Inputs:
        r_old - old position of the molecule
        v_old - old velocity of the molecule
        dt - change in time
        m_total - total mass
        epsilon - Potential energy at the equilibrium bond length (eV)
        r_0 - Distance at which the potential energy is zero
        r_a - Position of atom a
        r_b - Position of atom b

        Output:
        r_new - new position of the molecule

        """
        r_new = np.float_(r_old
                          ) + (np.float_(v_old) * np.float_(dt
                                                            )) + (np.float_(force)/np.float_(m_total
                                                                                             )) * (1/2 * np.float_(dt) ** 2)
        #print(v_old * dt + (force/ m_total) * 1 / 2 * dt ** 2)
        return r_new



def change_in_r(v_old, dt, m, force):
        """
        Gives the new position of the molecule.

        

        Inputs:
        v_old - old velocity of the molecule
        dt - change in time
        m_total - total mass
        epsilon - Potential energy at the equilibrium bond length (eV)
        r_0 - Distance at which the potential energy is zero

        Output:
        r_change - new position of the molecule

        """
##        print(v_old)
##        print(dt)
##        print(m)
##        print(force)
        change_r = np.float_(v_old
                             ) * np.float_(dt) + (np.float_(force)/np.float_(m)) * 1/2 * (np.float_(dt) ** (2))
        #change_r = v_old * dt + (force/m) * 1/2 * (dt ** (1/2))
        
        return change_r


def v_cm_new(v_old, dt, m_total, force):
        """
        Gives the new velocity of the molecule.

        Inputs:
        v_old - old velocity of the molecule
        dt - change in time
        m_total - total mass
        epsilon - Potential energy at the equilibrium bond length (eV)
        r_0 - Distance at which the potential energy is zero
        r_a - Position of atom a
        r_b - Position of atom b

        Output:
        v_new - new velocity of the molecule

        """
        v_new = np.float_(v_old) + ((np.float_(force) / np.float_(m_total)) * np.float_(dt))
        return v_new

def phi_new(phi_old, omega_old, dt, m, torque, a):
        """
        Gives the new phi of the molecule.

        Inputs:
        phi_old - old phi of the molecule
        omega_old - old omega of the molecule
        dt - change in time
        m_total - total mass
        epsilon - Potential energy at the equilibrium bond length (eV)
        r_0 - Distance at which the potential energy is zero
        r_a - Position of atom a
        r_b - Position of atom b

        Output:
        phi_new - new phi of the molecule



        """
        I = Inertia(m, a)

        phi_new = np.float_(phi_old) + np.float_(omega_old) * np.float_(dt) + ((np.float_(torque))/(np.float_(I))) * (1/2 * np.float_(dt) ** 2)
        
        if phi_new > (2 * math.pi):
                phi_ratio = phi_new/(2 * math.pi)
                phi_new = phi_new - ((2 * math.pi) * abs(math.floor(phi_ratio)))
        elif phi_new < 0:
                phi_ratio = phi_new/(2 * math.pi)
                phi_new = phi_new + ((2 * math.pi) *abs(math.floor(phi_ratio)))
        return phi_new 
 
def omega_new(omega_old, dt, m_total, torque, a):
        """
        Gives the new omega of the molecule.

        Inputs:
        omega_old - old omega of the molecule
        dt - change in time
        m_total - total mass
        epsilon - Potential energy at the equilibrium bond length (eV)
        r_0 - Distance at which the potential energy is zero
        r_a - Position of atom a
        r_b - Position of atom b

        Output:
        omega_new - new phi of the molecule

        """
        I = Inertia(m, a)
        
        omega_new = np.float_(omega_old) + (np.float_(torque)/(np.float_(I))) * np.float_(dt)
        return omega_new

def change_in_cart(r, phi): 
        """
        Convert polar coord. to cartisian coord.

        inputs:
        r - radius
        phi - angle

        Outputs:
        [x - position, y - position]

        """
        cart = [r * math.cos(phi), r * math.sin(phi)]
        return cart

def change_in_polar(x, y):
        polar = [(((x ** 2) + (y ** 2)) ** 0.5), math.atan2(y, x)]
        return polar

def atoms_postion(molecular_position, phi, a):
     """
     Atom Positions within the molecule 
     """
     [x, y] = molecular_position
     position_atom_1 = [np.float_(x) + np.float_(a/2) * np.float_(math.cos(phi - math.pi / 2)),
                        y + np.float_(a/2) * np.float_(math.sin(phi  - math.pi / 2))]
     position_atom_2 = [x + np.float_(a/2) * np.float_(math.cos(phi + math.pi / 2)),
                        y + np.float_(a/2) * np.float_(math.sin(phi  + math.pi / 2))]
     return [position_atom_1, position_atom_2]
"""
The order of the defined functions.
 - The centre of the molecules are origionally in cart coords.
 - Get the cart coords of the atoms in the molcules
 - Get the change in cart coods between the atoms
 - Convert the change of coords to polar
 - Calculate the new r, v, phi and omega
 - Convert polar cords to cart coords
 Polar
 ##### - r = [
  - phi = [0, (1 / 2 * math.pi)]
 Cartesian
  - x = [0, 0.001]
  - y = [0, 0.001]

"""

### Find the position and angle of atoms in the molecule

def find_atom_angle_and_position(a, r_new, phi_new):
        atom_angle = [phi_new - 1/2 * math.pi, phi_new + 1/2 * math.pi]
        position_a = [((a / 2) * math.cos(atom_angle[0])), ((a / 2) * math.sin(atom_angle[0]))]
        position_b = [((a / 2) * math.cos(atom_angle[1])), ((a / 2) * math.sin(atom_angle[1]))]
        return [position_a, atom_angle[0], position_b, atom_angle[1]]

### Find alpha, the angle between the atoms.
def find_alpha(atom_position_1, atom_position_2):
        """
        atom_position_1 = [x, y]
        atom_position_2 = [x, y]
        
        find_alpha([x_1, y_1], [x_2, y_2])  --> atan2(y_2 - y_1, x_2 - x_1]
        """
        [x_1, y_1] = atom_position_1
        [x_2, y_2] = atom_position_2
        alpha = math.atan2(y_2 - y_1, x_2 - x_1)
        return alpha


### Find r, the distance between atoms
def find_r(atom_position_1, atom_position_2):
        """
        atom_position_1 = [x, y]
        atom_position_2 = [x, y]
        
        find_alpha([x_1, y_1], [x_2, y_2])  -->  ((x_2 - x_1)**2 + (y_2 - y_1)**2) ** 0.5
        """
        [x_1, y_1] = atom_position_1
        [x_2, y_2] = atom_position_2
        r = (((x_2 - x_1) ** 2 + (y_2 - y_1) ** 2) ** 0.5)
        return r

 
def give_random():
        """
        give_random_position() --> [[x, y], phi]
        """
        x = ((random.randint(0, 10 ** 6)) * (10 ** - 9))
        y = ((random.randint(0, 10 ** 6)) * (10 ** - 9))
        phi = ((random.randint(0, 2 * 10 ** 6) * 10 ** - 9) * math.pi)
        return [[x, y], phi]

def give_random_molecule(total_number_molecules):
        """
            Import the number of molecules wanted in simulation
            give_random_molecule_position(10) -> molecule_positions with a length of 10

            Uses the give_random() function
        """
        i = 0
        molecule_position = []
        while i != total_number_molecules:
                molecule_position.append(give_random())
                i = i + 1
        return molecule_position




def give_atom_position(molecule_position_and_phi_list, a):
        """
                Input the list of molceules locations list and phi molecule list
                
                Output a list of atoms in the locations
        """
        i = 0
        atom_positions_in_molecule = []
        while i != len(molecule_position_and_phi_list):
                #print(find_atoms(molecule_position_and_phi_list[i][0], molecule_position_and_phi_list[i][1], a))
                atom_positions_in_molecule.append(find_atoms(molecule_position_and_phi_list[i][0], molecule_position_and_phi_list[i][1], a))
                i = i + 1
        return atom_positions_in_molecule






####
"""
                Note the place_position_grid is the function to the place the molecules, and record the major elements
                recorded in a list

                start code
                place_positions_grid(total_no, length_of_x_axis, length_of_y_axis, starting_phi, starting_velocity, a, rounded_to)


                Set up periodic boundary conditions

                def changes the cart to polar

                calculate the forces between the atoms
                        then add up all of the force with the direction at which it's applying

                        Calculate the torque acting on the molecule, need alpha for calculation

                calculate the kinematic equations and rotational equation to get new positions
                
                move the dt forward

                repeat for x amount of times

                make function to calculate instantanous pressure function

"""



"""
--------------------------------------------------------------------------------------------------------------------------------
Find the force function.


"""


##def sum_of_torque(list_of_distant, epsilon, r_0, a):
##     #    Max Counter
##     max_counter = len(list_of_distant[1][0])
##     #    Counter
##     i = 0
##     torque_int = 0
##     polar_dist_a = []
##     polar_dist_b = []
##     while max_counter > i:
##          polar_dist_a.append([list_of_distant[1][0][i], list_of_distant[1][1][i]])
##          polar_dist_b.append([list_of_distant[1][0][i], list_of_distant[2][1][i]])
##          torque_int = torque_int + Torque(a, epsilon, r_0,
##                                             polar_dist_a[i][0], polar_dist_b[i][0], polar_dist_a[i][1], polar_dist_b[i][1])
##          i = i + 1
##     #print('Sum of torque  = ' + str(torque_int))
##     return torque_int
     
##def unit_cell_graph(simulation_unit, n):
##     i = 0
##     x1 = []
##     y1 = []
##     x2 = []
##     y2 = []
##     while i < len(simulation_unit):
##          #print(simulation_unit[i][list_number])
##          x1.append(simulation_unit[i][1][0][0])
##          y1.append(simulation_unit[i][1][0][1])
##          x2.append(simulation_unit[i][1][1][0])
##          y2.append(simulation_unit[i][1][1][1])
##          
##          i = i + 1
##
##     plt.scatter(x1, y1)
##     plt.scatter(x2, y2)
##
##     #  plt.show()
##
##     plt.savefig(str(n) + ".png")
##
##     return None


##def work_of_cell(unit_data, epsilon, r_0):
##        #       Calculates the work of the unit
##        mol_position = []
##        i = 0
##
##        while i < len(unit_data):
##                mol_position.append([np.float_(unit_data[i][0][0]),np.float_(unit_data[i][0][1])])
##                i = i + 1
##
##        diff_r = []
##        
##        ii = 0
##        while ii < len(unit_data):
##                iii = 0
##                while iii < len(unit_data):
##                        if ii != iii:
##                                diff_r.append([mol_position[ii][0] - mol_position[iii][0],
##                                               mol_position[ii][1] - mol_position[iii][1]])
##                        iii = iii + 1
##                ii = ii + 1
##
##        iv = 0
##        mol_r = []
##        mol_phi = []
##        while iv < len(diff_r):
##                mol_r.append((diff_r[iv][0] ** 2 + diff_r[iv][1] ** 2) ** 0.5)
##                mol_phi.append(math.atan2(diff_r[iv][1], diff_r[iv][0]))
##                #print(iv)
##                iv = iv + 1
##
##        v = 0
##        mol_f = []
##
##        x_f = []
##        y_f = []
##
##        force_vector = []
##
##        while v < len(mol_r):
##                mol_f.append(LJ_Force(mol_r[v], epsilon, r_0))
##
##                x_f.append(mol_f[v] * math.cos(mol_phi[v]))
##                y_f.append(mol_f[v] * math.sin(mol_phi[v]))
##
##                force_vector.append([x_f[v], y_f[v]])
##                #print(v)
##
##                v = v + 1
##
##        vi = 0
##        w = 0
##        while vi < len(mol_r):
##                w = w + np.dot(mol_r, force_vector)
##                vi = vi + 1
##                #print(vi)
##        return w

     
"""
Let's put the functions together
-----------------------------------------------------------------------------------------------------------------------
        There are four main parts
        1. Set up simulation
        2. Energy Minimisation
                a. Give velocities corresponding to
                        T ~ 0.05 * T_target -> <k_0>
                b. Allow the system to equilibrate:
                        Will see some flow <U_0> -> <k_0>
        3. Equilibrium Proper
                a. Now re-set velocities to corresponding to target
                b. Allow the system to equilibrate.
        4. Run the molecular dynamics run -> "Production run"       

"""


class Diatomic2D:
        def __init__(self, total_no, length_x_axis, length_y_axis,
                     target_temp, a, m, epsilon, r_0, dt, total_repeats,
                     maximum_force, record_every_n, file_name):

                """#       Step 1 -        Set up simulation"""
                #       Start Temp
                start_temp = target_temp * 0.0005 # K

                boltz_con = 1.380649e-23
                #self.unit = place_positions_grid(total_no, (length_x_axis - length_x_axis/100), (length_y_axis - length_y_axis/100), start_velo, a)
                self.unit = place_positions_grid(total_no, length_x_axis, length_y_axis, a, m, start_temp)
                self.unit = periodic_boundary_conditions(self.unit, length_x_axis, "x", a)
                self.unit = periodic_boundary_conditions(self.unit, length_y_axis, "y", a)
                
                self.compare1 = []
                self.compare2 = []

                #       Change in dt Value
                self.origional_dt = dt
                self.check_time = 11
                
                #       Setup Counter
                n = 0
                self.total_dt = []
                self.sum_dt = 0
                file_counter = True

                #       Check Neihbours
                self.check_list = []

                #       Setup the pressure, temperature, kinetic energy
                self.KE = []
                self.U_P = []
                self.P = []
                self.temp = []
                self.P_per_temp = []
                self.area = length_x_axis * length_y_axis
                self.U_p = []
                                        #self.work = []
                self.total_energy = []
                self.change_U_p = 1000
                self.change_E_total = 2

                #       Find the different Phase of Simulation
                simulation_phase_list = ['Energy Minimisation',
                                         'Equilibrium Proper', 'Molecular Dynamics Run']
                
                self.phase_counter = simulation_phase_list[0]

                """#       Step 2 -        Energy Minimisation"""
                phase_counter_int = 0
                

                if self.phase_counter == 'Energy Minimisation':

                        """
        Change the while loop in the energy minimisation phase
                so that
                when the abs(U_p[-2]) - abs(U_p[-1]) < 1e-22


                        """
                        while n < 1000 and abs(self.change_U_p) > (0.01 * boltz_con * start_temp):
                                self.grid = periodic_boundary_outside(self.unit, length_x_axis, length_y_axis, a)
##                                £print('No of mol in grid = ' + str(len(self.grid)))

                                self.unit_int = find_specific_counter(self.grid, "Unit Cell")

                                #       Check if the boundary conditions are working
                                self.check = []
                                
                                n1 = 0
                                while n1 != len(self.unit_int):
                                        self.check.append(near_neighbours_int(self.grid, self.unit_int[n1], r_0 * 25))
                                        n1 = n1 + 1
                                self.check_list.append(self.check)

                                #Setup Counter i and introduce lists
                                i = 0
                                self.cdb = [] # cdb = dist_between
                                self.faw = [] # faw = force and work
                                self.work = []
                                self.sfl = []   # sfl = sum force list
                                self.stl = []   # stl = sum torque list
                                self.distance_between_unit = []
                                self.potential_energy = []

                                self.unit_dist_between = []
                                self.unit_work_n_force_total = []
                                self.unit_work = []
                                self.unit_force = []
                                self.unit_torque = []


                                #       Use of While loop in loop
                                while i < len(self.unit_int):
                                        self.cdb.append(distance_between_all_neighbours(self.grid, self.unit_int[i]))
                                        #       Distance between unit data
                                        self.distance_between_unit.append(dis_between_all_unit(self.unit, i))
                                        #print(self.cdb)
                                        self.faw.append(sum_of_force(self.cdb[i], epsilon, r_0))
                                        #print(sum_of_force(self.cdb[i], epsilon, r_0))
                                        self.sfl.append(self.faw[i][0])
                                        self.work.append(self.faw[i][1])
                                        #print(self.sfl)
                                        self.stl.append(sum_of_torque(self.cdb[i], epsilon, r_0, a))

                                        self.unit_dist_between.append(distance_between_all_neighbours(self.unit,
                                                                                                      i))
                                        self.unit_work_n_force_total.append(
                                                sum_of_force(self.unit_dist_between[i], epsilon, r_0))
                                        self.unit_force.append(self.unit_work_n_force_total[i][0])
                                        self.unit_work.append(self.unit_work_n_force_total[i][1])
                                        self.unit_torque.append(
                                                sum_of_torque(self.unit_dist_between[i], epsilon, r_0, a))

                                        #       Calculate the potential energy of the unit
                                        self.potential_energy.append(
                                                sum_potential_energy(self.unit_dist_between[i],
                                                                     epsilon, r_0))
                                                                            
                                        
                                        i = i + 1

##                                print('x  ' + str(self.cdb[1][0][0][0]) + ', y  '+ str(self.cdb[1][0][0][1]))
##                                print('Force  ' + str(self.sfl[-1]) + 'Torque  '+ str(self.stl[-1]))
                                i = 0
                                #       Potential Energy

                                
                                [dt, self.check_time] = check_force(dt, self.origional_dt,
                                                                    self.sfl, self.check_time, maximum_force)

                                i = 0
                                self.pp = []    # pp = polar position
                                self.cv = []    # cv = cart velocity
                                
                                while i < len(self.unit_int):
                                        self.pp.append(change_in_polar(self.grid[self.unit_int[i]][0][0],
                                                                       self.grid[self.unit_int[i]][0][1]))
                                        self.cv.append(change_in_cart(self.grid[self.unit_int[i]][2][0],
                                                                      self.grid[self.unit_int[i]][2][1]))
                                        i = i + 1    
                                
                                
                                #       Setip counter and introduce list and constats
                                i = 0

                                self.nrx = []   #       nrx = new r x
                                self.nvx = []   #       nvx = new v x
                                
                                self.nry = []   #       nry = new r y
                                self.nvy = []   #       nvy = new v y

                                self.nphi = []  #       nphi = new phi
                                self.nom = []   #       nom = now omega
                                self.nap = []   #       nap = new atom position

                                self.nrm = []   #       nrp = new r mol

                                self.crx = []   #       change r x
                                self.cry = []   #       change r y
                                self.cr = []   #       change r

                                self.unit_v_x = []       #         unit velocity x
                                self.unit_v_y = []       #         unit velocity y
                                self.unit_omega = []


                                #print("The mass is " + str(m))

                                while i < len(self.unit_int):
                                        #print(self.unit[i][0][0])
                                        #print(self.cv[i][0])
                                        #print(self.faw[i])
                                        #print(self.sfl[i])
                                        self.nrx.append(r_cm_new(self.unit[i][0][0], self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.nry.append(r_cm_new(self.unit[i][0][1], self.cv[i][1], dt, m, self.sfl[i][1]))
                                        self.nmolpos = [self.nrx[i], self.nry[i]]                       

                                        self.nvx.append(v_cm_new(self.cv[i][0], dt, m, self.sfl[i][0]))
                                        #print(self.sfl[i])
                                        #print(self.cv[i][0])
                                        self.nvy.append(v_cm_new(self.cv[i][1], dt, m, self.sfl[i][1]))

                                        self.nphi.append(phi_new(self.unit[i][3], self.unit[i][4], dt, m, self.stl[i], a))
                                        self.nom.append(omega_new(self.unit[i][4], dt, m, self.stl[i], a))

                                        self.nap.append(atoms_postion([self.nrx[i], self.nry[i]], self.nphi[i], a))

                                        self.nrm.append([self.nrx[i], self.nry[i]])

                                        self.crx.append(change_in_r(self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.cry.append(change_in_r(self.cv[i][1], dt, m, self.sfl[i][1]))
                                        self.cr.append([self.crx[i], self.cry[i]])

##                                        self.work.append(np.dot(self.cr[i], self.sfl[i]))

                                        self.unit_v_x.append(v_cm_new(self.cv[i][0], dt,
                                                                      m, self.unit_force[i][0]))
                                        #print(self.unit_v_x[i])
                                        self.unit_v_y.append(v_cm_new(self.cv[i][1], dt,
                                                                      m, self.unit_force[i][1]))

                                        self.unit_omega.append(omega_new(self.unit[i][4], dt,
                                                                         m, self.unit_torque[i], a))
                                        #print(self.unit_omega[i])
                                        
                                        

                                        i = i + 1
                        
                                #       Setip Counter
                                i = 0#

                                #

                                while i < len(self.unit_int):
                                        self.unit[i][0] = [self.nrx[i], self.nry[i]]
                                        self.unit[i][1] = self.nap[i]
                                        self.unit[i][2] = change_in_polar(0, 0)
                                        self.unit[i][3] = self.nphi[i]
                                        self.unit[i][4] = 0
                                        i = i + 1

                                #       Check if in boundary conditions
                                self.unit = periodic_boundary_conditions(self.unit, length_x_axis, "x", a)
                                self.unit = periodic_boundary_conditions(self.unit, length_y_axis, "y", a)

                                self.grid = periodic_boundary_outside(self.unit, length_x_axis, length_y_axis, a)

                                #print('cbd = ' + str(self.cdb[0][0]))
                                #print('sfl = ' + str(self.sfl[0]))

                                #       Calculates the work of the cell
                                self.unit_work = sum(self.unit_work)
                                print('Work ' + str(self.unit_work))

                                #       Calculate the potential energy
                                
                                self.U_p.append(sum(self.potential_energy))
                                print("Potential Energy of the System = " + str(self.U_p[-1]))

                                #       kinetic energy, temp, pressure
                                self.KE.append(np.float_(
                                        unit_diatomic_kinetic_energy([self.unit_v_x, self.unit_v_y],
                                                                     self.unit_omega, m, a)))
                                print("KE of the System: " + str(np.float_(self.KE[-1])))

                                #       Total ENergy
                                self.total_energy.append(abs(self.U_p[-1]) + abs(self.KE[-1]))
                                print("Total Energy of System " + str(self.total_energy[-1]))
                                
                                self.temp.append(np.float_(localised_temperature(self.KE[n])))
                                print("Temperature of the System: " + str(np.float_(self.temp[-1])))
                                self.P.append(np.float_(localised_pressure(self.temp[n], self.area,
                                                                           self.unit_work, total_no)))
                                                                           
                                print("Pressure of the System: " + str(self.P[-1]))


                                if n == 0:
                                        self.total_dt.append(dt)
                                else:
                                        self.total_dt.append(np.float_(self.total_dt[n-1]) + np.float_(dt))

                                self.sum_dt = self.sum_dt + np.float_(dt)
                                
                        
                                #print(str(self.compare1) + " = " + str(self.compare2))

                                if n > 3:
                                        self.change_U_p = (abs(np.float_(self.U_p[-2])) - abs(np.float_(self.U_p[-1])))
                                        self.change_E_total = abs(np.float_(self.total_energy[-2])) - abs(np.float_(self.total_energy[-1]))
                                        print('Change in U_p = ' + str(self.change_U_p))

                                n = n + 1
                #Avererage energy temp


                """#       Step 3 -        Setup Equilibrium proper"""

                nn = 0
                print(str(self.phase_counter) + ' is completed.')

                self.phase_counter = simulation_phase_list[1]

                self.change_U_p = 1000

                self.change_E_total = 2

                self.sum_dt = 0

                #       Change unit velocity
                self.unit = change_velocity_keep_positions(self.unit, target_temp, m, a)

                if self.phase_counter == 'Equilibrium Proper':
                        #
                        while nn < 1000 and abs(self.change_U_p) > 0.001 * boltz_con * target_temp:
                                self.grid = periodic_boundary_outside(self.unit, length_x_axis, length_y_axis, a)

                                self.unit_int = find_specific_counter(self.grid, "Unit Cell")

                                #       Check if the boundary conditions are working
                                self.check = []
                                
                                n1 = 0
                                while n1 != len(self.unit_int):
                                        self.check.append(near_neighbours_int(self.grid, self.unit_int[n1], r_0 * 25))
                                        n1 = n1 + 1
                                self.check_list.append(self.check)

                                #Setup Counter i and introduce lists
                                i = 0
                                self.cdb = [] # cdb = cart_dist_between
                                self.faw = [] # faw = force and work
                                self.work = []
                                self.sfl = []   # sfl = sum force list
                                self.stl = []   # stl = sum torque list
                                self.distance_between_unit = []
                                self.potential_energy = []

                                #       Use of While loop in loop
                                while i < len(self.unit_int):
                                        self.cdb.append(distance_between_all_neighbours(self.grid, self.unit_int[i]))
                                        #       Distance between unit data
                                        self.distance_between_unit.append(dis_between_all_unit(self.unit, i))
                                        #print(self.cdb)
                                        self.faw.append(sum_of_force(self.cdb[i], epsilon, r_0))
                                        #print(sum_of_force(self.cdb[i], epsilon, r_0))
                                        self.sfl.append(self.faw[i][0])
                                        self.work.append(self.faw[i][1])
                                        #print(self.sfl)
                                        self.stl.append(sum_of_torque(self.cdb[i], epsilon, r_0, a))

                                        #       Calculate the potential energy of the unit
                                        self.potential_energy.append(sum_potential_energy(self.cdb[i],
                                                                                          epsilon, r_0))
                                        
                                        i = i + 1
                                i = 0
                                #       Potential Energy

                                
                                [dt, self.check_time] = check_force(dt, self.origional_dt,
                                                                    self.sfl, self.check_time, maximum_force)

                                i = 0
                                self.pp = []    # pp = polar position
                                self.cv = []    # cv = cart velocity
                                
                                while i < len(self.unit_int):
                                        self.pp.append(change_in_polar(self.grid[self.unit_int[i]][0][0],
                                                                       self.grid[self.unit_int[i]][0][1]))
                                        self.cv.append(change_in_cart(self.grid[self.unit_int[i]][2][0],
                                                                      self.grid[self.unit_int[i]][2][1]))
                                        i = i + 1    
                                
                                
                                #       Setip counter and introduce list and constats
                                i = 0

                                self.nrx = []   #       nrx = new r x
                                self.nvx = []   #       nvx = new v x
                                
                                self.nry = []   #       nry = new r y
                                self.nvy = []   #       nvy = new v y

                                self.nphi = []  #       nphi = new phi
                                self.nom = []   #       nom = now omega
                                self.nap = []   #       nap = new atom position

                                self.nrm = []   #       nrp = new r mol

                                self.crx = []   #       change r x
                                self.cry = []   #       change r y
                                self.cr = []   #       change r

                                #print("The mass is " + str(m))

                                while i < len(self.unit_int):
                                        #print(self.unit[i][0][0])
                                        #print(self.cv[i][0])
                                        #print(self.faw[i])
                                        #print(self.sfl[i])
                                        self.nrx.append(r_cm_new(self.unit[i][0][0], self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.nry.append(r_cm_new(self.unit[i][0][1], self.cv[i][1], dt, m, self.sfl[i][1]))
                                        self.nmolpos = [self.nrx[i], self.nry[i]]                       

                                        self.nvx.append(v_cm_new(self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.nvy.append(v_cm_new(self.cv[i][1], dt, m, self.sfl[i][1]))

                                        self.nphi.append(phi_new(self.unit[i][3], self.unit[i][4], dt, m, self.stl[i], a))
                                        self.nom.append(omega_new(self.unit[i][4], dt, m, self.stl[i], a))

                                        self.nap.append(atoms_postion([self.nrx[i], self.nry[i]], self.nphi[i], a))

                                        self.nrm.append([self.nrx[i], self.nry[i]])

                                        self.crx.append(change_in_r(self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.cry.append(change_in_r(self.cv[i][1], dt, m, self.sfl[i][1]))
                                        self.cr.append([self.crx[i], self.cry[i]])


                                        #       Work
##                                        self.work.append(math.dot(self.cr[i], self.sfl[i]))

                                        i = i + 1
                        
                                #       Setip Counter
                                i = 0

                                while i < len(self.unit_int):
                                        self.unit[i][0] = [self.nrx[i], self.nry[i]]
                                        self.unit[i][1] = self.nap[i]
                                        self.unit[i][2] = change_in_polar(self.nvx[i], self.nvy[i])
                                        self.unit[i][3] = self.nphi[i]
                                        self.unit[i][4] = self.nom[i]
                                        i = i + 1

                                #       Check if in boundary conditions
                                self.unit = periodic_boundary_conditions(self.unit, length_x_axis, "x", a)
                                self.unit = periodic_boundary_conditions(self.unit, length_y_axis, "y", a)

                                self.grid = periodic_boundary_outside(self.unit, length_x_axis, length_y_axis, a)

                                #print('cbd = ' + str(self.cdb[0][0]))
                                #print('sfl = ' + str(self.sfl[0]))

                                #       Calculates the work of the cell
                                self.work = sum(self.work)
                                print('Work ' + str(self.work))

                                #       Calculate the potential energy
                                
                                self.U_p.append(sum(self.potential_energy))
                                #print("Potential Energy of the System = " + str(self.U_p[-1]))

                                #       kinetic energy, temp, pressure
                                self.KE.append(np.float_(diatomic_kinetic_energy(self.unit, m, a)))

                                #       Total ENergy
                                self.total_energy.append(abs(self.U_p[-1]) + abs(self.KE[-1]))
                                print("Total Energy of System " + str(self.total_energy[-1]))
                                
##                                print("KE of the System: " + str(np.float_(self.KE[-1])))
                                self.temp.append(np.float_(localised_temperature(self.KE[n])))
                                print("Temperature of the System: " + str(np.float_(self.temp[-1])))
                                self.P.append(np.float_(localised_pressure(self.temp[n], self.area,
                                                                           self.work, total_no)))
                                                                           
                                print("Pressure of the System: " + str(self.P[-1]))
                                
                                self.P_per_temp.append(np.float_(self.P[n]) / self.temp[n])

##                                print("Pressure/Temp = " + str(self.P_per_temp[-1]))



##                                print(len(self.cdb))

                                if nn == 0:
                                        self.total_dt.append(dt)
                                else:
                                        self.total_dt.append(np.float_(self.total_dt[n-1]) + np.float_(dt))

                                self.sum_dt = self.sum_dt + np.float_(dt)
                                
                                self.compare1.append(self.area * self.P[-1])
                                self.compare2.append(self.temp[-1] * (total_no * 1.66054e-24) * 8.314462)
                                
                                print(str(self.compare1[-1]) + " = " + str(self.compare2[-1])) 

                                nn = nn + 1
                                if nn > 3:
                                        self.change_U_p = abs(np.float_(self.U_p[-2])) - abs(np.float_(self.U_p[-1]))
                                        self.change_E_total = abs(np.float_(self.total_energy[-2])) - abs(np.float_(self.total_energy[-1]))
                                        print('Change in U_p = ' + str(self.change_U_p))

############################################################################


                
                #       Step 4 -        Molecular Dynamics Run 

                self.phase_counter = simulation_phase_list[2]

                self.sum_dt = 0

                nnn = 0
                

                if self.phase_counter == 'Molecular Dynamics Run':
                        #       Start writing csv
                        file_counter = save_data_csv(file_name, self.unit, self.P,
                                                     self.temp, self.total_dt, file_counter)
                        #Use While Loop
                        while nnn != total_repeats:
                                self.grid = periodic_boundary_outside(self.unit, length_x_axis, length_y_axis, a)

                                self.unit_int = find_specific_counter(self.grid, "Unit Cell")

                                #       Check if the boundary conditions are working
                                self.check = []
                                
                                n1 = 0
                                while n1 != len(self.unit_int):
                                        self.check.append(near_neighbours_int(self.grid, self.unit_int[n1], r_0 * 25))
                                        n1 = n1 + 1
                                self.check_list.append(self.check)

                                #Setup Counter i and introduce lists
                                i = 0
                                self.cdb = [] # cdb = cart_dist_between
                                self.faw = [] # faw = force and work
                                self.work = []
                                self.sfl = []   # sfl = sum force list
                                self.stl = []   # stl = sum torque list
                                self.distance_between_unit = []
                                self.potential_energy = []

                                self.unit_work_n_force =[]
                                self.unit_work = []
                                self.unit_force = []
                                self.unit_torque = []

                                #       Use of While loop in loop
                                while i < len(self.unit_int):
                                        self.cdb.append(distance_between_all_neighbours(self.grid, self.unit_int[i]))
                                        #       Distance between unit data
                                        self.distance_between_unit.append(dis_between_all_unit(self.unit, i))
                                        #print(self.cdb)
                                        self.faw.append(sum_of_force(self.cdb[i], epsilon, r_0))
                                        #print(sum_of_force(self.cdb[i], epsilon, r_0))
                                        self.sfl.append(self.faw[i][0])
                                        self.work.append(self.faw[i][1])
                                        #print(self.sfl)
                                        self.stl.append(sum_of_torque(self.cdb[i], epsilon, r_0, a))

                                        #       Calculate the potential energy of the unit
                                        self.potential_energy.append(sum_potential_energy(self.cdb[i],
                                                                                          epsilon, r_0))


                                        self.unit_work_n_force_total.append(
                                                sum_of_force(self.unit_dist_between[i], epsilon, r_0))
                                        self.unit_force.append(self.unit_work_n_force_total[i][0])
                                        self.unit_work.append(self.unit_work_n_force_total[i][1])
                                        self.unit_torque.append(
                                                sum_of_torque(self.unit_dist_between[i], epsilon, r_0, a))
                                        
                                        i = i + 1
                                i = 0
                                #       Potential Energy

                                
                                [dt, self.check_time] = check_force(dt, self.origional_dt,
                                                                    self.sfl, self.check_time, maximum_force)

                                i = 0
                                self.pp = []    # pp = polar position
                                self.cv = []    # cv = cart velocity
                                
                                while i < len(self.unit_int):
                                        self.pp.append(change_in_polar(self.grid[self.unit_int[i]][0][0],
                                                                       self.grid[self.unit_int[i]][0][1]))
                                        self.cv.append(change_in_cart(self.grid[self.unit_int[i]][2][0],
                                                                      self.grid[self.unit_int[i]][2][1]))
                                        i = i + 1    
                                
                                
                                #       Setip counter and introduce list and constats
                                i = 0

                                self.nrx = []   #       nrx = new r x
                                self.nvx = []   #       nvx = new v x
                                
                                self.nry = []   #       nry = new r y
                                self.nvy = []   #       nvy = new v y

                                self.nphi = []  #       nphi = new phi
                                self.nom = []   #       nom = now omega
                                self.nap = []   #       nap = new atom position

                                self.nrm = []   #       nrp = new r mol

                                self.crx = []   #       change r x
                                self.cry = []   #       change r y
                                self.cr = []   #       change r

                                self.unit_v_y = []
                                self.unit_v_x = []
                                self.unit_omega = []

                                #print("The mass is " + str(m))

                                while i < len(self.unit_int):
                                        #print(self.unit[i][0][0])
                                        #print(self.cv[i][0])
                                        #print(self.faw[i])
                                        #print(self.sfl[i])
                                        self.nrx.append(r_cm_new(self.unit[i][0][0], self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.nry.append(r_cm_new(self.unit[i][0][1], self.cv[i][1], dt, m, self.sfl[i][1]))
                                        self.nmolpos = [self.nrx[i], self.nry[i]]                       

                                        self.nvx.append(v_cm_new(self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.nvy.append(v_cm_new(self.cv[i][1], dt, m, self.sfl[i][1]))

                                        self.nphi.append(phi_new(self.unit[i][3], self.unit[i][4], dt, m, self.stl[i], a))
                                        self.nom.append(omega_new(self.unit[i][4], dt, m, self.stl[i], a))

                                        self.nap.append(atoms_postion([self.nrx[i], self.nry[i]], self.nphi[i], a))

                                        self.nrm.append([self.nrx[i], self.nry[i]])

                                        self.crx.append(change_in_r(self.cv[i][0], dt, m, self.sfl[i][0]))
                                        self.cry.append(change_in_r(self.cv[i][1], dt, m, self.sfl[i][1]))
                                        self.cr.append([self.crx[i], self.cry[i]])


                                        self.unit_v_x.append(v_cm_new(self.cv[i][0], dt,
                                                                      m, self.unit_force[i][0]))
                                        #print(self.unit_v_x[i])
                                        self.unit_v_y.append(v_cm_new(self.cv[i][1], dt,
                                                                      m, self.unit_force[i][1]))

                                        self.unit_omega.append(omega_new(self.unit[i][4], dt,
                                                                         m, self.unit_torque[i], a))

##                                        self.work.append(math.dot(self.cr[i], self.sfl[i]))

                                        i = i + 1
                        
                                #       Setip Counter
                                i = 0

                                while i < len(self.unit_int):
                                        self.unit[i][0] = [self.nrx[i], self.nry[i]]
                                        self.unit[i][1] = self.nap[i]
                                        self.unit[i][2] = change_in_polar(self.nvx[i], self.nvy[i])
                                        self.unit[i][3] = self.nphi[i]
                                        self.unit[i][4] = self.nom[i]
                                        i = i + 1

                                #       Check if in boundary conditions
                                self.unit = periodic_boundary_conditions(self.unit, length_x_axis, "x", a)
                                self.unit = periodic_boundary_conditions(self.unit, length_y_axis, "y", a)

                                self.grid = periodic_boundary_outside(self.unit, length_x_axis, length_y_axis, a)

                                #       Calculates the work of the cell
                                self.unit_work = sum(self.unit_work)
                                print('Work ' + str(self.unit_work))

                                #       Calculate the potential energy
                                
                                self.U_p.append(sum(self.potential_energy))
                                print("Potential Energy of the System = " + str(self.U_p[-1]))

                                #       kinetic energy, temp, pressure
                                self.KE.append(np.float_(
                                        unit_diatomic_kinetic_energy([self.unit_v_x, self.unit_v_y],
                                                                     self.unit_omega, m, a)))
                                print("KE of the System: " + str(np.float_(self.KE[-1])))

                                #       Total ENergy
                                self.total_energy.append(abs(self.U_p[-1]) + abs(self.KE[-1]))
                                print("Total Energy of System " + str(self.total_energy[-1]))
                                
                                self.temp.append(np.float_(np.float_(localised_temperature(self.KE[n]))))
                                print("Temperature of the System: " + str(np.float_(self.temp[-1])))
                                self.P.append(np.float_(localised_pressure(self.temp[n], self.area,
                                                                           self.unit_work, total_no)))
                                                                           
                                print("Pressure of the System: " + str(self.P[-1]))

####                                #       Calculates the work of the cell
####                                self.work = sum(self.work)
####                                print('Work ' + str(self.work))
####
####                                #       Calculate the potential energy
####                                
####                                self.U_p.append(sum(self.potential_energy))
####                                #print("Potential Energy of the System = " + str(self.U_p[-1]))
####
####                                #       kinetic energy, temp, pressure
####                                self.KE.append(np.float_(diatomic_kinetic_energy(self.unit, m, a)))
####
####                                #       Total ENergy
####                                self.total_energy.append(abs(self.U_p[-1]) + abs(self.KE[-1]))
####                                print("Total Energy of System " + str(self.total_energy[-1]))
####                                
####                                #print("KE of the System: " + str(np.float_(self.KE[-1])))
####                                self.temp.append(np.float_(localised_temperature(self.KE[-1])))
####                                #print("Temperature of the System: " + str(np.float_(self.temp[-1])))
####                                self.P.append(np.float_(localised_pressure(self.temp[-1], self.area,
####                                                                           self.work, total_no)))
####                                                                           
####                                print("Pressure of the System: " + str(self.P[-1]))
                                
                                self.P_per_temp.append(np.float_(self.P[n]) / self.temp[n])

                                #print("Pressure/Temp = " + str(self.P_per_temp[-1]))



                                print(len(self.cdb))

                                if nnn == 0:
                                        self.total_dt.append(dt)
                                else:
                                        self.total_dt.append(np.float_(self.total_dt[n-1]) + np.float_(dt))

                                self.sum_dt = np.float_(self.sum_dt + np.float_(dt))
                                
                                self.compare1.append(self.area * self.P[-1])
                                
                                self.compare2.append(self.temp[-1] * (total_no * 1.66054e-24) * 8.314462)
                                                #print(str(self.compare1) + " = " + str(self.compare2)) 

                                nnn = nnn + 1
                                print('Cycle Number Done, ' + str(nnn))

                                if nnn % record_every_n == 0:
                                        file_counter = save_data_csv(file_name, self.unit, self.P[-1],
                                                                     self.temp[-1], self.sum_dt, file_counter)


#       import time


                #unit_cell_graph(self.unit, 0, "Recent Simulation Unit Positions")
                #       Ideal Gas Law
                
"""
-------------------------------_______________________-----------------------------------
"""
        #       Diameter of argon is 0.4 nm
        #       

                        #1.2 * 10 ** (-10) # van der Waal radai 
#t_0 = time.clock()


no_of_particles = 100 #  0#12 ** 2#12 ** 2
epsilon  = 997  * (1.66054e-24 )
r_0 = (2 ** (1/6)) * (3.4 * 10 ** (-10))                 
a = 1e-10
m_kg = diatomic_mass/1000
#sim_boundary =  4e-9
dt = np.float_(1e-2015)
total_cycles = 1000
force_max = 1e-6
store_cycle = 1# 10 - 100 cycles
file_name = '100ParticleDiatomic2D_at_100K_5pd.csv'
start_T = 100 #K
m = m_kg

total_mass = no_of_particles * m_kg
density = 4e-5
sim_boundary = math.sqrt(total_mass/density)
sigma = (r_0/((2 ** (1/6))))

length_molecule_x = (a + sigma * 2) * math.sqrt(no_of_particles)

length_molecule_y = (2 * sigma) * math.sqrt(no_of_particles)

sb_x = (length_molecule_x * (100)/100)
sb_y = (length_molecule_y * (100)/100)

#print(sim_boundary)

##Dia2Dsim = Diatomic2D(no_of_particles, sb_x, sb_y, start_T, a, m_kg, epsilon,
##                      r_0, dt, total_cycles, force_max, store_cycle, file_name)

