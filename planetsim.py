"""
PlanetSim, a project to simulate astronomical bodies
Douglas Hull & Jack Manners
2/3/20
"""

import time
from particle3D import Particle3D as P3D
import numpy as np
import sys

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Value of G taken from CODATA 2018
G=6.67430e-11

class Planet(P3D):
    """
    Planet class, inherits the Particle3D class
    with additional planetary properties.

    Properties:
    primary(int) - index of primary body planet orbits around
    apo/per(float) - apoapsis/periapsis of planet
    """

    template = {'label':(str,1), 'position':(float,3), 'velocity':(float,3), 'mass':(float,1), 'primary':(int,1)}

    def __init__(self, label, pos, vel, mass, pri):
        super().__init__(label, pos, vel, mass)
        self.primary = pri
        self.apo = self.per = 0


class PlanetSim(object):
    """
    Class represents a planetary simulation.

    Properties:
    dt(float) - timestep (s)
    sim_time(float) - duration of simulation (s)
    t0(float) - initial time (s)
    time(float) - current time (s)
    forces(np.array[N,3](float) - current interplanetary forces
    """

    def __init__(self, planet_file, param_file):
        """
        Create simulation using data from files

        :param planet_file: file containing data for all planets
        :param param_file:  file containing simulation parameters
        """
        self.planet_list = Planet.from_file(planet_file)
        self.dt, self.sim_time, self.t0 = self.read_params(param_file)
        self.time = self.t0
        self.forces = np.zeros((len(self.planet_list),3))


    def read_params(self,param_file):
        """
        Read param file, extract parameters and return

        File in format: # <timestep> <duration> <initial time>

        :param param_file: parameter file name (string)
        """
        with open(param_file, 'r') as f:
            for line in f:
                if line[0] != '#':
                    #skip lines that don't start with '#'
                    continue
                #cuts line into list and drops '# char
                params = line.split()[1:]
                #convert params to floats
                params = [float(x) for x in params]
        f.close()
        return params

    def CoM_fix(self):
        """
        Correct for centre-of-mass motion of planets

        CoM has velocity V = P/M,
        where P is the total momentum of the planets, and M is the total mass.
        V is subtracted from all bodies.
        """
        P = np.zeros(3)
        M = 0
        for p in self.planet_list:
            P += (p.mass * p.velocity)
            M += p.mass
        v_com = P/M
        for p in self.planet_list:
            p.velocity -= v_com

    def force(self, p1, p2):
        """
        Returns force between two planets due to gravitational attraction.

        For planets p1 & p2, with masses m1,m2:
        F = G * m1 * m2 / R_12**2 , with direction r_12
        where r_12 is the vector from p1 to p2, R_12 is it's magnitude,
        and G is the gravitational constant.
        """

        r = p2.position - p1.position
        R = np.linalg.norm(r)
        return r*G*p1.mass*p2.mass*(R**-3)

    def potential(self, p1, p2):
        """
        Returns potential between two planets due to gravitational attraction.

        For planets p1 & p2, with masses m1,m2:
        U = -G * m1 * m2 / R_12**2
        where R_12 is the distance between them and G is the gravitational constant.
        """
        r = p2.position - p1.position
        R = np.linalg.norm(r)
        return -G*p1.mass*p2.mass/R

    def update_forces(self):
        #Calculate force interactions between every pair
        n = len(self.planet_list)
        self.forces = np.zeros((n,3))
        for i,p_i in enumerate(self.planet_list):
            for j,p_j in enumerate(self.planet_list):
                if j > i:
                    f = self.force(p_i, p_j)
                    self.forces[i] += f
                    self.forces[j] -= f

    def update_all(self, method, forces):
        #Calls 'p.method(dt,force[p])' for each planet p
        #Used to call time-integration methods on all planets
        for i,p_i in enumerate(self.planet_list):
            getattr(p_i, method)(self.dt, forces[i])

    def calc_apsides(self):
        # Updates the apo/peri-apsis if the planet's orbital distance
        # is greater/lesser the stored values.
        for p in self.planet_list:
            r = np.linalg.norm(p.position - self.planet_list[p.primary].position)
            if r < p.per:
                p.per = r
            elif r > p.apo:
                p.apo = r

    def count_period(self, pos):
        # Calculate the period by counting the minima of the orbital distance.
        # A minimum is found as a point lower than the points on either side.
        # The period T = dt * (steps_for_N_orbits / N) / (seconds/day)
        minima_steps = []
        T=0
        for i in range(len(pos)-2):
            if pos[i] - pos[i+1] > 0:
                if pos[i+2] - pos[i+1] > 0:
                    minima_steps.append(i+1)
        N = len(minima_steps)
        if N > 0:
            T = self.dt*((minima_steps[-1])/N)/(60*60*24)
        return T

    def estimated_period(self,p):
        # Estimate for the period using Kepler's 3rd Law
        # Used to give a first estimate for the curve fitting method
        return 2*np.pi*((((p.apo+p.per)/2)**3/(G*(p.mass+self.planet_list[p.primary].mass)))**(1/2))/(60*60*24)

    def curve_period(self,p,pos):
        # Fits the general function pfunc to the data stored in pos
        # Estimate_period is used to give the optimizer a starting point
        try:
            xdata = np.array([k for k in range(len(pos))])
            ydata = pos
            pe = self.estimated_period(p)*(60*60*24)/self.dt
            def pfunc(x,a,b,c): return a * np.absolute(np.sin(x*np.pi/b))**c

            params, _ = curve_fit(pfunc, xdata, ydata, [2,pe,1], bounds = (0.1,np.inf))

            return (params[1])*self.dt/(60*60*24)
        except RuntimeError:
            # Catch the error thrown by scipy if the curve fitting fails
            return 0


    def energy(self):
        # Sums the kinetic energy of every planet
        # with the potential energy of every pair.
        E = 0
        for i,p_i in enumerate(self.planet_list):
            E += p_i.kinetic_energy()
            for j,p_j in enumerate(self.planet_list):
                if j > i:
                    E += self.potential(p_i,p_j)
        return E

    def run(self, output):
        numstep = round(self.sim_time/self.dt)
        N = len(self.planet_list)

        self.CoM_fix()
        # Find initial forces
        self.update_forces()

        # Set up vmd_file and output string
        vmd_file = open(output, 'w')
        vmd_header = str(N) + '\n' + "Time = "

        E0 = self.energy()
        E=[]

        # Save initial positions
        positions = np.zeros((N-1,3))
        for j,p in enumerate(self.planet_list[1:]):
            positions[j] = p.position-self.planet_list[p.primary].position
            # Initial value for apsides
            r = np.linalg.norm(positions[j])
            p.apo = p.per = r
        pos =[[] for _ in self.planet_list[1:]]

        # Main time loop:
        for i in range(numstep):

            # Write position information
            if (i%10)==0:
                vmd_file.write(vmd_header + str(i*self.dt) + '\n')
                for p in self.planet_list:
                    vmd_file.write(str(p) + '\n')
                print(round(100*i/numstep), '%', end='\r')
                #Calculate percentage change of the total energy
            E.append(100*(self.energy()-E0)/E0)

            # Update particle position
            self.update_all('leap_pos2nd', self.forces)
            # Save last forces for average
            old_forces = np.copy(self.forces)
            # Calculate new forces
            self.update_forces()
            # Update particle velocity using average of old & new forces
            avg_force = (old_forces + self.forces)*0.5
            self.update_all('leap_velocity', avg_force)
            self.calc_apsides()

            for j,p in enumerate(self.planet_list[1:]):
                # Save (scaled) distance of planet from initial position in orbit to calculate period
                dist = np.linalg.norm(p.position-self.planet_list[p.primary].position-positions[j])
                pos[j].append(dist/np.linalg.norm(positions[j]))

            # Step forward time
            self.time += self.dt

        # After simulation:
        vmd_file.close()

        #Plot energy
        xdata = [self.dt*x for x in range(len(E))]
        plt.plot(xdata,E)
        plt.ylim(top=0.0005)
        plt.ylabel("Change in Energy(%)")
        plt.xlabel("Time(s)")
        plt.title("Percentage Change in Energy")
        plt.savefig("energy_fluctuations")

        # Calculate period
        for j,p in enumerate(self.planet_list[1:]):
            p.TCurve = self.curve_period(p,pos[j])
            p.TCount = self.count_period(pos[j])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Not enough arguments given.")
        quit()
    planet_f,param_f,traj_f = sys.argv[1:]

    t0 = time.time()
    S = PlanetSim(planet_f,param_f)
    S.run(traj_f)
    t1 = time.time()
    print('Simulation completed in {:f}s: dt = {}s, time = {}s'.format(t1-t0, S.dt, S.time))

    def print_table(table):
        #Neater console output of information
        longest_cols = [(max([len(str(row[i])) for row in table]) + 3) for i in range(len(table[0]))]
        row_format = "".join(["{:>" + str(longest_col) + "}" for longest_col in longest_cols])
        for row in table:
            print(row_format.format(*row))

    headers = ["Planet","Apoapsis(m)","Periapsis(m)","T-count(days)","T-curve(days)"]
    table = [[p.label,
                int(round(p.apo)),
                int(round(p.per)),
                    round(p.TCount,3),
                    round(p.TCurve,3)]
                for p in S.planet_list[1:]]

    table.insert(0,headers)
    print_table(table)
