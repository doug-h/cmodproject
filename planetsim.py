from particle3D import Particle3D as P3D
import numpy as np
import sys

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Value of G taken from CODATA 2018
G=6.67430e-11

class Planet(P3D):

    template = {'label':(str,1), 'position':(float,3), 'velocity':(float,3), 'mass':(float,1), 'primary':(int,1)}

    def __init__(self, label, pos, vel, mass, pri):
        super().__init__(label, pos*1e3, vel*1e3, mass)
        self.primary = pri
        self.apo = self.per = -1
        self.period = 0
        self.period2 = 0

class PlanetSim(object):
    def __init__(self, planet_list):
        """
        Initialise a PlanetSim instance

        :param particle_list: list[Particle3D instances]
        """
        self.planet_list = planet_list
        self.forces = np.zeros((len(planet_list),3))


    def read_params(self,fname):
        """
        Read param file, extract parameters and save to sim. instance
        
        File in format: # <timestep> <duration> <initial time>
            
        :param fname: parameter file name (string)
        """
        with open(fname, 'r') as f:
            for line in f:
                #skip lines that don't start with '#'
                if line[0] != '#': 
                    continue 
                #cuts line into list and drops '# char
                params = line.split()[1:] 
                #convert params to floats
                params = [float(x) for x in params] 
        #save params to each variable
        self.dt, self.sim_time, self.time, = params

    def CoMfix(self):
        """
        Correct for centre-of-mass motion of planets
        
        CoM has velocity V = P/M, 
        where P is the total momentum of the planets, and M is the total mass.
        V is subtrated from all bodies.
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
        #Calulate pairwise force interactions
        n = len(self.planet_list)
        self.forces = np.zeros((n,3))
        for i,p_i in enumerate(self.planet_list):
            for j,p_j in enumerate(self.planet_list):
                if j > i:
                    f = self.force(p_i, p_j)
                    self.forces[i] += f
                    self.forces[j] -= f

    def update_all(self, method, forces):
        for i,p_i in enumerate(self.planet_list):
            getattr(p_i, method)(self.dt, forces[i])

    def calc_apsides(self):
        for p in self.planet_list:
            r = np.linalg.norm(p.position - self.planet_list[p.primary].position)
            if p.per < 0:
                p.per = r
            elif r < p.per:
                p.per = r
            if p.apo < 0:
                p.apo = r
            elif r > p.apo:
                p.apo = r

    def estimated_period(self,p):
        #estimate for the period using Kepler's 3rd Law
        return 2*np.pi*((((p.apo+p.per)/2)**3/(G*(p.mass+self.planet_list[p.primary].mass)))**(1/2))/(60*60*24)
    
    def calc_period(self,pos):
        plotsx = [0]
        plotsi = [0]
        P=0
        for i in range(len(pos)-2):
            if pos[i]-pos[i+1]>0:
                if pos[i+2]-pos[i+1]>0:
                    plotsx += [pos[i+1]]
                    plotsi += [i+1]
        if len(plotsi)>1:
            P=((plotsi[-1]-plotsi[0])/(len(plotsi)-1))
            
        return plotsi,plotsx,P

    def energy(self):
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

        self.CoMfix() #remove velocity of COM
        self.update_forces() #find initial forces

        self.calc_apsides()

        vmd_file = open(output, 'w')

        vmd_header = str(N) + '\n' + "Time = "

        E0 = self.energy()
        E=[]


        s = np.zeros((N-1,3))
        for j,p in enumerate(self.planet_list[1:]):
            s[j] = p.position-self.planet_list[p.primary].position
        pos =[[] for _ in self.planet_list[1:]]
        for i in range(numstep):
            for j,p in enumerate(self.planet_list[1:]):
                pos[j] += [np.linalg.norm(p.position-self.planet_list[p.primary].position-s[j])/np.linalg.norm(s[j])]

            
            vmd_file.write(vmd_header + str(i*self.dt) + '\n')
            for p in self.planet_list:
                vmd_file.write(str(p) + '\n')
            #print(round(100*i/numstep), '%')

            # save last forces for average
            old_forces = np.copy(self.forces)

            # Update particle position
            self.update_all('leap_pos2nd', self.forces)

            # Calculate new forces
            self.update_forces()

            # Update particle velocity using average of old & new forces
            avg_force = (old_forces + self.forces)*0.5
            self.update_all('leap_velocity', avg_force)

            # Step forward time
            self.time += self.dt

            self.calc_apsides()

            E+=[self.energy()-E0]
        xdata = np.array([k for k in range(numstep)])

        e_fig = plt.figure()
        plt.scatter(xdata,E)
        plt.title("E(t)-E(t0)")


        vmd_file.close()
        
        
        fig = plt.figure(figsize=(12, 10))
        for j,p in enumerate(self.planet_list[1:]):
            xdata = np.array([k for k in range(len(pos[j]))])
            ydata = pos[j]
            n=round(np.sqrt(N))
            while n**2<len(self.planet_list[1:]):
                n+=1
            fig.add_subplot(n,n,j+1)
            plt.scatter(xdata,ydata,s=1)
            pers = self.calc_period(pos[j])
            plt.scatter(*pers[:-1])
            plt.title(p.label)
            pe = self.estimated_period(p)*(60*60*24)/self.dt
            def pfunc(x,a,b,c): return a * np.absolute(np.sin(x*np.pi/b))**c
            popt, pcov = curve_fit(pfunc, xdata,ydata, [2,pe,1], bounds = (0.1,np.inf))
            p.period = pers[-1]
            p.period2=(popt[1])*self.dt/(60*60*24)
            plt.plot(xdata, pfunc(xdata, *popt), 'g',
                label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
            plt.legend()
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("not enough args")
        quit()
    planet_f,param_f,traj_f = sys.argv[1:]
        
    planets = Planet.from_file(planet_f)
    S = PlanetSim(planets)
    S.read_params(param_f)
    S.run(traj_f)
    for p in S.planet_list:
        print(p.label,p.apo,p.per,S.estimated_period(p),p.period,p.period2)
