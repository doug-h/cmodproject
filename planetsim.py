from particle3D import Particle3D as P3D
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

G=6.67430e-11

class Planet(P3D):

    template = {'label':(str,1), 'position':(float,3), 'velocity':(float,3), 'mass':(float,1), 'primary':(int,1)}

    def __init__(self, label, pos, vel, mass, pri):
        super().__init__(label, pos*1e3, vel*1e3, mass)
        self.primary = pri
        self.apo = self.per = -1
        self.period = 0

class PlanetSim(object):
    def __init__(self, planet_list):
        """
        Initialise a PlanetSim instance

        :param particle_list: list[Particle3D instances]
        :param dt: timestep as float
        :param sim_time: simulation length as float
        :param t0: initial time as float
        """
        self.planet_list = planet_list
        self.forces = np.zeros((len(planet_list),3))
        self.dt = 60*60
        self.time = 0
        self.sim_time = 60*60*24*365


    def read_params(self,fname):
        with open(fname, 'r') as f:
            for linenumber, line in enumerate(f):
                if line[0] != '#': #skip lines unless they start with '#'
                    continue
                params = line.split()[1:] #cuts line into list and drops '# char
                params = [float(x) for x in params]

        self.dt, self.sim_time, self.time, = params

    def COMfix(self):
        P = np.zeros(3)
        M = 0
        for p in self.planet_list:
            P += (p.mass * p.velocity)
            M += p.mass
        v_com = P/M
        for p in self.planet_list:
            p.velocity -= v_com

    def force(self, p0, p1):
        r_01 = p1.position - p0.position
        m_01 = np.linalg.norm(r_01)
        return r_01*G*p0.mass*p1.mass*(m_01**-3)

    def potential(self, p0, p1):
        r_01 = p1.position - p0.position
        m_01 = np.linalg.norm(r_01)
        return -G*p0.mass*p1.mass/m_01

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

        self.COMfix() #remove velocity of COM
        self.update_forces() #find initial forces

        self.calc_apsides()

        vmd_file = open('out/' + output, 'w')

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

            if i%24==0:
                vmd_file.write(vmd_header + str(i*self.dt) + '\n')
                for p in self.planet_list:
                    vmd_file.write(str(p) + '\n')
                print(round(100*i/numstep), '%')

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
            fig.add_subplot(n,n,j+1)
            plt.scatter(xdata,ydata,s=1)
            plt.title(p.label)
            pe = self.estimated_period(p)*(60*60*24)/self.dt
            def pfunc(x,a,b,c): return a * np.absolute(np.sin(x*np.pi/b))**c
            popt, pcov = curve_fit(pfunc, xdata,ydata, [2,pe,1])
            p.period=(popt[1])*self.dt/(60*60*24)
            plt.plot(xdata, pfunc(xdata, *popt), 'g',
                label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
            plt.legend()
        plt.show()

if __name__ == "__main__":
    planets = Planet.from_file("all.txt")
    S = PlanetSim(planets)
    S.read_params("params.txt")
    S.run("orbit.xyz")
    for p in S.planet_list:
        print(p.label,p.apo,p.per,S.estimated_period(p),p.period)
