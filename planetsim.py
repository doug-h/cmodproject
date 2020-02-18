from particle3D import Particle3D as P3D
import numpy as np

G=6.67430e-11

class Planet(P3D):
    template = {'label':(str,1), 'position':(float,3), 'velocity':(float,3), 'mass':(float,1), 'primary':(int,1)}
    def __init__(self, label, pos, vel, mass, pri):
        super().__init__(label, pos, vel, mass)
        self.primary = pri
        self.per = 1e100
        self.apo = 0
        self.esti_period = 0
        self.period = 0
        self.init_pos = np.copy(pos)*1e3

class PlanetSim(object):
    def __init__(self, planet_list, dt=60*60, sim_time=60*60*24*365*6, t0 = 0):
        """
        Initialise a PlanetSim instance

        :param particle_list: list[Particle3D instances]
        :param dt: timestep as float
        :param sim_time: simulation length as float
        :param t0: initial time as float
        """
        self.planet_list = planet_list
        self.forces = np.zeros((len(planet_list),3))
        self.dt = dt
        self.t0 = t0
        self.time = t0
        self.sim_time = sim_time

        self.tp = int(np.log10(sim_time) // 1) + 1
        self.pos_file = None

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
                if r < p.per:
                    p.per = r
                if r > p.apo:
                    p.apo = r

    def esti_period(self):
        for p in self.planet_list:
            p.esti_period = 2*np.pi*((p.apo**3/(G*self.planet_list[p.primary].mass))**(1/2))/(60*60*24)

    def calc_period(self,pos):
        for i,p in enumerate(self.planet_list):
            j0=0
            for j in range(len(pos[i])-2):
                if pos[i][j+1] < pos[i][j] and pos[i][j+2] > pos[i][j+1]:
                    if j0:
                        T = (j+1-j0)/24
                        print(T)
                        p.period += T
                        if p.period != T:
                            p.period /= 2
                    j0 = j+1

    def kinetic_energy(self):
        K = 0
        for p in self.planet_list:
            K += p.kinetic_energy()
        return K

    def write_output(self, file):
        for p in self.planet_list:
            file.write('{:0{tp}d} {} \n'.format(self.time, p, tp = self.tp))

    def run(self, output):
        self.time = self.t0
        numstep = round(self.sim_time/self.dt)

        self.COMfix() #remove velocity of COM

        self.update_forces() #find initial forces

        self.calc_apsides()

        self.esti_period()

        pos_file = open('out/' + output, 'w')

        self.write_output(pos_file)

        test = [[] for _ in range(len(self.planet_list))]

        for i in range(numstep):
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

            for j,p in enumerate(self.planet_list):
                test[j] += [np.linalg.norm(p.position-p.init_pos-self.planet_list[p.primary].position)]
        self.calc_period(test)


if __name__ == "__main__":
    S = PlanetSim(Planet.from_file("all.txt"))
    S.run("orbit.dat")

    for p in S.planet_list:
        print(p.label,p.apo,p.per,p.esti_period,p.period)
