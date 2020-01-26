from particle3D import Particle3D as P3D
import numpy as np

G=6.67430e-11

class NBodySim(object):
    def __init__(self, particle_list, dt=10, sim_time=315576, t0 = 0):
          self.particle_list = particle_list
          self.forces = np.zeros((len(particle_list),3))
          self.dt = dt
          self.t0 = t0
          self.time = t0
          self.sim_time = sim_time
          self.tp = int(np.log10(sim_time) // 1) + 1
          self.pos_file = None
    """
     def COMfix(self):
          update particle_list

        """
    def force(self, p0, p1):
        r_01 = p1.position - p0.position
        m_01 = np.linalg.norm(r_01)
        return r_01*G*p0.mass*p1.mass*(m_01**-3)

    def update_forces(self):
        #Calulate pairwise force interactions
        n = len(self.particle_list)
        self.forces = np.zeros((n,3))
        for i,p_i in enumerate(self.particle_list):
            for j,p_j in enumerate(self.particle_list):
                if j > i:
                    f = self.force(p_i, p_j)
                    self.forces[i] += f
                    self.forces[j] -= f

    def update_all(self, method, forces):
        for i,p_i in enumerate(self.particle_list):
            getattr(p_i, method)(self.dt, forces[i])

    def kinetic_energy(self):
        K = 0
        for p in self.particle_list:
            K += p.kinetic_energy()
        return K


    def write_output(self, file):
        for p in self.particle_list:
            file.write('{:0{tp}d} {} \n'.format(self.time, p, tp = self.tp))


    def run(self, output):
        self.time = self.t0
        numstep = round(self.sim_time/self.dt)

        self.update_forces() #find initial forces

        pos_file = open('out/' + output, 'w')

        self.write_output(pos_file)

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

            if (i%3600==0):
                self.write_output(pos_file)

        pos_file.close()


if __name__ == "__main__":
    S = NBodySim(P3D.from_file("bodies.txt"))
    S.update_forces()
    S.run("orbit.dat")
    for p in S.particle_list:
        print(p)
