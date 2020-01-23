from Particle3D import Particle3D as P3D
import numpy as np

G=1

class NBodySim(object):
    def __init__(self, particle_list, dt=0.01, sim_time=10):
          self.particle_list = particle_list
          self.forces = np.zeros((len(particle_list),3))
    """
     def COMfix(self):
          update particle_list
        
        """  
    def force(self, p0, p1):
        r_01 = p1.position - p0.position
        m_01 = r_01@r_01
        return G*p0.mass*p1.mass*m_01**(-3/2)*r_01

    def update_forces(self):
        n = len(self.particle_list)
        self.forces = np.zeros((n,3))
        for i,pi in enumerate(self.particle_list):
            for j,pj in enumerate(self.particle_list):
                if j > i:
                    f = self.force(pi, pj)
                    self.forces[i] += f
                    self.forces[j] -= f


    def kinetic_energy(self):
        K = 0
        for p in self.particle_list:
            K += p.kinetic_energy()
        return K

    """def run(self, output):
        for i in range(numstep):
        # Update particle position
        p1.leap_pos2nd(dt, force)
        p2.leap_pos2nd(dt, -force)

        # Calculate force
        force_new = force_ms(p1, p2)

        # Update particle velocity using average of old & new forces
        avg_force = (force+force_new)*0.5
        p1.leap_velocity(dt, avg_force)
        p2.leap_velocity(dt,-avg_force)

        # Step forward time
        force = force_new
        time += dt"""

     
S = NBodySim(P3D.from_file("bodies.txt"), "param.txt")
print(S.kinetic_energy())
print(S.forces)
S.update_forces()
print(S.forces)