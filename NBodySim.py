from Particle3D import Particle3D as P3D
import numpy as np

class NBodySim(object):
    def __init__(self, particle_file, param_file):
          self.particle_list = [P3D("p"+str(1), np.array([0, 0, 0]), np.array([1, 1, 1]), 3),
                                P3D("p"+str(1), np.array([0, 0, 0]), np.array([1, 2, 1]), 1)]
    """
     def COMfix(self):
          update particle_list
          
     def generic force:
     
     def generic energy:
    """
    def kinetic_energy(self):
        K = 0
        for p in self.particle_list:
            K += p.kinetic_energy()
        return K

     
S = NBodySim("particle.txt", "param.txt")
print(S.kinetic_energy())