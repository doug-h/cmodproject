"""
Particle3D, a class to describe 3D particles
"""

import numpy as np

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    mass(float) - particle mass
    label(string) - particle identifier
    position(numpy.array[float]) - position in xyz
    velocity(numpy.array[float]) - velocity in xyz

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates

    Static Methods:
    * create particle(s) from file
    * vector displacement between two particles
    """

    def __init__(self, label, pos, vel, mass):
        """
        Initialise a Particle3D instance

        :param label: label as string
        :param pos: position as numpy.array[float]
        :param vel: velocity as numpy.array[float]
        :param mass: mass as float
        """

        self.label = label
        self.position = pos
        self.velocity = vel
        self.mass = mass

    def __str__(self):
        """
        XYZ-compatible string representation
        Outputs as "<label> <x-pos> <y-pos> <z-pos>"
        """
        LXYZ = ([self.label]
                + [" " + str(x_i) for x_i in self.position])
        return "".join(LXYZ)

    def kinetic_energy(self):
        """
        :return: kinetic energy as float: 1/2*mass*|vel|^2
        """
        return 0.5*self.mass*(self.velocity@self.velocity)

    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m

        :param dt: timestep as float
        :param force: force on particle as numpy.array(float)
        """
        self.velocity += (dt/self.mass) * force

    def leap_pos1st(self, dt):
        """
        First-order position update,
        r(t+dt) = r(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position += dt*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as numpy.array(float)
        """
        self.position += dt*self.velocity + ((0.5/self.mass)*(dt**2))*force

    @staticmethod
    def from_file(fname):
        """
        Creates particle(s) from file
        Lines should start with '#' character then parameters in order:
        # <label> <x-pos> <y-pos> <z-pos> <x-vel> <y-vel> <z-vel> <mass>

        :param fname: filename as string

        :return: particle instance if file only contains one particle.
                 array of particle instances if file contatins multiple.
        """
        par_data = ['<label>', '<x-pos>', '<y-pos>', '<z-pos>', '<x-vel>', '<y-vel>', '<z-vel>', '<mass>']
        # only used to give more helpful error messages to user
        particles = []
        with open(fname, 'r') as f:
            for linenumber, line in enumerate(f):

                if line[0] != '#': #skip lines unless they start with '#'
                    continue 

                params = line.split()[1:] #cuts line into list and drops '# char

                if (len(params) != 8): #skip line if it has the wrong number of arguments
                    print("Warning - Line",linenumber + 1, "of", fname, 'was skipped:',
                          len(params), "arguments given when 8 were expected.")
                    print( 'Use form "#', ', '.join(par_data)+'".')
                    print()
                    continue 

                label = str(params[0])
                try:
                    for i in range(1, len(params)):
                        params[i] = float(params[i])
                except ValueError:
                    print("Warning - Line",linenumber + 1, "of", fname, 'was skipped:',
                          "unable to convert argument", i, "to float.")
                    print()
                    continue #skip line if it contains data the program doesn't understand

                else:
                    position = np.array(params[1:4])*1e3
                    velocity = np.array(params[4:7])*1e3
                    mass = float(params[7])
                    particles.append(Particle3D(label, position, velocity, mass))

        f.close() #not strictly necessary as 'with open()' closes automatically
        
        if len(particles) == 0:
            print("Warning - No valid particles were found in", fname)
            print()
        return particles



if __name__ == "__main__":
    #Create particles for testing functions
    """for i in range(10):
        p = Particle3D("p"+str(i),
                       np.array([2*i, 3*i, 5*i]),
                       np.array([7*i, 9*i, 11*i]),
                       13*i)
        print(p)
    """
    x = Particle3D.from_file("bodies.txt")
    for p in x:
        print(p)
