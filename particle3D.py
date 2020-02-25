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
    * second order position update

    Static Methods:
    * create particle(s) from file
    """
    template = {'label':(str,1), 'position':(float,3), 'velocity':(float,3), 'mass':(float,1)}
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
                + [str(x_i) for x_i in self.position])
        return " ".join(LXYZ)

    def kinetic_energy(self):
        """
        :return: kinetic energy as float: 1/2*mass*|vel|^2
        """
        return 0.5*self.mass*((np.linalg.norm(self.velocity))**2)

    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m

        :param dt: timestep as float
        :param force: force on particle as numpy.array(float)
        """
        self.velocity += (dt/self.mass) * force

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)/m

        :param dt: timestep as float
        :param force: force as numpy.array(float)
        """
        self.position += dt*self.velocity + ((0.5/self.mass)*(dt**2))*force

    @classmethod
    def from_file(cls,fname):
        particles = []
        with open(fname, 'r') as f:
            for linenumber, line in enumerate(f):
                if line[0] != '#': #skip lines unless they start with '#'
                    continue
                params = line.split()[1:] #cuts line into list and drops '# char
                args = []
                for key in cls.template:
                    d_type, d_len = cls.template[key]
                    try:
                        if d_len == 1:
                            args.append(d_type(params.pop(0)))
                        else:
                            args.append(np.array([d_type(params.pop(0)) for _ in range(d_len)]))
                    except ValueError:
                        print("Warning - Line",linenumber + 1, "of", fname, ':',
                          "unable to convert argument", key, "to", d_type, ".")
                    except IndexError:
                        print("Warning - Line",linenumber + 1, "of", fname, ':',
                          "not enough arguments given.")
                try:
                    p = cls(*args)
                    particles.append(p)
                except TypeError:
                    print("Error:", '"' + args[0] + '"', "failed.")
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
    x = Particle3D.from_file("all.txt")
    for p in x:
        print(p)
