"""
Particle3D, a class to describe 3D particles
Douglas Hull & Jack Manners
2/3/20
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
    * formatted output string
    * kinetic energy
    * first-order velocity update
    * second order position update

    Class Method:
    * create particle(s) from file
    """

    # Data template for reading from file
    template = {'label': (str, 1),
                'position': (float, 3),
                'velocity': (float, 3),
                'mass': (float, 1)}

    def __init__(self, label, pos, vel, mass):
        """
        Initialise a Particle3D instance

        :param label: label as string
        :param pos: position as numpy.array[float]
        :param vel: velocity as numpy.array[float]
        :param mass: mass as float
        """
        self.label = label
        self.position = pos * 1e3
        self.velocity = vel * 1e3
        self.mass = mass

    def __str__(self):
        """
        VMD-compatible string representation
        Outputs as "<label> <x-pos> <y-pos> <z-pos>"
        """
        LXYZ = ([self.label] + [str(x_i) for x_i in self.position])
        return " ".join(LXYZ)

    def kinetic_energy(self):
        """
        Kinetic energy as float: 1/2*mass*|vel|^2
        """
        return 0.5 * self.mass * ((np.linalg.norm(self.velocity))**2)

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m

        :param dt: timestep as float
        :param force: force on particle as numpy.array[float]
        """
        self.velocity += (dt / self.mass) * force

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)/m

        :param dt: timestep as float
        :param force: force as numpy.array(float)
        """
        self.position += dt * self.velocity + \
            ((0.5 / self.mass) * (dt**2)) * force

    @classmethod
    def from_file(cls, fname):
        """
        Creates a list of particles using data from a file.
        The data's structure is given by the class variable 'template':
        {'label':str,
        'position':[float,3],
        'velocity':[float,3],
        'mass':[float,1]}

        :param fname: name of file containing particle data
        """

        particles = []
        with open(fname, 'r') as f:
            for linenumber, line in enumerate(f):

                if line[0] != '#':
                    # Skips lines unless they start with '#'
                    continue

                # Cuts line into list and drops '#' char
                params = line.split()[1:]
                args = []

                # Runs through each parameter in the template:
                for key in cls.template:
                    # Gets data type and length of the parameter
                    d_type, d_len = cls.template[key]
                    try:
                        if d_len == 1:
                            # Converts a single read-in value to the type
                            # given in 'template'
                            arg = d_type(params.pop(0))
                            args.append(arg)
                        else:
                            # Converts multiple read-in values to the type
                            # given in 'template' and stores in an np.array
                            arg = np.array([d_type(params.pop(0))
                                            for _ in range(d_len)])
                            args.append(arg)

                    # Catching invalid input
                    except ValueError:
                        print("Warning - Line", linenumber + 1, "of", fname, ':',
                              "unable to convert argument", key, "to", d_type, ".")
                    except IndexError:
                        print("Warning - Line", linenumber + 1, "of", fname, ':',
                              "not enough arguments given.")
                # try:
                    # Creates a particle using args
                p = cls(*args)
                particles.append(p)
                # except TypeError:
                #    print("Error:", '"' + args[0] + '"', "failed.")
        f.close()
        return particles


if __name__ == "__main__":
    # Create particles for testing
    x = Particle3D.from_file("all.txt")
    for p in x:
        print(p)
