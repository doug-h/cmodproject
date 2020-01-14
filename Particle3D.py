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
        params = []
        with open(fname, 'r') as f:
            # only look at lines beginning with '#', and remove newline char
            params = [line.rstrip() for line in f if line[0] == '#']
            # seperate parameters and remove '#' character
            params = [p.split()[1:] for p in params]

        f.close() #not strictly necessary as 'with open()' closes automatically

        if len(params) == 0:
            print("No particles found.")
            print("Particles should be in form: '#' <x> <y> <z> ...")
            return None
        else:
            particles = []
            for p in params:
                label = p[0]
                position = np.array([float(r_i) for r_i in p[1:4]])
                velocity = np.array([float(v_i) for v_i in p[4:7]])
                mass = float(p[7])
                particles.append(Particle3D(label, position, velocity, mass))
            if len(particles) == 1:
                return particles[0]
            else:
                return particles

    @staticmethod
    def seperation(p1, p2):
        """
        Finds seperation between particles;
        vector p2.pos - p1.pos

        :param p1,p2: particles as Particle3D instances
        """
        return p2.position - p1.position


if __name__ == "__main__":
    #Create particles for testing functions
    for i in range(10):
        p = Particle3D("p"+str(i),
                       np.array([2*i, 3*i, 5*i]),
                       np.array([7*i, 9*i, 11*i]),
                       13*i)
        print(p)

    x = Particle3D.from_file("param.txt")
    print(x[0],x[1])
