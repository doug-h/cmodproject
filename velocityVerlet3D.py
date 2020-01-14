"""
Velocity verlet time integration of
two particles interacting via Morse potential.

Outputs files containing particle seperation and
total energy of system over time, as well as plots of both,
are generated in /out/xxx.dat, /out/xxx.png.

"""


import sys
import numpy as np
from Particle3D import Particle3D as P3D
import os
import matplotlib.pyplot as plt


def pot_energy_m(p1, p2, a, D_e, r_e):
    """
    Method to return potential energy of pair particles with Morse potential.

    :params p1, p2: Particle3D instances
    :param a: curvature of potential minimum
    :param D_e: depth of potential minimum
    :param r_e: position of potential minimum

    :return: potential energy of pair as float
    """
    r_12_m = np.linalg.norm(P3D.seperation(p1, p2))
    potential = D_e * ((1 - np.exp(-a * (r_12_m-r_e)))**2 - 1)
    return potential


def force_m(p1, p2, a, D_e, r_e):
    """
    Method to return the force between pair particles with Morse potential.
    Due to symmetry, the force on p2 is equal and opposite to the force on p1.

    :params p1, p2: Particle3D instances
    :param a: curvature of potential minimum
    :param D_e: depth of potential minimum
    :param r_e: position of potential minimum

    :return: force acting on p1 as Numpy.array[float]
    """
    r_12 = P3D.seperation(p1, p2)
    r_12_m = np.linalg.norm(r_12)

    z = np.exp(-a * (r_12_m-r_e))
    force = 2 * a * D_e * (1-z) * z
    return (force/r_12_m) * r_12


# Begin main code
def main():
    # Read name of output files from command line
    if len(sys.argv) != 2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <parameter file>")
        quit()
    else:
        param_file = open(sys.argv[1], 'r')

    # Read sim paramters from file
    # Take lines beginning with '@' and remove newline char
    params = [line.rstrip() for line in param_file if line[0] == '@']
    # Seperate parameters, convert to floats and remove '@' char
    params = [[float(v) for v in p.split(" ")[1:]] for p in params]
    param_file.close()

    p1, p2, *_ = P3D.from_file(sys.argv[1])
    dt, sim_time, t0, *_ = params[0]
    r_e, a, D_e, *_ = params[1]

    # Set up simulation
    time = t0
    numstep = round(sim_time/dt)  # so simulation-time independent of step size

    # find number of decimal places of dt so time is outputted correctly
    time_prec = len(str(dt).split('.')[-1])

    # Physical constants taken from CODATA
    eV = 1.602176634e-19  # (J)
    amu = 1.66053906660e-27  # (kg)
    ang = 1e-10  # (m)

    # Specific force/potential functions with user-entered properties
    def force_ms(p1, p2):
        return force_m(p1, p2, a, D_e, r_e)

    def pot_energy_ms(p1, p2):
        return pot_energy_m(p1, p2, a, D_e, r_e)

    # Write out initial conditions
    energy = (p1.kinetic_energy()
              + p2.kinetic_energy()
              + pot_energy_ms(p1, p2))

    # Open output files
    sep_file = open("out/seperation.dat", "w")
    energy_file = open("out/energy.dat", "w")

    # write initial seperation/energy
    r_12 = np.linalg.norm(P3D.seperation(p1, p2))
    sep_file.write("{:#.{tp}f} {:.8f}\n".format(time, r_12, tp=time_prec))
    energy_file.write("{:#.{tp}f} {:.14f}\n".format(time, energy, tp=time_prec))

    # Get initial force
    force = force_ms(p1, p2)

    # Initialise data lists for plotting later
    time_list = [time]
    sep_list = [r_12]
    energy_list = [energy]


    # Start the time integration loop
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
        time += dt

        # Output particle information
        energy = (p1.kinetic_energy()
                  + p2.kinetic_energy()
                  + pot_energy_ms(p1, p2))

        r_12 = np.linalg.norm(P3D.seperation(p1, p2))

        sep_file.write("{:#.{tp}f} {:.8f}\n".format(time, r_12, tp=time_prec))

        energy_file.write("{:#.{tp}f} {:.14f}\n".format(time, energy, tp=time_prec))

        # Append information to data lists
        time_list.append(time)
        sep_list.append(r_12)
        energy_list.append(energy)

    # Post-simulation:
    # Close output file
    sep_file.close()
    energy_file.close()

    # Convert units for plotting
    u = ang*(amu/eV)**(0.5)
    time_list = [t*u*1e15 for t in time_list]

    # Save plot of particles' seperation
    fig_s, ax_s = plt.subplots()
    ax_s.set_title('Velocity Verlet: Particle Seperation vs Time')
    ax_s.set_xlabel('Time (fs)')
    ax_s.set_ylabel('Seperation (A)')
    ax_s.plot(time_list, sep_list, label="r_12")
    # add line showing equilibrium bond distance
    ax_s.axhline(y=r_e, color='g', label="r_e")
    ax_s.grid(which='both')
    ax_s.legend(loc=4)
    #ax_s.set_xlim(0, time_list[-1])
    fig_s.savefig("out/seperation.png")


    # Save plot of total energy
    fig_e, ax_e = plt.subplots()
    ax_e.set_title('Velocity Verlet: Total Energy vs Time')
    ax_e.set_xlabel('Time (fs)')
    ax_e.set_ylabel("Energy (eV)")
    ax_e.set_xlim(0, time_list[-1])
    ax_e.plot(time_list, energy_list, label='Total Energy')
    #ax_e.legend(loc=4)
    ax_e.ticklabel_format(useOffset=False)
    ax_e.grid(which='both')
    fig_e.tight_layout()
    fig_e.savefig("out/energy.png")
    #plt.show()

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
