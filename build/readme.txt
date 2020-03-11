
Usage
python planetsim.py planets.txt params.txt orbit.xyz
or
python3 planetsim.py planets.txt params.txt orbit.xyz

Orbital information is printed to the terminal.
Period is estimated most accurately using counting method, but
 for incomplete orbits a curve-fitting method gives an estimate.
 This is inaccurate for simulations much shorter than a full orbit.
Energy fluctuation (% change from initial value) is output as a .png file,
 can be hard to read for very long simulations with small time-steps.
