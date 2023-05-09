# cmodproject
Code I co-wrote for the Computer Modelling (2019) course.
Simulates the solar system to predict orbital information.

Usage:
From within the build directory, type
```
python planetsim.py planets.txt params.txt orbit.xyz
```


Orbital information is printed to the terminal.
Position data for the planets is saved to orbit.xyz

Period is estimated most accurately using counting method, but
 for incomplete orbits a curve-fitting method gives an estimate.
 This is inaccurate for simulations much shorter than a full orbit.
