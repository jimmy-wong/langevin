### Source file structure

####src/cubic_spline.c

calculate the splined value in the gap of the grids.

####src/inertia.c

calculate the inertia tensor of the compound nuclei.

####src/dissipative.c

calculate the dissipative tensor of the compound nuclei.

####src/random.c

generate the random force on the nuclei.

####src/shape.c

class shape: calculate the surface normal speed and other shape
related terms.

####src/store.c

store the current PES.

####langevin.c

main program

### Program illustration

This program is used to simulate the dynamics of the fission process
on the potential energy surface (PES) in the framework of Langevin
equation. Currently, we only consider the temperature independent PES
because the excitation energy of the compound nuclei is very low. The
temperature of low excitation energy is at most 1MeV, which will not
bring too much correction to the PES. Also we perform the temperature
independent inertia and dissipative tensor on this surface.