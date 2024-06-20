ifort -c -O2 ~/PepBD/Helpers/datatypes.f90
ifort -c -O2 ~/PepBD/Helpers/sys_vars.f90
ifort -c -O2 ~/PepBD/Main/input.f90
ifort -c -O2 ~/PepBD/Helpers/math.f90
ifort -c -O2 ~/PepBD/Helpers/utilities.f90
ifort -c -O2 ~/PepBD/Helpers/pdbfile.f90
ifort -c -O2 ~/PepBD/Helpers/database.f90 -diag-disable 8291
ifort -c -O2 ~/PepBD/Helpers/surface_area.f90
ifort -c -O2 ~/PepBD/Helpers/energy_calculation.f90
ifort -c -O2 ~/PepBD/Helpers/advanced_function.f90
ifort -c -O2 ~/PepBD/Helpers/optimization_techniques.f90
ifort -c -O2 ~/PepBD/Main/PepBD.f90
ifort -O2 -o main *.o -mkl=sequential