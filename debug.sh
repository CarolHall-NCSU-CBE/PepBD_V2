ifort -c -g -check all -fpe0 -warn -traceback -debug extended datatypes.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended sys_vars.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended input.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended math.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended utilities.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended pdbfile.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended database.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended surface_area.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended energy_calculation.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended advanced_function.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended optimization_techniques.f90
ifort -c -g -check all -fpe0 -warn -traceback -debug extended PepBD.f90
ifort -g -check all -fpe0 -warn -traceback -debug extended -o main *.o -mkl=sequential
