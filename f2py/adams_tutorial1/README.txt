Before starting the tutorial, clean the directory by

> rm *.exe *.mod *.so

The file "trig_subroutines_module.f90" contains a 
Fortran _module_ called "trig_subroutines_md".

That module contains a couple of _subroutines_ called 
"compute_area" and "compute_total_area". 

Now, those subroutines can be called either by a Fortran
program, or imported as a module into Python by f2py.

The file "trig_program.f90" is a Fortran program which calls
both subroutines. It can be compiled and exected by (e.g.)

> gfortran -o test.exe trig_subroutines_module.f90 trig_program.f90
> ./test.exe 

(note that after the first command, a .mod file will have been 
created)

To use f2py, by contrast, we first have to comment out the module 
part of "trig_subroutines_module.f90", leaving a file with only 
the submodules. This is done in "trig_subroutines_for_python.f90".
We can then use those subroutines in Python by (e.g.)

> f2py -c trig_subroutines_for_python.f90 -m trig
> python 
>>> import trig
>>> trig.compute_area(10)
>>> trig.compute_total_area([1,2,3])
>>> exit()

The tutorial is finished :) 
