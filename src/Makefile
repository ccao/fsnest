include ../make.sys

objects = input.o flib.o  interpolate_bands.o calc_occ.o compute_fsnest.o output_fsnest.o

LIBS = ../modules/mod.a

$(LIBS):
	cd ../modules
	make mod.a
	cd ..

fsnest: fsnest.f90 $(objects) $(LIBS)
	$(F90) $(F90FLAGS) -fpp -o fsnest.x fsnest.f90 $(objects) $(LIBS) -I$(MKLINCLUDE) -I../modules -I. $(LAPACK95LIBS) $(MKLLIBS)

clean:
	rm -rf *.o *.mod *.x
