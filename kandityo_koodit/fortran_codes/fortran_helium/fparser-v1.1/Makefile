F90          = g95
F90FLAGS     = -g

DEST         = testfp
BASIS        = parameters.o fparser.o testfp.o
MODULES      = parameters.mod fparser.mod

.SUFFIXES:
.SUFFIXES: .f90 .mod .o

.f90.o: 
	$(F90) $(F90FLAGS) -c $<

.f90.mod:
	$(F90) $(F90FLAGS) -c $<

testfp: $(BASIS) $(MODULES)
	$(F90) -o $(DEST) $(F90FLAGS) $(BASIS)

clean:
	rm -fr *.o *.mod testfp

testfp.o: parameters.o fparser.o

fparser.o: parameters.o






