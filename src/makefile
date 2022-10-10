FC=mpifort
FLAGS=-O3 -fopenmp -lfftw3
OBJECTS=  cdata.o derive.o energy.o output.o initial.o run.o 

%.o: %.f90
	$(FC) $(FLAGS) -c $<

cdata.o: cdata.f90 parameter.local
initial.o: initial.f90 cdata.o energy.o 
output.o: output.f90 cdata.o
energy.o: energy.f90 cdata.o
derive.o: derive.f90 cdata.o initial.o

default:run
run: $(OBJECTS) 
	$(FC) $(FLAGS) -o run $(OBJECTS) 
	rm -rf *.o *.mod
clean:
	rm -f *.o *.mod 

full-clean:clean
	rm -f *.txt 


