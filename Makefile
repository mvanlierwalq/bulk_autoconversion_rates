#Makefile for MCMC code

FC = gfortran
#FCFLAGS = -O0 -g -fbacktrace -fdefault-real-8 -fdefault-double-8
FCFLAGS = -O2 -march=native -fdefault-real-8 -fdefault-double-8
FLIBS = -lnetcdf -lnetcdff
#LPATH = -L/home/marcus/local/lib
#IPATH = -I/home/marcus/local/include
LPATH = -L/usr/local/lib
IPATH = -I/usr/local/include
OBJS = run_autoconversion_mod.o
MAIN = driver.f90

#..Main compilation
driver.x : $(MAIN) $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LPATH) $(FLIBS) $(IPATH)
#	$(FC) $(FCFLAGS) ForwardOperator.o fwd_model.o main.f90 -o mcmc.x $(LPATH) $(FLIBS) $(IPATH)

#..object files
#%.o : %.f90
#	$(FC) $(FCFLAGS) -c $< $(LPATH) $(FLIBS) $(IPATH)
#ForwardOperator.o : ForwardOperator.f90
#	$(FC) $(FCFLAGS) -c ForwardOperator.f90 $(LPATH) $(FLIBS) $(IPATH)
#

run_autoconversion_mod.o : run_autoconversion_mod.f90
	$(FC) $(FCFLAGS) -c run_autoconversion_mod.f90 $(LPATH) $(FLIBS) $(IPATH)

#..clean
clean :
	rm *.o *.mod
