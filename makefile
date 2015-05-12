FCOMP    = ifort
OPTS     =  -c -traceback -CB
LINKOPTS =  -o
OBJS =  Modules.o NS_Solver.o PoissonSolver.o GCM.o

NS_Solver:$(OBJS)
	$(FCOMP) $(LINKOPTS) NS.out $(OBJS)

clean:
	rm  *.o *.mod NS.out

.SUFFIXES: .o .F .c .f .swp .f90

Modules.o: Modules.F90
	$(FCOMP) $(OPTS) Modules.F90

NS_Solver.o: NS_Solver.F90
	$(FCOMP) $(OPTS) NS_Solver.F90
	
PoissonSolver.o: PoissonSolver.F90
	$(FCOMP) $(OPTS) PoissonSolver.F90
	
GCM.o: GCM.F90
	$(FCOMP) $(OPTS) GCM.F90
