FC = gfortran
F_OPT = -O2 -fbounds-check -fbacktrace -ffpe-trap=invalid
TARGET = a.out
OBJS = main.o analysis.o matrix.o element_C3D8.o solver.o io.o util.o

$(TARGET):$(OBJS)
	$(FC) $(F_OPT) -o $(TARGET) $(OBJS)

main.o: main.f90
	$(FC) -c $(F_OPT) main.f90
analysis.o: analysis.f90
	$(FC) -c $(F_OPT) analysis.f90
matrix.o: matrix.f90
	$(FC) -c $(F_OPT) matrix.f90
element_C3D8.o: element_C3D8.f90
	$(FC) -c $(F_OPT) element_C3D8.f90
solver.o: solver.f90
	$(FC) -c $(F_OPT) solver.f90
io.o: io.f90
	$(FC) -c $(F_OPT) io.f90
util.o: util.f90
	$(FC) -c $(F_OPT) util.f90

clean:
	rm *.o *.mod a.out

