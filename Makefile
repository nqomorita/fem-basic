FC       = gfortran
FFLAGS   = -O2 -fbounds-check -fbacktrace -ffpe-trap=invalid
LDFLAGS  =
LIBS     =
INCLUDE  = -I ./include
MOD_DIR  = -J ./include
BIN_DIR  = ./bin
SRC_DIR  = ./src
OBJ_DIR  = ./build
BIN_LIST = a.out
TARGET   = $(addprefix $(BIN_DIR)/, $(BIN_LIST))
SRC_LIST = main.f90 analysis.f90 update.f90 matrix.f90 shape_C3D8.f90 element_C3D8.f90 solver.f90 io.f90 util.f90
SOURCES  = $(addprefix $(SRC_DIR)/, $(SRC_LIST))
OBJS     = $(subst $(SRC_DIR), $(OBJ_DIR), $(SOURCES:.f90=.o))
RM       = rm

$(TARGET): $(OBJS) $(LIBS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) $(MOD_DIR) -o $@ -c $<

clean:
	$(RM) $(OBJS) $(TARGET) ./include/*.mod

.PHONY: clean
