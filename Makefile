# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -g -I. -IClasses1D -Imain

# Output directory
BIN_DIR := bin

# Shared source files
CLASSES1D_SRCS := $(wildcard Classes1D/*.cpp)
GEOMETRY_SRCS := $(wildcard Classes1D/Geometry/*.cpp)
THERMOPHYSICAL_PROPERTIES_SRCS := $(wildcard Classes1D/Thermophysical_Properties/*.cpp)
SOLVERSETTINGS_SRCS := $(wildcard Classes1D/SolverSettings/*.cpp)
MESH1D_SRCS := $(wildcard Classes1D/Mesh1D/*.cpp)
FIELD1D_SRCS := $(wildcard Classes1D/Field1D/*.cpp)
INITIALCONDITIONS_SRCS := $(wildcard Classes1D/InitialConditions/*.cpp)
BOUNDARIES_SRCS := $(wildcard Classes1D/Boundaries/*.cpp)
BOUNDARYCONDITIONS_SRCS := $(wildcard Classes1D/BoundaryConditions/*.cpp)
FILEWRITER1D_SRCS := $(wildcard Classes1D/FileWriter1D/*.cpp)
CONVERGENCE_SRCS := $(wildcard Classes1D/Convergence/*.cpp)
ANALYTICALSOLUTION1D_SRCS := $(wildcard Classes1D/AnalyticalSolution1D/*.cpp)
PLOTTER_SRCS := $(wildcard Classes1D/Plotter/*.cpp)
PRINTER1D_SRCS := $(wildcard Classes1D/Printer1D/*.cpp)
IMPLICIT_SOLVERS1D_SRCS := $(wildcard Classes1D/ImplicitSolvers1D/*.cpp) 

# Headers used as prerequisites for rebuilds
ESSENTIALHEADERS_SRCS := $(wildcard EssentialHeaders/*.hpp)

# Specific implicit solver source files
COEFFICIENTMATRIX_SRCS := Classes1D/Matrices/CoefficientMatrix.cpp
TRIDIAGONALMATRIX_SRCS := Classes1D/Matrices/TridiagonalMatrix.cpp
RHS_SRCS := Classes1D/Matrices/RHS.cpp

GAUSSIANELIMINATION_SRCS := Classes1D/LinearSolvers/GaussianElimination.cpp
THOMASALGORITHM_SRCS := Classes1D/LinearSolvers/ThomasAlgorithm.cpp


# Common shared sources for the single program
COMMON_SRCS := $(CLASSES1D_SRCS) \
               $(GEOMETRY_SRCS) \
               $(THERMOPHYSICAL_PROPERTIES_SRCS) \
               $(SOLVERSETTINGS_SRCS) \
               $(MESH1D_SRCS) \
               $(FIELD1D_SRCS) \
               $(INITIALCONDITIONS_SRCS) \
               $(BOUNDARIES_SRCS) \
               $(BOUNDARYCONDITIONS_SRCS) \
               $(PRINTER1D_SRCS) \
               $(FILEWRITER1D_SRCS) \
               $(ANALYTICALSOLUTION1D_SRCS) \
               $(CONVERGENCE_SRCS) \
               $(PLOTTER_SRCS) \
               $(RHS_SRCS) \
               $(IMPLICIT_SOLVERS1D_SRCS)


# Program registry
PROGRAMS := MAIN

# Main source file
MAIN_SRC := main/main.cpp

# Additional Optional source files required by both methods
MAIN_EXTRA_SRCS := \
    $(COEFFICIENTMATRIX_SRCS) \
    $(TRIDIAGONALMATRIX_SRCS) \
    $(GAUSSIANELIMINATION_SRCS) \
    $(THOMASALGORITHM_SRCS) \
    

MAIN_EXE := $(BIN_DIR)/main

# Aggregate lists
ALL_EXES := $(foreach prog,$(PROGRAMS),$($(prog)_EXE))

.PHONY: all clean main

all: $(ALL_EXES)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

define build_program
$($(1)_EXE): $($(1)_SRC) $(COMMON_SRCS) $($(1)_EXTRA_SRCS) $(ESSENTIALHEADERS_SRCS)
	$(CXX) $(CXXFLAGS) -o $$@ $($(1)_SRC) $(COMMON_SRCS) $($(1)_EXTRA_SRCS)
endef

# Map the single program to the generic build structure
MAIN_SRC := $(MAIN_SRC)
MAIN_EXTRA_SRCS := $(MAIN_EXTRA_SRCS)
MAIN_EXE := $(MAIN_EXE)

$(foreach prog,$(PROGRAMS),$(eval $(call build_program,$(prog))))

main: $(MAIN_EXE)

clean:
	rm -rf $(BIN_DIR)/*
