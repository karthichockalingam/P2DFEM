-include makefile.local

# Use the MFEM install directory
MFEM_INSTALL_DIR ?= ../mfem-debug
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

EXECUTABLES = batree

OP_SRC_FILES = EChemOperator.cpp
OP_INC_FILES = EChemOperator.hpp
EQ_SRC_FILES = SolidConcentration.cpp \
               ElectrolyteConcentration.cpp \
               SolidPotential.cpp \
               ElectrolytePotential.cpp
EQ_INC_FILES = SolidConcentration.hpp \
               ElectrolyteConcentration.hpp \
               SolidPotential.hpp \
               ElectrolytePotential.hpp \
               Equation.hpp
CF_INC_FILES = ExchangeCurrentCoefficient.hpp \
               ReactionCurrentCoefficient.hpp \
               OpenCircuitPotentialCoefficient.hpp \
               OverPotentialCoefficient.hpp
CL_INC_FILES = LGM50.hpp

SRC_FILES = $(addprefix equations/, $(EQ_SRC_FILES)) \
            $(addprefix operators/, $(OP_SRC_FILES)) \
            constants.cpp batree.cpp
INC_FILES = $(addprefix equations/, $(EQ_INC_FILES)) \
            $(addprefix operators/, $(OP_INC_FILES)) \
            $(addprefix coefficients/, $(CF_INC_FILES)) \
            $(addprefix cells/, $(CL_INC_FILES)) \
            constants.hpp

.PHONY: all clean

all: $(EXECUTABLES)

batree: $(SRC_FILES) $(INC_FILES) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) -I. $(SRC_FILES) -o $@ $(MFEM_LIBS)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean:
	rm -f *.o *~ $(EXECUTABLES)
	rm -rf *.dSYM *.TVD.*breakpoints
