# ============================================================
# NERSC / MPI Makefile
# Build current production executables
# ============================================================

CXX := CC

CXXFLAGS := -std=c++17 -O2 -DNDEBUG -DUSE_MPI -MMD -MP
LDFLAGS  :=

BIN_DIR   := bin
BUILD_DIR := build

INCLUDES := \
  -I. \
  -I./include \
  -I./include/eigen-3.4.0 \
  -I./include/nlohmann \
  -I./Common \
  -I./LinearAlgebra \
  -I./MeasureEngines \
  -I./PhysStruct \
  -I./Source \
  -I./Source/Models \
  -I./Source/Parameters \
  -I./Util \
  -I./Util/Config \
  -I./Util/IO

FABRIC_CFLAGS := $(shell pkg-config --cflags libfabric 2>/dev/null)
FABRIC_LIBS   := $(shell pkg-config --libs libfabric 2>/dev/null)

ifeq ($(strip $(FABRIC_LIBS)),)
  LIBFABRIC_LIB := $(firstword \
      $(wildcard /opt/cray/pe/libfabric/*/lib64) \
      $(wildcard /opt/cray/pe/libfabric/*/lib))
  ifneq ($(strip $(LIBFABRIC_LIB)),)
    FABRIC_LIBS := -L$(LIBFABRIC_LIB) -lfabric
    FABRIC_RPATH := -Wl,-rpath,$(LIBFABRIC_LIB) -Wl,-rpath-link,$(LIBFABRIC_LIB)
  endif
endif

CXXFLAGS += $(FABRIC_CFLAGS)
LDFLAGS  += $(FABRIC_LIBS) $(FABRIC_RPATH)

EXEC_SRCS := \
  Executables/calculate_fermi.cpp \
  Executables/calculate_chi.cpp \
  Executables/calculate_orbital_moment.cpp \
  Executables/calculate_berry_curvature.cpp \
  Executables/calculate_magnetic_band.cpp \
  Executables/calculate_magnetic_dos.cpp

EXEC_NAMES := $(notdir $(EXEC_SRCS:.cpp=))
EXEC_BINS  := $(addprefix $(BIN_DIR)/,$(EXEC_NAMES))
EXEC_OBJS  := $(addprefix $(BUILD_DIR)/,$(EXEC_SRCS:.cpp=.o))
DEPS       := $(EXEC_OBJS:.o=.d)

.PHONY: all exec clean print help $(EXEC_NAMES)

all: exec

exec: $(EXEC_BINS)

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN_DIR)/%: $(BUILD_DIR)/Executables/%.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(EXEC_NAMES):
	@$(MAKE) $(BIN_DIR)/$@

print:
	@echo "CXX       = $(CXX)"
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "LDFLAGS   = $(LDFLAGS)"
	@echo "INCLUDES  = $(INCLUDES)"
	@echo "EXEC_BINS = $(EXEC_BINS)"

clean:
	rm -rf $(BIN_DIR) $(BUILD_DIR)

help:
	@echo "Targets:"
	@echo "  make                  build all production executables"
	@echo "  make calculate_fermi  build calculate_fermi only"
	@echo "  make calculate_chi    build calculate_chi only"
	@echo "  make calculate_orbital_moment"
	@echo "  make calculate_berry_curvature"
	@echo "  make calculate_magnetic_band"
	@echo "  make calculate_magnetic_dos"
	@echo "  make print            show compiler flags"
	@echo "  make clean            remove bin/ and build/"

-include $(DEPS)
