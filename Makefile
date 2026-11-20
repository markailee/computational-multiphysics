PETSC_DIR ?= /home/lwa/program/petsc-3.21.0
include $(PETSC_DIR)/lib/petsc/conf/variables

CXX := mpic++
TARGET := multiphysics
BUILD_DIR := build
SOURCES := $(wildcard src/main.cpp src/mesh/*.cpp src/fvm/*.cpp src/flow/*.cpp src/acoustic/*.cpp)
OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPENDENCIES := $(OBJECTS:.o=.d)
PROJECT_FLAGS := -std=c++17 -O3 -Wall -Wextra -Wpedantic -Isrc -I$(PETSC_DIR)/include

PYTHON ?= python3

.PHONY: all clean meshes validate validate-turbulence validate-acoustic cylinder-aeolian
all: $(TARGET)

meshes:
	$(PYTHON) cases/generate_meshes.py
	$(PYTHON) cases/generate_acoustic_meshes.py

validate: $(TARGET)
	./$(TARGET) cases/cavity
	./$(TARGET) cases/backstep
	./$(TARGET) cases/cylinder
	$(PYTHON) cases/validate.py
	$(PYTHON) cases/validate_turbulence.py --run
	$(PYTHON) cases/validate_acoustic.py --run

validate-turbulence: $(TARGET)
	$(PYTHON) cases/validate_turbulence.py --run

validate-acoustic: $(TARGET)
	$(PYTHON) cases/validate_acoustic.py --run

cylinder-aeolian: $(TARGET)
	$(PYTHON) cases/run_cylinder_aeolian.py --run

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(PETSC_KSP_LIB)

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(PROJECT_FLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

-include $(DEPENDENCIES)
