# Simple makefile to build a petsc application
.PHONY: all clean inspect

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c

# Directories
BUILD_DIR = build

# Source and object files
SRCS = main.c
OBJS = $(SRCS:%.c=$(BUILD_DIR)/%.o)

# Final executable
TARGET = $(BUILD_DIR)/main

# Petsc variables
include ${PETSC_DIR}/lib/petsc/conf/variables

# Default rule
all: $(TARGET)

# Link final program
# NOTE: CLINKER uses -O0 and -g3 by default... so change to higher optimization
# in production...
$(TARGET): $(OBJS)
	${CLINKER} -o $@ $^ ${PETSC_LIB}

# Compile .c to build/*.o
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	 ${CC} -c $< -o $@ ${PETSC_CC_INCLUDES}

# Create build directory if missing
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Cleanup
clean:
	rm -rf $(BUILD_DIR)

# Inspect variables 
inspect:
	@echo "CC:" $(CC)
	@echo "CFLAGS:" $(CFLAGS)
	@echo "CLINKER:" ${CLINKER}
	@echo "PETSC_LIB:" ${PETSC_LIB}
	@echo "PETSC_CC_INCLUDES:" ${PETSC_CC_INCLUDES}
