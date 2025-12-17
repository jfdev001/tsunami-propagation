# Simple makefile to build a petsc application

# Compiler settings
export PKG_CONFIG_PATH := $(HOME)/petsc-build/v3.22/lib/pkgconfig:$(PKG_CONFIG_PATH)
PETSC_CFLAGS := `pkg-config --cflags petsc`
PETSC_LIBS := `pkg-config --libs petsc`

CC := mpicc
CFLAGS := -Wall -O2 

# Directories
BUILD_DIR = build

# Source and object files
SRCS = main.c
OBJS = $(SRCS:%.c=$(BUILD_DIR)/%.o)

# Final executable
TARGET = $(BUILD_DIR)/main

# Default rule
all: $(TARGET)

# Link final program
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(PETSC_LIBS) -o $@ $^

# Compile .c to build/*.o
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(PETSC_CFLAGS) -c $< -o $@

# Create build directory if missing
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Cleanup
clean:
	rm -rf $(BUILD_DIR)

# Inspect variables 
inspect:
	@echo $(CC)
	@echo $(CFLAGS)
	@echo $(PETSC_CFLAGS)
	@echo $(PETSC_LIBS)

.PHONY: all clean inspect

