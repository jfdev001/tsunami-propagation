# Simple makefile to build a petsc application

# Compiler settings
PETSC_PKGCONFIG_PATH = $(HOME)/petsc-build/v3.22/lib/pkgconfig/
CC = mpicc
CFLAGS = -Wall -O2 $(PETSC_CFLAGS)
PETSC_CFLAGS = $(shell pkg-config --with-path=$(PETSC_PKGCONFIG_PATH) --cflags petsc)
LIBS = $(shell pkg-config --with-path=$(PETSC_PKGCONFIG_PATH) --libs petsc)

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
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^

# Compile .c to build/*.o
$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@

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
	@echo $(LIBS)

.PHONY: all clean inspect petsc 

