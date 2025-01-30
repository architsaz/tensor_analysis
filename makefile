# Compiler and flags
CC = gcc
 CFLAGS = -Wall -Wextra -Wshadow -Wconversion -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L -g -fsanitize=address -DDEBUG -Iinclude  -I/dagon1/achitsaz/mylib/include
# CFLAGS = -Wall -Wextra -Wshadow -Wconversion -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L  -g -fsanitize=address -Iinclude
# CFLAGS = -pg -Iinclude # the run program to generate the gmon.out and open it with command gprof program gmon.out > report.txt
# CFLAGS = -O3 -march=native -ftree-vectorize -fopt-info-vec -Wall -Wextra -fstack-protector-strong -D_FORTIFY_SOURCE=2 -Iinclude -I/dagon1/achitsaz/mylib/include
 MYLIB_DIR = /dagon1/achitsaz/mylib/lib
 LIB_DIR = lib
# LDFLAGS = -lm
  LDFLAGS = -L$(MYLIB_DIR) -loperator -lmesh -lCGsolver -lCRSmat -lListDS -lm# -L$(LOC_LIB_DIR) -lListDS  # Linker flags to specify libraries

# Directories
SRC_DIR = src
BIN_DIR = bin

# Target executable
TARGET = $(BIN_DIR)/EXEC_tensor

# Source and object files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BIN_DIR)/%.o)

# Default target
all: $(TARGET)

# Linking step
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Compilation step for each .c file
$(BIN_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean the build
clean:
	rm -rf $(BIN_DIR) 

# Phony targets
.PHONY: all clean

# generate static library by object file in the bin directory
# ar rcs lib/ListDS.a bin/*.o
