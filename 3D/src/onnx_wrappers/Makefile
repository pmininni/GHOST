# ==========================================
# Compiler Setup
# ==========================================
CC       := mpicc

# Directory structure
SRC_DIR  := src

# ==========================================
# Flags & Path Configurations
# ==========================================
# Path to ONNX Runtime
ONNX_DIR := /home/Ning.Wang/onnxruntime

CFLAGS   := -Wall -Wextra -O3 -I$(ONNX_DIR)/include -I/usr/lib/gcc/x86_64-redhat-linux/11/include

# Default Rule
all: onnx_wrapper.o

# ==========================================
# Compilation Rules
# ==========================================
# Rule for compiling C source files from $(SRC_DIR) into object files
%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# ==========================================
# Clean & Helper Commands
# ==========================================
clean:
	rm -f *.o

.PHONY: all clean
