# Compiler and flags
CXX = gcc
CXXFLAGS = -std=c99 -Wall -Wextra -I./src

# Directories
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
TARGET := $(BIN_DIR)/euler_2D

# Find all .cpp files in src/
SRCS := $(shell find $(SRC_DIR) -name '*.c')
# Convert src/file.cpp to obj/file.o
OBJS := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# Default rule
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) -o $@

# Compile rule: src/foo/bar.cpp -> obj/foo/bar.o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) results/*.csv results/frames/*.csv video_frames/*.png

.PHONY: all clean
