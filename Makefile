# Targets:
#   make                 solver (build/bin/euler_2D) and renderer (build/bin/render)
#   make test            build and run all tests (regression, physics, snapshot I/O, render smoke test, docs)
#   make benches         build the benchmarks; make bench also runs them
#   make clean           delete build/
#   make clean-results   delete results/ and *.mp4 in the top directory
# OPT=-O3 or BUILD=build-debug can be set on the command line.

# Compiler and flags
CC = gcc
OPT ?= -O2
# No FMA contraction: keeps results bitwise reproducible across optimization levels
CFLAGS = -std=c99 $(OPT) -ffp-contract=off -Wall -Wextra -MMD -MP -I./src
LDLIBS = -lm

# All build outputs live under BUILD (override to keep builds with different OPT apart)
SRC_DIR := src
BUILD ?= build
OBJ_DIR := $(BUILD)/obj
BIN_DIR := $(BUILD)/bin
TARGET := $(BIN_DIR)/euler_2D

# Solver library: everything in src/ except main.c
SRCS := $(shell find $(SRC_DIR) -name '*.c')
LIB_SRCS := $(filter-out $(SRC_DIR)/main.c, $(SRCS))
LIB_OBJS := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(LIB_SRCS))
MAIN_OBJ := $(OBJ_DIR)/main.o

TEST_UTIL_OBJ := $(OBJ_DIR)/tests/test_util.o
TESTS := $(BIN_DIR)/test_regression $(BIN_DIR)/test_physics $(BIN_DIR)/test_io
BENCHES := $(BIN_DIR)/bench $(BIN_DIR)/bench_compat $(BIN_DIR)/bench_layout
TOOLS := $(BIN_DIR)/render

# Default rule: solver and tools
all: $(TARGET) $(TOOLS)

# Link
$(TARGET): $(LIB_OBJS) $(MAIN_OBJ)
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@ $(LDLIBS)

$(BIN_DIR)/test_%: $(OBJ_DIR)/tests/test_%.o $(TEST_UTIL_OBJ) $(LIB_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@ $(LDLIBS)

$(BIN_DIR)/bench: $(OBJ_DIR)/bench/bench.o $(LIB_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@ $(LDLIBS)

$(BIN_DIR)/bench_%: $(OBJ_DIR)/bench/bench_%.o $(LIB_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@ $(LDLIBS)

$(BIN_DIR)/render: $(OBJ_DIR)/tools/render.o $(LIB_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@ $(LDLIBS)

# Compile rules
$(OBJ_DIR)/tests/%.o: tests/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -I./tests -c $< -o $@

$(OBJ_DIR)/bench/%.o: bench/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/tools/%.o: tools/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

tests: $(TESTS)
benches: $(BENCHES)
tools: $(TOOLS)

test: $(TESTS) $(TARGET) $(TOOLS)
	@$(BIN_DIR)/test_regression
	@$(BIN_DIR)/test_physics
	@$(BIN_DIR)/test_io $(BUILD)
	@sh tests/render_smoke.sh $(BIN_DIR) $(BUILD)
	@sh tests/check_docs.sh $(BIN_DIR)

bench: $(BENCHES)
	@set -e; for b in $(BENCHES); do $$b; done

# Removes build outputs only; simulation output is removed by clean-results
clean:
	rm -rf $(BUILD)

clean-results:
	rm -rf results *.mp4

-include $(shell find $(OBJ_DIR) -name '*.d' 2>/dev/null)

.PHONY: all clean clean-results test bench tests benches tools
# Keep objects built through pattern-rule chains (tests, benches, tools)
.SECONDARY:
