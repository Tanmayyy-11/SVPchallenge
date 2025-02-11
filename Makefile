CXX = g++
CFLAGS = -std=c++20 -Wall -Wextra -O2
LDFLAGS = -lgmp -lmpfr -lfplll -lcrypto
INCLUDE = -I/opt/homebrew/opt/openssl/include -I/opt/homebrew/include
LIBS = -L/opt/homebrew/opt/openssl/lib -L/opt/homebrew/lib

# Source files
SRC = $(wildcard *.cpp) rng.c
OBJ = $(SRC:.cpp=.o)
OBJ := $(OBJ:.c=.o)  # Ensure rng.c is also compiled

# Output executable
TARGET = output_file

# Environment variables for GMP and MPFR
export CPLUS_INCLUDE_PATH := /opt/homebrew/opt/mpfr/include:/opt/homebrew/opt/gmp/include:$(CPLUS_INCLUDE_PATH)
export LIBRARY_PATH := /opt/homebrew/opt/mpfr/lib:/opt/homebrew/opt/gmp/lib:$(LIBRARY_PATH)

# Default build rule
all: $(TARGET)

# Compile C++ source files
%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Compile C source files (like rng.c)
%.o: %.c
	gcc -c $(INCLUDE) $< -o $@

# Link everything into the final executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) $(LIBS) -o $(TARGET)

# Clean up generated files
clean:
	rm -f $(OBJ) $(TARGET)
