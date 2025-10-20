# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -g -Wall
LDFLAGS = -lgtest -lpthread

# Source files
SRC = testing.cpp implementation.cpp
TARGET = test_runner

# Default target
all: $(TARGET)

# Build the test runner
$(TARGET): $(SRC)
    $(CXX) $(CXXFLAGS) $(SRC) $(LDFLAGS) -o $(TARGET)

# Run the test runner
run: $(TARGET)
    ./$(TARGET)

# Install dependencies (Linux only)
install:
    sudo apt update
    sudo apt install -y g++ cmake libgtest-dev
    cd /usr/src/gtest && sudo cmake . && sudo make
    sudo cp /usr/src/gtest/*.a /usr/lib

# Clean up build files
clean:
    rm -f $(TARGET)

# Help message
help:
    @echo "Makefile Usage:"
    @echo "  make          - Build the test runner"
    @echo "  make run      - Run the test runner"
    @echo "  make install  - Install dependencies (Linux only)"
    @echo "  make clean    - Clean up build files"
    @echo "  make help     - Show this help message"
