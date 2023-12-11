# Makefile

# Compiler and flags
CXX = g++
CXXFLAGS = -g -std=c++17 -Wall

# Target executable
TARGET = raytracer.exe

# Source files
SRCS = Raytracer.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Rule to compile C++ source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Phony target to run the program with any text file
run: $(TARGET)
	@input_file=$$(ls *.txt | head -n 1); \
	if [ -n "$$input_file" ]; then \
		./$(TARGET) "$$input_file"; \
	else \
		echo "No text files found in the current directory."; \
	fi

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
