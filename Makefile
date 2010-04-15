
VMMLIB_UNIT_TESTS =\
    tests/unit_test.cpp\
    tests/vector_test.cpp\
    tests/matrix_test.cpp\
    tests/quaternion_test.cpp\
    tests/qr_decomposition_test.cpp\
    tests/svd_test.cpp\
    tests/lapack_svd_test.cpp\
    tests/lapack_linear_least_squares_test.cpp\
    tests/vmmlib_unit_tests_main.cpp\

VMMLIB_UNIT_TESTS_OBJECTS = ${VMMLIB_UNIT_TESTS:%.cpp=%.o}  

CXXFLAGS += -I. -Iinclude 

# on mac we want to use the frameworks, not the unix style libs 
ARCH = $(shell uname)
ifeq "$(ARCH)" "Darwin"
CXXFLAGS += -framework Accelerate -DVMMLIB_USE_LAPACK
LDFLAGS += -framework Accelerate

else

endif

all: vmmlib_unit_tests

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@  

vmmlib_unit_tests: $(VMMLIB_UNIT_TESTS_OBJECTS)
	 $(CXX) $(LDFLAGS) $(VMMLIB_UNIT_TESTS_OBJECTS) -o $@  

clean:
	rm -rf $(VMMLIB_UNIT_TESTS_OBJECTS) vmmlib_unit_tests

