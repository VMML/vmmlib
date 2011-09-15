
VMMLIB_UNIT_TESTS =\
    tests/unit_test.cpp\
    tests/vector_test.cpp\
    tests/matrix_test.cpp\
    tests/quaternion_test.cpp\
    tests/qr_decomposition_test.cpp\
    tests/svd_test.cpp\
    tests/lapack_svd_test.cpp\
    tests/lapack_linear_least_squares_test.cpp\
    tests/lapack_gaussian_elimination_test.cpp\
    tests/vmmlib_unit_tests_main.cpp\
    tests/lapack_sym_eigs_test.cpp\
    tests/tensor3_test.cpp \
    tests/tensor3_iterator_test.cpp \
    tests/tucker3_tensor_test.cpp \
    tests/cp3_tensor_test.cpp \
    tests/t3_hosvd_test.cpp \
    tests/t3_hooi_test.cpp \
    tests/t3_hopm_test.cpp \
    tests/t3_ihopm_test.cpp \
    tests/matrix_pseudoinverse_test.cpp \
    tests/blas_dgemm_test.cpp \

VMMLIB_UNIT_TESTS_OBJECTS = ${VMMLIB_UNIT_TESTS:%.cpp=%.o}  

CXXFLAGS += -I. -Iinclude 

# on mac we want to use the frameworks, not the unix style libs 
ARCH = $(shell uname)
ifeq "$(ARCH)" "Darwin"
CXXFLAGS += -framework Accelerate -DVMMLIB_USE_LAPACK
LDFLAGS += -framework Accelerate

else
CXXFLAGS += -DVMMLIB_USE_LAPACK 
LDFLAGS +=
LIBS += -lclapack -lf2c

endif

all: vmmlib_unit_tests

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@  

vmmlib_unit_tests: $(VMMLIB_UNIT_TESTS_OBJECTS)
ifeq "$(ARCH)" "Darwin"
	 $(CXX) $(LDFLAGS) $(VMMLIB_UNIT_TESTS_OBJECTS) -o $@  
else
	 $(CXX) $(LDFLAGS) $(VMMLIB_UNIT_TESTS_OBJECTS) -o $@  $(LIBS)
endif

clean:
	rm -rf $(VMMLIB_UNIT_TESTS_OBJECTS) vmmlib_unit_tests

