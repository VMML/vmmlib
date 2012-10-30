
VMMLIB_UNIT_TESTS =\
    tests/unit_test.cpp\
    tests/unit_test_globals.cpp\
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
    tests/qtucker3_tensor_test.cpp \
    tests/tucker3_exporter_importer_test.cpp \
    tests/cp3_tensor_test.cpp \
    tests/t3_hosvd_test.cpp \
    tests/t3_hooi_test.cpp \
    tests/t3_hopm_test.cpp \
    tests/t3_ihopm_test.cpp \
    tests/t3_ttm_test.cpp \
    tests/matrix_pseudoinverse_test.cpp \
    tests/blas_dgemm_test.cpp \
    tests/blas_dot_test.cpp \
    tests/blas_daxpy_test.cpp \
    tests/tensor4_test.cpp \
    tests/t4_converter_test.cpp \

VMMLIB_UNIT_TESTS_OBJECTS = ${VMMLIB_UNIT_TESTS:%.cpp=%.o}  

CXXFLAGS += -I. -Iinclude -Itests -include stdint.h

# Mac OS X specific stuff 
# on mac we want to use the frameworks, not the unix style libs 
ARCH = $(shell uname)
ifeq "$(ARCH)" "Darwin"
CXXFLAGS += -DVMMLIB_USE_LAPACK
LDFLAGS += -framework Accelerate -fopenmp

else
# Linux specific stuff

CXXFLAGS += -include f2c.h -include f2c_fix.h 

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
  LIBDIR=$(DESTDIR)/usr/lib64
else
  LIBDIR=$(DESTDIR)/usr/lib
endif

CXXFLAGS += -DVMMLIB_USE_LAPACK
LDFLAGS +=  -fopenmp

# adjust libs depending on your LAPACK and BLAS distribution
LIBS += -lblas -llapack -lf2c



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

