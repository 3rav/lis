# Optional out-of-tree plugin: not linked into the normal LIS library.
# Run after configuring and building LIS.
LIS_ROOT ?= ../..
CC ?= cc
CPPFLAGS += -DHAVE_CONFIG_H -I$(LIS_ROOT)/include -I.
CFLAGS ?= -O2 -std=c99 -Wall -Wextra
SLU_CFLAGS := $(shell pkg-config --cflags superlu)
SLU_LIBS := $(shell pkg-config --libs superlu) -lblas -lm
LIS_LIB := $(LIS_ROOT)/src/.libs/liblis.a

all: benchsuperlu testsuperlu_direct testsuperlu_singular

benchsuperlu: benchsuperlu.c lis_superlu_direct.c lis_superlu_direct.h $(LIS_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SLU_CFLAGS) -o $@ benchsuperlu.c lis_superlu_direct.c $(LIS_LIB) $(SLU_LIBS)

testsuperlu_direct: testsuperlu_direct.c lis_superlu_direct.c lis_superlu_direct.h $(LIS_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SLU_CFLAGS) -o $@ testsuperlu_direct.c lis_superlu_direct.c $(LIS_LIB) $(SLU_LIBS)

testsuperlu_singular: testsuperlu_singular.c lis_superlu_direct.c lis_superlu_direct.h $(LIS_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SLU_CFLAGS) -o $@ testsuperlu_singular.c lis_superlu_direct.c $(LIS_LIB) $(SLU_LIBS)

check: all
	OMP_NUM_THREADS=1 ./testsuperlu_direct
	OMP_NUM_THREADS=1 ./testsuperlu_singular
	OMP_NUM_THREADS=1 ./benchsuperlu $(LIS_ROOT)/test/testmat2.mtx $(LIS_ROOT)/test/benchilutp_rhs_100.mtx

clean:
	rm -f benchsuperlu testsuperlu_direct testsuperlu_singular
.PHONY: all check clean
