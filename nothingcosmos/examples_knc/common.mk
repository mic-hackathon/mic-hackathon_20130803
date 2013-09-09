
TASK_CXX=../tasksys.cpp
TASK_LIB=-lpthread  
TASK_OBJ=objs/tasksys.o

#CXX=g++
#CC=gcc
CXX=icpc
CC=icc

#CXXFLAGS=-Iobjs/ -O2
#CCFLAGS=-Iobjs/ -O2
#CXXFLAGS=-Iobjs/ -O2  -DISPC_USE_CILK
#CXXFLAGS=-Iobjs/ -O2  -DISPC_USE_OMP -openmp
#CXXFLAGS=-Iobjs/ -O2  -DISPC_USE_TBB_TASK_GROUP -tbb -std=c++0x
#CXXFLAGS=-Iobjs/ -O2  -DISPC_USE_TBB_PARALLEL_FOR -tbb -std=c++0x
#CXXFLAGS=-Iobjs/ -O2  -DISPC_USE_CILK
#CXXFLAGS=-Iobjs/ -I/usr/include/i386-linux-gnu -O2  -DISPC_USE_PTHREADS_FULLY_SUBSCRIBED
CXXFLAGS=-Iobjs/ -O2 -m64
CCFLAGS=-Iobjs/  -O2 -m64
ISPC=ispc -O2 --arch=x86-64 $(ISPC_FLAGS)

#ISPC=ispc -O2 --arch=x86-64 $(ISPC_FLAGS)
LIBS=-lm $(TASK_LIB) -lstdc++
ISPC_OBJS=$(addprefix objs/, $(ISPC_SRC:.ispc=)_ispc.o $(ISPC_SRC:.ispc=)_ispc_sse2.o \
	$(ISPC_SRC:.ispc=)_ispc_sse4.o $(ISPC_SRC:.ispc=)_ispc_avx.o)
ISPC_HEADER=objs/$(ISPC_SRC:.ispc=_ispc.h)

CPP_OBJS=$(addprefix objs/, $(CPP_SRC:.cpp=.o))
CC_OBJS=$(addprefix objs/, $(CC_SRC:.c=.o))
OBJS=$(CPP_OBJS) $(CC_OBJS) $(TASK_OBJ) $(ISPC_OBJS)

default: $(EXAMPLE)

all: $(EXAMPLE) $(EXAMPLE)-sse4 $(EXAMPLE)-generic16 $(EXAMPLE)-scalar

.PHONY: dirs clean

dirs:
	/bin/mkdir -p objs/

objs/%.cpp objs/%.o objs/%.h: dirs

clean:
	/bin/rm -rf objs *~ $(EXAMPLE) $(EXAMPLE)-sse4 $(EXAMPLE)-generic16

$(EXAMPLE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

objs/%.o: %.cpp dirs $(ISPC_HEADER)
	$(CXX) $< $(CXXFLAGS) -c -o $@

objs/%.o: %.c dirs $(ISPC_HEADER)
	$(CC) $< $(CCFLAGS) -c -o $@

objs/%.o: ../%.cpp dirs
	$(CXX) $< $(CXXFLAGS) -c -o $@

objs/$(EXAMPLE).o: objs/$(EXAMPLE)_ispc.h

objs/%_ispc.h objs/%_ispc.o objs/%_ispc_sse2.o objs/%_ispc_sse4.o objs/%_ispc_avx.o: %.ispc
	$(ISPC) --target=$(ISPC_TARGETS) $< -o objs/$*_ispc.o -h objs/$*_ispc.h

objs/$(ISPC_SRC:.ispc=)_sse4.cpp: $(ISPC_SRC)
	$(ISPC) $< -o $@ --target=generic-4 --emit-c++ --c++-include-file=sse4.h

objs/$(ISPC_SRC:.ispc=)_sse4.o: objs/$(ISPC_SRC:.ispc=)_sse4.cpp
	$(CXX) -I../intrinsics -msse4.2 $< $(CXXFLAGS) -c -o $@

$(EXAMPLE)-sse4: $(CPP_OBJS) objs/$(ISPC_SRC:.ispc=)_sse4.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

objs/$(ISPC_SRC:.ispc=)_generic16.cpp: $(ISPC_SRC)
	$(ISPC) $< -o $@ --target=generic-16 --emit-c++ --c++-include-file=generic-16.h

objs/$(ISPC_SRC:.ispc=)_generic16.o: objs/$(ISPC_SRC:.ispc=)_generic16.cpp
	$(CXX) -I../intrinsics $< $(CXXFLAGS) -c -o $@

$(EXAMPLE)-generic16: $(CPP_OBJS) objs/$(ISPC_SRC:.ispc=)_generic16.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

objs/$(ISPC_SRC:.ispc=)_scalar.o: $(ISPC_SRC)
	$(ISPC) $< -o $@ --target=generic-1

$(EXAMPLE)-scalar: $(CPP_OBJS) objs/$(ISPC_SRC:.ispc=)_scalar.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
