#!/bin/bash
bin="icpc -O2 -mmic -lpthread -Iobjs -I../intrinsics -I/opt/intel/composer_xe_2013/include $*"
echo $bin
exec $bin

#-DISPC_USE_OMP -openmp \
#-DISPC_USE_TBB_TASK_GROUP -std=c++0x -tbb \
#-DISPC_USE_TBB_PARALLEL -std=c++0x -tbb \
