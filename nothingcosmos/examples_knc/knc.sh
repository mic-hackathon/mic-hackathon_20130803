#!/bin/bash
name=`basename $1 .ispc`
echo $name
echo ispc -O2 --emit-c++ -h ${name}_ispc.h --target=generic-16 --c++-include-file=knc.h -I../intrinsics $name.ispc -o ${name}_generic16_knc.cpp
#ispc -O2 --emit-c++ -h ${name}_ispc.h --target=generic-16 --c++-include-file=knc.h -I../intrinsics $name.ispc -o ${name}_generic16_knc.cpp
rm ${name}_generic*.cpp
ispc -O2 --emit-c++ -h ${name}_ispc.h --target=generic-16 --c++-include-file=knc.h -I../intrinsics $name.ispc -o ${name}_generic16_knc.cpp
