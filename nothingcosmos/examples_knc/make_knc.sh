#!/bin/bash

name=`basename $1 .ispc`
shift
../knc.sh $name.ispc
../compile_knc.sh *.cpp ../tasksys.cpp -o ${name}_knc.out $*
