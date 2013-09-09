#!/bin/bash

build() {
  dir=$1
  shift
  name=$1
  shift

  pushd $dir
    #../make_knc.sh $name 
    ../make_knc.sh $name  -DISPC_USE_OMP -openmp
  popd
}

build aobench ao
build mandelbrot_tasks mandelbrot
build options options
build rt rt
build volume_rendering volume

