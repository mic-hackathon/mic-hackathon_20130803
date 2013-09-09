#!/bin/bash

flag="$1"

function build() {
  if [ "$flag" == "clean" ] ; then
    make clean
    return
  fi
  if [ "$flag" == "run" ] ; then
    ./run.sh
    return
  fi

  make clean
  make
  ./run.sh
}

function run () {
  pushd $1
  build
  popd
}

run aobench
run deferred
run gmres
run mandelbrot
run mandelbrot_tasks
run noise
run options
run perfbench
run rt
run sort
run stencil
run volume_rendering
