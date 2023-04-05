#!/bin/bash

for i in 1 4 8 16
do
  export OMP_NUM_THREADS=$i
  ./aload
  echo " "
done
