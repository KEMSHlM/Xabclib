#!/bin/bash

  export OMP_NUM_THREADS=16
  numactl -c 0,1,2,3 -i 0,1,2,3 ./gload

