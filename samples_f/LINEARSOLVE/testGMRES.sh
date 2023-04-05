#!/bin/bash

for i in 1 4 8 16
do

  cp OPENATI_POLICY_INPUT.0.GMRES.$i OPENATI_POLICY_INPUT.0

  export OMP_NUM_THREADS=$i
  if [ $i -lt 5 ]
  then
    numactl -c 0 -i 0 ./gload
  elif [ $i -lt 9 ]
  then
    numactl -c 1,2,3 -i 1,2,3 ./gload
  else
    numactl -c 0,1,2,3 -i 0,1,2,3 ./gload
  fi
  echo ""
  mv OPENATI_POLICY_REPORT.0 OPENATI_POLICY_REPORT.0.$i

done
