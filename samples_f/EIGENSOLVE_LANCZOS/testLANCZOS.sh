#!/bin/bash

echo '301' > Input.param

for i in 1 4 8 16
do

  cp OPENATI_POLICY_INPUT.0.LANCZOS.$i OPENATI_POLICY_INPUT.0

  export OMP_NUM_THREADS=$i
  if [ $i -lt 5 ]
  then
    numactl -c 0 -i 0 ./lload
  elif [ $i -lt 9 ]
  then
    numactl -c 1,2,3 -i 1,2,3 ./lload
  else
    numactl -c 0,1,2,3 -i 0,1,2,3 ./lload
  fi
  echo ""

done
