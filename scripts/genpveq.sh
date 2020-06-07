#!/usr/bin/env bash

# Generate PV equations

# Constants
AXIS1=1
AXIS2=2
MAX_PV_EQ=16

for i in $AXIS1 $AXIS2; do
    for j in $(eval echo "{0..$MAX_PV_EQ}"); do
        if [[ "$j" == "3" || "$j" == "11" ]]; then
            continue
        fi
        python3 ./scripts/equations.py $i $j
    done
done
