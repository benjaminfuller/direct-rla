
#!/bin/bash

N_list=(456000 1000000 1679423 5007954  7963659 16140044)
#N_list=(456000 1679423 16140044)
PLAIN_COMPARISON=$1

# N_list=(16140044 7963659 5007954 1000000 456000)
#N_list=(6000 5000 4000)

for N in "${N_list[@]}"; do
    echo "Running duplicate audit for N=$N"
    python3 Inverse-Paper-Simulation/duplicate_audit_math.py "$N" "$PLAIN_COMPARISON"
done
