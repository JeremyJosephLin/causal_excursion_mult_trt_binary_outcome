#!/bin/bash
########################################################

### Builds the job index
# Create a sequential range
setting_values=`seq 1 4`
nsim=1000
nsetting=871

# Launch the job and then remove the temporarily created qsub file.
for i in $setting_values
    do 
   	sbatch sim_0723.sub $i $nsim $nsetting
done