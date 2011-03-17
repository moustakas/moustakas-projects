#!/bin/bash

# Setup inputs

startn=$1 #Starting process
nproc=$2 #Number of processes to run

for ((i=0;i<$nproc;i+=1)); do
    proc=$[$startn+$i]
    command="batch_ages_gandalf_specfit, /solar, $proc, $nproc"
    outfile="${HOME}/batchlogs/ages_gandalf_specfit_process."$proc
    outerrfile="${HOME}/batchlogs/ages_gandalf_specfit_error."$proc
    echo $command | (nice -19 idl) > $outfile 2> $outerrfile & 
done


