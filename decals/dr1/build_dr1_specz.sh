#!/bin/tcsh
#PBS -S /bin/tcsh
#PBS -N buildspecz
#PBS -j oe 
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -o buildspecz.$JOB_ID
#PBS -A desi
#PBS -q regular

# submit this using:
# qsub build_dr1_specz.sh

echo "build_dr1_specz" | idl
