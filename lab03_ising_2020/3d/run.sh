#!/bin/bash

# Set binary and output filenames
binfile=./ising3d
outputfile=data20

# Set parameters
L=20 # system size
Nwarmup=1000
Nsample=50000
Tmin=1.5
Tmax=6.0
dT=0.1

# Get a random seed
seed=$(od -N3 -t u2 /dev/urandom | awk '{print $2}')

# Run the binary with specified input parameters
$binfile $seed $L $Nwarmup $Nsample $Tmin $Tmax $dT $outputfile
