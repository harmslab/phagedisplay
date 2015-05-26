#!/bin/bash

echo "Testing fastq import"
../runPipeline.py from-fastq-input -f *.fastq -r rounds.txt -k key-file.txt

echo "Testing pickle import"
../runPipeline.py from-pickle-input -p test.pickle

echo "Test basic simulation (will take ~10 min)"
../runSimulation.py test-simulation -l 8

echo "Testing simulation import (illumina reads)"
../runPipeline.py from-sim-illumina-input -s test-simulation

echo "Testing simulation import (all samples)"
../runPipeline.py from-sim-all-input -s test-simulation -a


