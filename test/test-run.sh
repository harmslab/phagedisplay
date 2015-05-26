#!/bin/bash

echo "Testing fastq import"
../runPipeline.py from-fastq-input -f *.fastq -r rounds.txt -k key-file.txt

echo "Testing pickle import"
../runPipeline.py from-pickle-input -p test.pickle

echo "Testing simulation import (illumina reads)"
../runPipeline.py from-sim-illumina-input -s output.simulation

echo "Testing simulation import (all samples)"
../runPipeline.py from-sim-all-input -s output.simulation -a

echo "Test basic simulation"
../runSimulation.py test-simulation

