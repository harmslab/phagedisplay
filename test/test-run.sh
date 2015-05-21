#!/bin/bash

# Should create a directory called "from-fastq-output" containing the processed experimental input
../runPipeline.py from-fastq-output -f *.fastq -r rounds.txt -k key-file.txt
