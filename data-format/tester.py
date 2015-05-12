#!/usr/bin/env python3
__description__ = \
"""
Basic tests of read/write parsing functionality for the PhageDisplayExperiment
structure.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import os, shutil

print("Create test environment (test-results directory)")

def createDummyFile(name):
    f = open(name,'w')
    f.write(name)
    f.close()

try: 
    shutil.rmtree("test-results")
except FileNotFoundError:
    pass
os.mkdir('test-results')

os.chdir("test-results")

createDummyFile("test1.pickle")
createDummyFile("test2.pickle")
createDummyFile("test1.fastq")
createDummyFile("test2.fastq")

print("Try to {:s}...".format("import module"),end="")
import phageDisplayExperiment
print(" success.")

print("Try to {:s}...".format("create instance"),end="")
p = phageDisplayExperiment.PhageDisplayExperiment()
print(" success.")

print("Try to {:s}...".format("add arbitrary property"),end="")
p.addProperty("arbitrary",(1,2,3,"test"))
print(" success.")

print("Try to {:s}...".format("add bad property"),end="")
try:
    p.addProperty("date","a bad date entry")
    print(" failure.  Property added!")
except ValueError:
    print(" success.  Property rejected.")

print("Try to {:s}...".format("add files with and without type specified"),end="")
p.addFile("test1.pickle",1)
p.addFile("test2.pickle",27,file_type='pickle')
p.addFile("test1.fastq",1)
p.addFile("test2.fastq",5,file_type='fastq')
print(" success.")

print("Try to {:s}...".format("write out data"),end="")
p.write("test-dir")
print(" success.")

print("Try to {:s}...".format("read data we just wrote out"),end="")
g = phageDisplayExperiment.PhageDisplayExperiment()
g.read('test-dir')
print(" success.")

print("Try to {:s}...".format("write out data we read earlier"),end="")
g.write('test-dir_copied')
print(" success.")

print("Look in test-results directory for details.")
