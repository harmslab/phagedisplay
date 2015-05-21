__description__ = \
"""
Framework for simulating phage display experiments with explicit sampling of 
sequences according to frequency and affinity.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-27"

__all__ = ["masterProcessor","fastqListProcessor","fastqToCountsProcessor"]

from processors.masterProcessor import MasterProcessor
from processors.fastqListProcessor import FastqListProcessor
from processors.fastqToCountsProcessor import FastqToCountsProcessor
