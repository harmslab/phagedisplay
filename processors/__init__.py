__description__ = \
"""
Framework for simulating phage display experiments with explicit sampling of 
sequences according to frequency and affinity.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-27"

__all__ = ["BaseProcessor","MasterProcessor","FastqListProcessor","FastqToCountsProcessor"]

from .base import BaseProcessor
from .master import MasterProcessor
from .fastq_list import FastqListProcessor
from .fastq_to_counts import FastqToCountsProcessor
