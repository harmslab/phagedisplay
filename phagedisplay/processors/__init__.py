__description__ = \
"""
Framework for holding onto experimental and simulated phage display 
experiments and processing them in a consistent way.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-27"

__all__ = ["base",
           "master",
           "fastq_list",
           "fastq_to_counts",
           "pickle_dict",
           #"sequencing_error",
           "regress_enrichment"]

from .base import BaseProcessor
from .master import MasterProcessor
from .fastq_list import FastqListProcessor
from .fastq_to_counts import FastqToCountsProcessor
from .pickle_dict import PickleDictProcessor
#from .sequencing_error import SequencingErrorProcessor
from .regress_enrichment import RegressEnrichmentProcessor
