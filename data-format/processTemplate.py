

import processors

fastq_list = [None,"test-files/r1.fastq","test-files/r2.fastq","test-files/r3.fastq"]

m = processors.MasterProcessor(expt_name="tester")
m.create()

a = processors.FastqListProcessor(expt_name="raw-fastq-files")
m.addProcessor(a)
m.process(file_list=fastq_list)

b = processors.FastqToCountsProcessor(expt_name="raw-counts")
m.addProcessor(b)
m.process()
