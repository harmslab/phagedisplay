

import experiments

fastq_list = [None,"r1.fastq","r2.fastq","r3.fastq"]

t = experiments.TotalContainer()
t.create(expt_name="test")

t.addSubContainer(experiments.FastqListContainer,
                  expt_name="raw-fastq-files")
t.process(file_list=fastq_list)

t.addSubContainer(experiments.FastqToCountsContainer,
                  expt_name="raw-counts")
t.process()
