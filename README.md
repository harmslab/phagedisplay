# Phage display analysis tools

Computational pipeline for analyzing phage display experiments.

## Quick start

The phage display package will take a set of compressed fastq files, tracing rounds of
enrichment in a phage display experiment, and then identify clusters that enriched
over the experiment.  

After installation, the main pipeline should be accessible by typing `phage` on the command line.

```
phage test_run -f round0.fastq.gz round1.fastq.gz round2.fastq.gz
```

This will create a directory called `test_run` that has four levels of processing.

```
test_run/
    00_raw-fastq-files/
        round0.fastq.gz
        round1.fastq.gz
        round2.fastq.gz
    01_raw-counts/
        bad-counts.pickle
        good-counts.pickle
    02_binding-polynomial/
        human-readable-summary.txt        
    03_clustering/
        clusters/
            cluster_stats.txt
            summary_count.txt
            000_cluster.csv
            000_cluster.pdf
            ....
```

Each processing level will also have `log.txt`, `info.json` and `save.pickle`.  
These files are used by the library and (generally) don't need to be read by a 
human being.

#### raw-fastq-files
This holds the original fastq files.

#### raw-counts
This holds the quality-filtered counts of each sequence.  Each fastq sequence is
passed through a quality control filter that makes sure the 3' end of the sequence
has the proper phage fusion sequence and that there are no stop codons in the 
peptide region of the phage sequence.  The .pickle files are dictionaries that can
be read in by:

```
import pickle
counts = pickle.load(open("good-counts.pickle","rb"))
```

The dictionaries have the structure:

```
counts = {"SEQUENCEHERE":[round0_counts,round1_counts,round2_counts],...}
```

`good-counts.pickle` holds the sequences that passed quality control.
`bad-counts.pickle` holds the sequences that failed quality control.

#### binding-polynomial
Sequences are treated with a binding polynomial, where the relative counts of
each sequence between two rounds (usually rounds 1 and 2) are used to determine
which sequence enrich (log(K) > 0) or deplete (log(K) < 0) during those rounds.
The `human-readable-summary.txt` file holds sequences, enrichments, and counts.

#### clustering
Sequences are clustered by blosum62 distance (by default) using dbscan
(by default).  By default, sequences are only included in the clustering if their
logK values are > 0. The `clusters` folder holds the output of this clustering.  
`cluster_stats.txt` holds information about the clustering itself.  `summary_count.txt`
holds number of sequences in each cluster.  (The `-1` cluster has sequences that 
could not be clustered).  The `xxx_cluster.csv` and `xxx_cluster.pdf` files describe
each cluster.  The `.csv` file holds the sequences; the `.pdf` shows a sequence logo.

The `weight_clusters_by_K.py` script (in the `phagedisplay/util/` directory) can be
used to weight clsuters by the `K` values of all cluster members.
 
## Developers

Git must be installed to clone and contribute to this project

### Setting up for development

1. Clone this repo locally:

```
git clone https://github.com/harmslab/phagedisplay
```

2. Navigate to this directory, and install (softly) this python package with 

```
python setup.py develop
```

### Developing

1. Fork this repository on Github
2. Start a branch locally from master

```
git checkout -B <branch-name>
```

3. Make changes and commit to that branch.

```
git commit -a -m "<commit message>"
```

4. Push to your fork.

```
git push <remote-fork-url-or-alias> <branch-name>
```

5. Pull request into this master repository.

## Users

Clone this repo locally:

```
git clone https://github.com/harmslab/phagedisplay
```

Navigate to this directory, and install this python package with 

```
python setup.py install
```


