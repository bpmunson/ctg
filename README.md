# Compositional and Time-Course aware Genetic interaction analysis (CTG) (Beta)

The CTG package is an analysis framework for quantifying genetic interactions from dual knockout screens.

# Installation

CTG requires the following software

* Python 3.6
* [numpy](https://docs.scipy.org/doc/)
* [scipy](https://docs.scipy.org/doc/)
* [pandas](http://pandas.pydata.org/)
* [pysam](https://pysam.readthedocs.io/en/latest/api.html)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [statsmodels](https://www.statsmodels.org/stable/index.html)
* [configargparse](https://github.com/bw2/ConfigArgParse)

It is recommended to install these dependencies with conda from the [Anaconda](https://conda.io/docs/user-guide/install/download.html) distrubution of Python 3.6.

```bash 
# Optionally create a new virtual environment 
conda create -n <name>
source activate <name>

# Update 
conda update -n base -c defaults conda

# Install dependencies
conda install -c bioconda pysam bowtie2
conda install pandas numpy scipy statsmodels configargparse
```

# Install the CTG package
```bash
pip install git+https://github.com/bpmunson/ctg
```


Analysis Pipeline
-----------------

Analyzing genetic interactions with CTG consists of three main steps: counting constructs, aggregating counts, and scoring interactions.

#. **Count  from Raw Sequencing Reads**.  Given fastq files containing paired end sequencing data and a library defintion specifying the expected guide sequences, tabulate counts of construct pairs by aligning the reads with bowtie2.

```bash
ctg count --help
```

#. **Aggregate Replicate Timecourse Counts**. After construct counts have been quantified per replicate timepoint. The counts must be aggregated into sample level.
```bash
ctg aggregate --help
```

#. **Score genetic interaction**.  Once timepoint counts have been aggregated for each replicate of a sample, genetic interaction scoring is performed.  First, double knockout fitness is estimated as the growth rate in counts for each mutant. Individual gene fitnesses are estimated from the dual construct fitnesses and interaction scores are calculated at the gene level.
```bash 
ctg score --help
```

Getting Started
---------------
See the examples directory for how to get started as well as the command line utility.

```bash
ctg --help
```

References
----------
1. Shen, John Paul, et al. "Combinatorial CRISPR–Cas9 screens for de novo mapping of genetic interactions." Nature methods 14.6 (2017): 573.

