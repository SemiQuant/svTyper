# svTyper

A tool for structural variant (SV) analysis and visualization.

## Description

svTyper is a comprehensive tool for analyzing and visualizing structural variants in genomic data. It provides detailed analysis of different types of structural variants including deletions, insertions, duplications, and translocations.

## Prerequisites

Before installing svTyper, ensure you have the following tools installed:

### Required Software

* Python 3.7 or higher
* R 4.0 or higher
* BWA (Burrows-Wheeler Aligner)
* SAMtools
* BCFtools

### Required R packages

* ggplot2
* dplyr
* tidyr
* gridExtra

### Required Python packages

* pandas
* numpy
* pysam
* click
* tabulate

## Installation

```bash
# Clone the repository
git clone https://github.com/SemiQuant/svTyper.git

# Change to project directory
cd svTyper

# Install Python dependencies
pip install -r requirements.txt

# Install R dependencies
R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr', 'gridExtra'))"
```

## Usage

```bash
svTyper --input input.bam --output results/
```

### Options

| Option | Description |
|--------|-------------|
| --input | Input BAM/CRAM/VCF file |
| --output | Output directory |
| --threads | Number of threads to use |
| --verbose | Show detailed progress |

## Output

The tool generates:

* Detailed SV analysis reports in tabular format
* Visualization plots for SV distribution
* Quality metrics for detected variants