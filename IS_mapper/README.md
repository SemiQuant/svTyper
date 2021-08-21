# ISMapper

This program takes paired end Illumina short read sequence data, an IS query of interest and a reference genome or assembly and reports the locations of the IS query in the reference genome or the assembly.
For a more in depth description of how the program works, see 'Method' below.

For support, add an issue to GitHubs issue tracker, or you can email the author at hawkey dot jane at gmail dot com.

The paper can be found as a pre-print on Biorxiv here: http://biorxiv.org/content/early/2015/03/10/016345

## Dependencies
* Python v2.7.5
* BioPython v1.63 - http://biopython.org/wiki/Main_Page
* BWA v0.7.5a - http://bio-bwa.sourceforge.net/
* Samtools v0.1.19 - http://samtools.sourceforge.net/
* Bedtools v2.20.1 - http://bedtools.readthedocs.org/en/latest/content/installation.html
* BLAST+ v2.2.28 - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

## Installation
Install Python and its dependencies first.

Download ISMapper from GitHub:
 ```
git clone https://github.com/jhawkey/IS_mapper/
```

Install with pip:

```
pip install IS_mapper/
```

Testing ISMapper is installed:

```
ismap --version
compiled_table.py -h
```

### Running a test case

To check that all dependencies are working correctly and outputs are as expected:

Download _Streptococcus suis_ P1/7 reads from SRA

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225612/ERR225612_*.fastq.gz
```

Run ISMapper, using the IS query and reference genbank in test/inputs

```
ismap --reads ERR225612_*.fastq.gz --queries IS_mapper/test/inputs/ISSsu3.fasta --typingRef IS_mapper/test/inputs/S_suis_P17.gbk --log --runtype typing --output S_suis
```

There are some example output files in test/example_results to compare your results files to.

## Method

ISMapper finds locations of an IS query in short read data using a series of mapping steps.

First, all reads for an isolate are mapped to the IS query using BWA. The reads we are interested in are the ones whose pairs have mapped, but are unmapped themselves (so therefore must be flanking the IS query). These reads are selected using SamTools.

These unmapped pairs are selected and placed into two groups - left end (flanking the left end of the IS query) and right end (flanking the right end of the IS query).
Each of these groups are indpendently mapped (with BWA) to either a) the reference genome, which may or may not have a particular IS query location present or b) an assembly of the isolate in question.
From this mapping, Bedtools is used to find the depth at each position in the reference genome. We are looking for large peaks of left end and right end reads that indicate a possibly IS query location. Positions where the depth falls below 6 are eliminated from further analysis, as these may be incorrect alignments. (This threshold can be changed using the --cutoff flag, see 'Advanced options for ismap'.)

Overlapping or close regions are merged using the Bedtools merge function to prevent multiple hits representing the same region in the final output file. (These defaults can be changed, see 'Advanced options for ismap'.)

Once peaks have been selected, Bedtools is used to find regions which intersect (indicating a new position that is not present in the reference) and which regions are closest (indicating a position that is probably known in the reference).

These positions are then further analysed and tabulated into the _table.txt file if they are considered to be accurate. Any hits which do not make it into the _table.txt file are moved to _removedHits.txt, showing their position and which bed file they come from (intersect or closest) so they can be investigated further if required.

## Usage

There are two possible options for running ISMapper, depending on what reference genome you would like to compare to. The typing option looks for your IS query locations in your short read data, and compares these locations to a reference genome (which may or may not have that particular location in it).
Alternately, ISMapper can also be used for assembly improvement if the improvement option is chosen. Instead of a reference genome being supplied, the assembly of the short read data can be supplied instead. ISMapper then looks for your IS query locations in the assembly, indicating which contigs may be able to be joined together based on the available IS evidence.

### Typing

Input files:
* short read data in fastq format (can be gzipped)
* IS query/queries of interest in fasta format
* reference genome to compare to in genbank format

Basic usage for ISMapper:

Multiple read sets can be supplied as /path/to/reads/*.fastq.gz. They can also be supplied one after the other, separated by spaces, after --reads. ISMapper will pair the reads together.  
Multiple IS queries can also be supplied, seperated by spaces, after --queries. ISMapper will run queries sequentially in the same output folder.  

`
ismap --reads [isolateA_1.fastq.gz] [isolateA_2.fastq.gz] [isolateB_1.fastq.gz] [isolateB_2.fastq.gz] --queries IS_query.fasta --typingRef reference_genome.gbk --runtype typing --output prefix_out
`

Once ISMapper has finished running, for each isolate there will be multiple output files, the most interesting of which is the *_table.txt file, showing each location in the reference genome where there is a copy of your IS query in your isolate.

Output files:
The final _table.txt file contains the following columns:  
* region - region name  
* orientation - directon of IS in this location (F = forward, R = reverse)  
* x - left most position in the genome where the IS is located (does not indicate orientation)  
* y - right most position in the genome where the IS is located (does not indicate orientation)  
* gap - distance between x and y (small gaps usually indicate the overlap of the left and right ends and usually represent the DR the IS makes when it inserts)  
* call - either Known (in the reference) or Novel (not in the reference). A * next to either of these calls indicates that the call is imprecise; ie: the gap size is larger than expected. A ? next to either of these calls indicates that the call is unconfident; ie: one side (left or right) is high coverage, whilst the other side is low coverage and likely next to a repeat region. The call Possible related IS occurs when the gap between the two hits is approximately the right size for the query IS, however the BLAST results fall below a threshold of 80% ID and 80% coverage.
* %ID - percent match of the sequence between x and y to the IS query if a Known position  
* %Cov - percent coverage of the sequence between x and y to the IS query if a Known position  
* left_gene - information about the left most feature to the IS location (default is locus_tag, gene, product for CDS features or locus_tag and product for tRNA and rRNA features)  
* left_strand - directon of the left most feature to the IS location  
* left_distance - distance of the IS location from the start codon of the left most feature  
* right_gene - information about the right most feature to the IS location (default is locus_tag, gene, product for CDS features or locus_tag and product for tRNA and rRNA features)  
* right_strand - direction of the right most feature to the IS location  
* right_distance - distance of the IS location from the start codon of the right most feature   
* functional_prediction - states 'Gene interrupted' if both the left and right flanking genes are the same, indicating that the IS is in the middle of this feature


The individual _table.txt files for each isolate can be compiled together to generate one large table showing all possible IS query locations in all isolates as well as the reference genome by using the compiled_table script. Options required are shown below, further options are detailed in the section 'Other options for compiled_table'.

`
compiled_table.py --tables *_table.txt --reference_gbk reference_genome.gbk --seq IS_query.fasta --output compiled_table_out.txt
`

This final compiled table has a list of isolates compiled together in the first column, with a header showing the different IS locations.
The top row will always be the reference genome (the same one used in the original analysis).
A - sign indicates that this particular position is not present in the isolate, while a +, * or ? sign indicates that it is present, with the caveats about * (imprecise) and ? (uncofindant) calls as above.
The final seven rows indicate:
* orientation of the IS (either F or R)
* left flanking gene ID (default is locus_tag, or the first in the list set by `--cds`, `--trna` or `--rrna`)
* distance from the IS position to the left flanking gene
* left flanking gene strand (either 1 for forward or -1 for reverse)
* information about the left flanking gene in the feature annotation (default is gene name and product, or the remainder of the options in `--cds`, `--trna` or `--rrna`)
* right flanking gene ID, distance, strand and information


### Improvement

Input files:
* short read sequence data in fastq format (can be gzipped)
* IS query/queries of interest in fasta format
* assemblies for each of your isolates, either in fasta or genbank format

Basic usage:

When supplying assemblies for each isolate, it is vital that they all share the same file name structure, with only the isolate name differentiating them.  

Eg:  
isolateA_assembly.fasta matches to isolateA_1.fastq.gz and isolateA_2.fastq.gz  
isolateB_assembly.fasta matches to isolateB_1.fastq.gz and isolateB_2.fastq.gz  
NOT  
isolateA_assembly.fasta matches to isolateA_1.fastq.gz and isolateA_2.fastq.gz  
isoalteB_contigs.fasta matches to isolateB_1.fastq.gz and isolateB_2.fastq.gz  

Multiple read sets can be supplied one after the other, separated by spaces, after --reads. ISMapper will pair the reads together.  
Multiple IS queries can also be supplied, seperated by spaces, after --queries. ISMapper will run queries sequentially in the same output folder.

`
ismap --reads [isolateA_1.fastq.gz] [isolateA_2.fastq.gz] [isolateB_1.fastq.gz] [isolateB_2.fastq.gz] --queries IS_query1.fasta IS_query2.fasta --assemblies [isolateA_assembly.fasta] [isolateB_assembly.fasta] --assemblyid _assembly --extension .fasta --type fasta --runtype improvement --output prefix_out
`

Once ISMapper has finished running, multiple output files will be generated, the most interesting of which will be the *_table.txt file. The table file contains the names of contigs that have either a left or a right end in them.

## Common Issues

When running ISMapper on typing mode with a reference genome in Genbank format, check to see if all CDS/tRNA/rRNA features have the locus_tag qualifier. This is the default setting for ISMapper when choosing what to print in the final table. If the reference genome doesn't have locus_tag's, ISMapper will throw an error in the final table create step. This can be avoided by changing the setting of `--cds` (or `--trna` or `--rrna`) flags to db_xref gene product (or whatever other identifiers are present in your genbank file).

## Adavnced Options for ismap

`--forward` and `--reverse` are used to designate which reads are forward and reverse only used if NOT inMiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise default is _1, i.e. expect forward reads as sample_1.fastq.gz)

`--cutoff` is used to determine the read depth at a position that will be called as a 'true' peak, and thus make it to the next round of analysis. Default is currently 6 - if you're having issues detecting locations that you feel like should be there, lowering this cutoff may assist in finding them.  

`--min_range` and `--max_range` are used to determine the gap size between a left and right end for calling a known hit. The current default setting is 0.2 and 1.1 (so between 20% and 110% of the actual size - allows for the detection of known truncated hits).  

`--merging` is the distance between two blocks that will be merged by Bedtools. The default is 100. This setting helps remove duplicate peaks in the same region that are not overlapping, but still belong to the same IS query location/

`--min_clip` and `--max_clip` are used to determine the size of soft clipped sections of reads that are including in the final mapping step. The default size is currently 10 bp and 30 bp, which is ideal for reads between 100 - 150bp. To improve the detection of the exact target site duplication size, it is sometimes helpful to increase the size of `--max_clip`. However large values of `--max_clip` can increase the amount of noise in the final mapping step and cause nonsensical results.

`--a`, `--T` and `--t` are flags that are passed to BWA. --a will turn on all alignment reporting in BWA, and --T is used to give an integer mapping score to BWA to determine what alignments are kept. These options may be useful in finding IS query positions that are next to repeated elements as BWA will report all hits for the read not just the best random hit. However, using these options may cause noise and confusion in the final output files. `--t` is used to supply more threads to BWA if required.

`--cds`, `--trna` and `--rrna` are used to specify what qualifiers will be looked for in the reference genbank when determining genes flanking the IS query location. Defaults are locus_tag gene product.

`--igv` turns on the creation of a trackline and hovertext display for the BED file for viewing in IGV.

`--chr_name` specifies the chromosome name for the IGV BED file - must match the genome name loaded into IGV. The default is 'not_specified', which will pull the genome accession from the reference genbank.

`--log` turns on the log file.

`--directory` sets an output directory for the output files (defualt is the directory where ISMapper is being run).

`--temp` turns on keeping the temporary files instead of deleting them once the run has completed.

`--bam` turns on keeping the final sorted and indexed BAM files of flanking reads for comparison against the reference genome (by default these files are deleted to save disk space).


## Other options for compiled_table

`--gap` determines the overlap between nearby positions to be called as the same position in the final compiled tables output. Default is 0 (so positions must be exactly the same in different isolates to compile together), however increasing this number can simplify the final output

`--cds`, `--trna` and `--rrna` are used to specify what qualifiers will be looked for in the reference genbank when determining genes flanking the IS query location. Defaults are locus_tag and product.



## Running ISMapper without installing  

ISMapper can be run directly from its directory without installing it via pip. To do so, ismap.py needs the path to the folder that contains all the scripts supplied to the argument `--path`.

eg:  
`ismap.py --reads x_1.fastq.gz x_2.fastq.gz --queries is_query.fasta --path /path/to/IS_mapper/scripts/`

If using the slurm_ismap.py script, `--path` can be supplied inside `--other_args`, but slurm_ismap.py also requires the path to the actual ismap.py script as well.

eg:  
`slurm_ismap.py --reads x_2.fastq.gz x_2.fastq.gz --queries is_query.fasta --script /path/to/IS_mapper/scripts/ismap.py --other_args "--path /path/to/IS_mapper/scripts"`  

## Running multiple jobs on a SLURM queing system

There is no need to set the --output flag in other_args for this script, it sets it itself using the name of the readset.  
Also, make sure the query sequence is already index with bwa (using `bwa mem query.fasta`) to prevent corruption of the index once multiple jobs start running.  

When running in improvement mode, ensure that --assemblies and --assembliesid is passed to the slurm script, but --extension and --type are passed to the ismap script through --other_args.

```
python slurm_ismap.py -h
usage: slurm_ismap.py [-h] [--walltime WALLTIME] [--memory MEMORY]
                      [--rundir RUNDIR] --script SCRIPT --query QUERY --reads
                      READS [READS ...] [--forward FORWARD]
                      [--reverse REVERSE]
                      [--assemblies ASSEMBLIES [ASSEMBLIES ...]]
                      [--assemblyid ASSEMBLYID] --runtype RUNTYPE
                      [--logprefix LOGPREFIX] [--other_args OTHER_ARGS]

Submit ISMapper jobs to SLURM

optional arguments:
  -h, --help            show this help message and exit
  --walltime WALLTIME   Amount of wall time. Default 1 hr
  --memory MEMORY       Amount of memory (in MB). Default is 16gb
  --rundir RUNDIR       Directory to run in. Default is current directory
  --script SCRIPT       Location of ISMapper script, ismap.py
  --queries QUERIES         Path to IS queries.
  --reads READS [READS ...]
                        Paired end read files in fastq.gz format
  --forward FORWARD     Identifier for forward reads if not in MiSeq format
                        (default _1)
  --reverse REVERSE     Identifier for reverse reads if not in MiSeq format
                        (default _2)
  --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        Contig assemblies, one for each read set (If using
                        improvement option)
  --assemblyid ASSEMBLYID
                        Identifier for assemblies eg: sampleName_contigs
                        (specify _contigs) or sampleName_assembly (specify
                        _assembly). Do not specify extension.
  --runtype RUNTYPE     Runtype for the program, either improvement or typing
  --logprefix LOGPREFIX
                        Creates a prefix for the log file (default is just
                        sample name)
  --other_args OTHER_ARGS
                        String containing all other arguments to pass to
                        ISMapper
```
