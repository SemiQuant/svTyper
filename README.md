# svTyper

Program to find structural variants from Illumina sequencing. 

# Usage

|     Flag                  |     Explanation                                                                                                                   |     Default    |     Required    |   |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------|----------------|-----------------|---|
|     -t\|--threads         |     Threads to use                                                                                                                |     4          |     N           |   |
|     -r\|--reference       |     Full path to reference fasta                                                                                                  |     NA         |     Y           |   |
|     -d\|--read_dir        |     Full path folder containing reads                                                                                             |     NA         |     Y           |   |
|     -r1\|--read_1         |     Read 1 name                                                                                                                   |     NA         |     Y           |   |
|     -r2\|--read_2         |     Read 2 name                                                                                                                   |     NA         |     Y           |   |
|     -c\|--chr             |     Reference name                                                                                                                |     NA         |     Y           |   |
| -g\|--gene_name     | Main gene name                                                                                                            | Rv0678  | N        |             |            |
|     -s\|--gene_str        |     Main gene start                                                                                                               |     778990     |     N           |   |
|     -e\|--gene_end       |     Main gene end                                                                                                                 |     779487     |     N           |   |
|     -i\|--is_search       |     Search for IS6110 sequences?                                                                                                  |     NA         |     N           |   |
|     -rg\|--Ref_gbk        |     If --is_search then provide genome gbk                                                                                        |     NA         |     N           |   |
|     -is\|--is6110         |     If --is_search then provide IS6110 (or other) fasta                                                                           |     NA         |     N           |   |
|                           |                                                                                                                                   |                |                 |   |
|     Single Use Options    |                                                                                                                                   |     NA         |                 |   |
|     -dl\|--container      |      download the   singularity container to this path and exit (then run using 'singularity run   ~/wgs_pipeline_latest.sif')    |     NA         |     N           |   |
|     -p\|--plots           |      dir_to_outputs make   plots of the 'Split_gene' function output                                                              |     NA         |     N           |   |

## basic usage

Get the singularity container
```
singularity pull library://semiquan7/default/sv_pipe:latest
# or singularity pull library://semiquan7/default/sv_pipe:sha256.ba30c0474c8b1e0784cbcfeebbd7f299e729cfcf1e9d9c76aac405a6af0c0047
```

Get the script and reference files
```
git clone --recursive https://github.com/SemiQuant/svTyper.git
```

Run main
```
singularity run /path/to/container/sv_pipe.sif \
  /path/to/script/sv_pipeline.sh \
  --threads 4 \
  --reference "/path/to/references/NC_000962.3.fasta" \
  --read_dir "/path/to/reads/" \
  --read_1 "/path/to/reads/read_R1_001.fastq.gz" \
  --read_2 "/path/to/reads/read_R2_001.fastq.gz" \
  --is_search \
  --Ref_gbk "/path/to/references/NC_000962.3.gb" \
  --is6110 "/path/to/references/IS6110.fasta"
```

# Outputs

|     Filename                                     |     Description                                                 |
|--------------------------------------------------|-----------------------------------------------------------------|
|     ${sample}_svs.bam                            |     BWA mapped bam                                              |
|     ${sample}_svs_stats                          |     Mapping stats                                               |
|     ${sample}_sample.discordants.unsorted.bam    |     Discordant alignments                                       |
|     ${sample}_sample.splitters.unsorted.bam      |     Split alignments                                            |
|     ${sample}_lumpy_svtyper_SVs.vcf              |     Lumpy structural variants                                   |
|     ${sample}_lumpy_SVs.denovo.vcf               |     Lumpy structural variants denovo corrected                  |
|     ${sample}_gridss.vcf.gz                      |     Gridss structural variants                                  |
|     ${sample}_splits.coverage.info.txt           |     Split read coverage stats over gene of interest             |
|     ${sample}_splits.coverage.tsv                |     Read and split read coverage depth over gene of interest    |
|     ${sample}_IS6110_table.txt                   |     IS6110 insertion sites                                      |
|     splits.html                                  |     Interactive plot of coverage and split read coverage        |


[Example output limited data, update]("https://svtype.netlify.app/")



### TODO

* Add mask file (NC_000962.3.nonUniquelyMappable.bed) so we can ignore repetative regions and better detect SVs throught the genome
* Make a single use option for the 'Split_gene' function
* Fix output naming and folders
* Filtering of vcfs is not correct, it does not account for the length of the variant