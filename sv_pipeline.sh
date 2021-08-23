#!/bin/bash
usage () { #echo -e to read \n
  echo "
Usage Options
  -t|--threads Default:4
  -r|--reference = full path to reference fasta
  -d|--read_dir = full path folder contiaing reads
  -r1|--read_1 = read 1 name
  -r2|--read_2 = read 2 name
  -c|--chr = reference name
  -s|--gene_str = main gene start (Default: 778990 [Rv0678])
  -e|--gene_send = main gene end (Default: 779487 [Rv0678])
  -i|--is_search = search for IS6110 sequences?
  -rg|--Ref_gbk = if --is_search then provide genome gbk
  -is|--is6110 = if --is_search then provide IS6110 (or other) fasta
  
Single Use Options
  -dl|--container = download the singularity container to this path and exit (then run using 'singularity run ~/wgs_pipeline_latest.sif')
  -p|--plots = dir_to_outputs, make plots of the 'Split_gene' function output
  
TODO
  Add mask file (NC_000962.3.nonUniquelyMappable.bed) so we can ignore repetative regions and better detect SVs throught the genome
  Make a single use option for the 'Split_gene' function
  Fix denovo for gridss
"
}


if [ $# == 0 ]
then
    usage
    exit 1
fi


declare_globals () {
    while [[ "$#" -gt 0 ]]
    do
      case $1 in
          -t|--threads)
          threads="$2"
          ;;
          -r|--reference)
          Ref_name="$2"
          ;;
          -d|--read_dir)
          Data="$2"
          ;;
          -r1|--read_1)
          R1="$2"
          ;;
          -r2|--read_2)
          R2="$2"
          ;;
          -c|--chr)
          chr="$2"
          ;;
          -s|--gene_str)
          gene_start="$2"
          ;;
          -e|--gene_send)
          gene_end="$2"
          ;;
          -p|--plots)
          plots="$2"
          ;;
          -i|--is_search)
          is_search="y"
          ;;
          -rg|--Ref_gbk)
          Ref_gbk="$2"
          ;;
          -is|--is6110)
          is6110="$2"
          ;;
          -fr|--forward)
          forward="$2"
          ;;
          -rr|--reverse)
          reverse="$2"
          ;;
      esac
        shift
    done
}

declare_globals "$@"



if [[ ! -z $container ]]
then
    cd "$container"
    singularity pull library://semiquan7/default/sv_pipe:latest
    exit 0
fi


if [[ ! -z $plots ]]
then
    # Get script dir, posix version from stack
    a="/$0"; a=${a%/*}; a=${a:-.}; a=${a#/}/; Script_dir=$(cd $a; pwd)
    
    cd "$plots"
    ${Script_dir}/splt_plots.R "$plots"
    
    exit 0
fi



# set defaults
threads=${threads:-4}
gene_start=${gene_start:-778990}
gene_end=${gene_end:-779487}
chr=${chr:-"NC_000962.3"}
ram=$(echo $threads*2 | bc)
export PERL5LIB=/usr/local/bin/scalpel/
forward=${forward:-"R1_001"}
reverse=${reverse:-"R2_001"}

bwa_map () {
  if [[ -e ${1}.btw ]]
    then
      echo "$1 index exists"
    else
      bwa index "$1"
  fi
  
  if [[ -e "${1/.f*/.fai}" ]]
  then
    echo "$1 fai exists"
  else
    samtools faidx "$1"
  fi
  
  bwa mem -M -R "@RG\\tID:${4}\\tSM:${4}\\tPL:Illumina" -t "$6" "$1" "$2" "$3" > "${7}/${4}_svs.sam"
  samtools view -b -h -@ $6 --reference "$1" "${7}/${4}_svs.sam" | samtools sort - -o "${5}/${4}_svs.bam"
  rm "${7}/${4}_svs.sam"
  samtools flagstat -@ $6 "${5}/${4}_svs.bam" > "${5}/${4}_svs_stats"
  
}


LargeSVs () {
  samtools index "$1"
  #LUMPY
  # Extract the discordant paired-end alignments.
  samtools view -b -F 1294 "$1" > "${2}/${3}_sample.discordants.unsorted.bam"
  # Extract the split-read alignments
  samtools view -h "$1" \
    | extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > "${2}/${3}_sample.splitters.unsorted.bam"

  mkdir "${2}/lumpy_tmp"
  
  samtools index "${2}/${3}_sample.splitters.unsorted.bam"
  samtools index "${2}/${3}_sample.discordants.unsorted.bam"

  lumpyexpress \
    -B "$1" \
    -T "${2}/lumpy_tmp" \
    -S "${2}/${3}_sample.splitters.unsorted.bam" \
    -D "${2}/${3}_sample.discordants.unsorted.bam" \
    -o "${2}/${3}_lumpy_SVs.tmp.vcf"
    

  vcf-sort "${2}/${3}_lumpy_SVs.tmp.vcf" > "${4}/${3}_lumpy_SVs.vcf"


  svtyper \
	  -i "${4}/${3}_lumpy_SVs.vcf" \
	  -B "$1" \
	  -l ${1}.json > "${2}/${3}_lumpy_svtyper_SVs.tmp.vcf"
	  
	 vcf-sort "${2}/${3}_lumpy_svtyper_SVs.tmp.vcf" > "${4}/${3}_lumpy_svtyper_SVs.vcf"

    # -S "${2}/${3}_sample.splitters.unsorted.bam" \
  

    bgzip "${4}/${3}_lumpy_SVs.vcf"
    bgzip "${4}/${3}_lumpy_svtyper_SVs.vcf"
    tabix "${4}/${3}_lumpy_SVs.vcf.gz"
    tabix "${4}/${3}_lumpy_svtyper_SVs.vcf.gz"

}


LargeSV_denovo () {
  #check all large variants from lumy using denovo assembly
  # for each SV in the lumpy file, do a local alignment and then check the coverage?
  zcat "$2" > "${2}.tmp"

  vcfToBedpe \
    -i "${2}.tmp" \
    -o "${7}.sv.bedpe"

  rm "${2}.tmp"

  awk '! /\#/ {print $1"\t"$2"\t"$3"\t"$9}' "${7}.sv.bedpe" > "${7}.sv.bed.tmp.tmp"
  awk '! /\#/ {print $1"\t"$5"\t"$6"\t"$10}' "${7}.sv.bedpe" >> "${7}.sv.bed.tmp.tmp"
  
  # sometime there can be a neg, but not supposed to be
  awk -v OFS='\t' '{print $1, ($2<0)? 0 : $2, $3 }' "${7}.sv.bed.tmp.tmp" > "${7}.sv.bed.tmp"

  #not sure about this - the problem is the bed file and the confidence in the ends, should just make something with pysam that does what i need.
  sortBed -i "${7}.sv.bed.tmp" > "${7}.sv.bed"
  rm "${7}.sv.bed.tmp"
  #not sure about this - the problem is the bed file and the confidence in the ends, should just make something with pysam that does what i need.

  scalpel-discovery --single --bam "$1" --bed "${7}.sv.bed" --ref "$3" --kmer 25 sc \
    --dir ${6} --numprocs $5

  scalpel-export --single --db "${6}/variants.db.dir" --bed "${7}.sv.bed" --ref "$3" \
    --variant-type "all" --max-ins-size 10000 --max-del-size 10000 --min-alt-count 10 \
    --min-del-size 3 --min-ins-size 3 \
    --min-coverage 10 --min-vaf 0.2 > "${4}/${7}_${8}_SVs.denovo.vcf"
   # NOTE: The database.db file can be found in the output directory for the single operation
   bgzip "${4}/${7}_${8}_SVs.denovo.vcf"
   tabix "${4}/${7}_${8}_SVs.denovo.vcf.gz"
}

Split_gene () {
    samtools depth -a -r "$1":${2}-${3} "${6}/${4}_sample.splitters.unsorted.bam" | awk -v tmp="${4}_split" '{print tmp"\t"$0}' > "${5}/${4}_splits.coverage.tsv"
    echo "Median split reads on gene" >> "${5}/${4}_splits.coverage.info.txt"
    awk ' { a[i++]=$3; } END { print a[int(i/2)]; }' "${5}/${4}_splits.coverage.tsv" >> "${5}/${4}_splits.coverage.info.txt"
    echo "Positions with ≥1 split read" >> "${5}/${4}_splits.coverage.info.txt"
    awk '{ if ($3!=0) print $0_}' "${5}/${4}_splits.coverage.info.txt" | wc - l >> "${5}/${4}_splits.coverage.info.txt"
    samtools depth -a -r "$1":${2}-${3} "${7}" | awk -v tmp="${4}_all" '{print tmp"\t"$0}' >> "${5}/${4}_splits.coverage.tsv"
}


command -v bwa >/dev/null 2>&1 || { echo >&2 "I require bwa but it's not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "I require samtools but it's not installed. Aborting."; exit 1; }
command -v lumpyexpress >/dev/null 2>&1 || { echo >&2 "I require lumpy but it's not installed. Aborting."; exit 1; }
command -v lumpy >/dev/null 2>&1 || { echo >&2 "I require lumpy but it's not installed. Aborting."; exit 1; }
command -v scalpel-discovery >/dev/null 2>&1 || { echo >&2 "I require scalpel but it's not installed. Aborting."; exit 1; }
command -v scalpel-export >/dev/null 2>&1 || { echo >&2 "I require scalpel but it's not installed. Aborting."; exit 1; }
command -v svtyper >/dev/null 2>&1 || { echo >&2 "I require svtyper but it's not installed. Aborting."; exit 1; }

 
# Temp="/wynton/home/ernst/jdlim/B107_wgs/tmp"
# R1="TB-B107-056_S1_L001_R1_001.fastq.gz"
# Data="/wynton/home/ernst/jdlim/B107_wgs/out_test"
# Ref_name=""
# threads=1
# chr="gi|448814763|ref|NC_000962.3|"
# gene_start=778990
# gene_end=779487


# R2="${nme}_R2_001.fastq.gz"
sample="${R1/_R1_001.fastq.gz/}"
bam="${Data}/${sample}_svs.bam"
R1="${Data}/${R1}"
R2="${Data}/${R2}"


Temp="${Data}/${sample}"
mkdir "$Temp"

bwa_map "$Ref_name" "$R1" "$R2" "$sample" "$Data" $threads "$Temp"
LargeSVs "$bam" "$Temp" "$sample" "$Data"
Split_gene "$chr" $gene_start $gene_end "$sample" "$Data" "$Temp" "$bam"
LargeSV_denovo "$bam" "${Data}/${sample}_lumpy_svtyper_SVs.vcf.gz" "$Ref_name" "$Data" $threads "$Temp" "$sample" "lumpy"


gridss \
  --reference "$Ref_name" \
  --output "${Data}/${sample}_gridss.vcf.gz" \
  --assembly "${Data}/${sample}_gridss.bam" \
  --threads $threads \
  --workingdir "$Temp" \
  --jvmheap "${ram}g" \
  --jar /usr/bin/gridss.jar \
  "$bam"

# --blacklist <blacklist.bed>
# --jar /usr/bin/gridss.jar

# LargeSV_denovo "$bam" "${Data}/${sample}_gridss.vcf.gz" "$Ref_name" "$Data" $threads "$Temp" "$sample" "gridss"


if [[ ! -z $is_search ]]
then      
  #IS6110 - this should ideally be done to the reference strain closes to the input
    /usr/bin/IS_mapper/scripts/ismap_edited.py \
      --path "/usr/bin/IS_mapper/scripts/" \
      --reads "$R1" "$R2" \
      --queries "$is6110" --typingRef "$Ref_gbk" --runtype typing --output "${sample}.IS_6110" \
      --forward "$forward" \
      --reverse "$reverse" \
      --t $threads
      # --cutoff 5 --max_range 1.2 --max_clip 50 --log --bam
      
      mkdir "${Data}/IS6110_${sample}"
      mv *IS6110* "${Data}/IS6110_${sample}"
fi
