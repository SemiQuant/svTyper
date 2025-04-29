#!/bin/bash

# Set Java environment
if [ -d "$HOME/miniconda3" ]; then
    export JAVA_HOME="$HOME/miniconda3/lib/jvm"
    export PATH="$JAVA_HOME/bin:$PATH"
    # Ensure proper Java library path
    export LD_LIBRARY_PATH="$JAVA_HOME/lib:$JAVA_HOME/lib/server:$LD_LIBRARY_PATH"
fi

# Check for required tools
check_tools() {
    for tool in bwa samtools gridss Rscript samplot; do
        if ! command -v "$tool" &> /dev/null; then
            echo "Error: $tool is not installed or not in PATH."
            exit 1
        fi
    done
    
    # Verify R packages
    Rscript -e 'if(!all(c("rmarkdown", "knitr", "plotly", "tidyverse", "RColorBrewer", "vcfR") %in% installed.packages()[,"Package"])) { quit(status=1) }'
    if [ $? -ne 0 ]; then
        echo "Error: Required R packages are missing. Please install them using:"
        echo "R -e 'install.packages(c(\"rmarkdown\", \"knitr\", \"plotly\", \"tidyverse\", \"RColorBrewer\", \"vcfR\"))'"
        exit 1
    fi
}

# Default values
a="/$0"; a=${a%/*}; a=${a:-.}; a=${a#/}/; script_dir=$(cd $a; pwd)
threads=4
ram=12
Data="./output"
Ref_name="${script_dir}/NC_000962.3.fasta"
Query_file="${script_dir}/IS6110.fasta"
Gb_file="${script_dir}/NC_000962.3.gb"
chr="NC_000962.3"
gene_start=778990
gene_end=779487
gene_name="Rv0678"


# Usage function
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR -r REFERENCE -q QUERY -g GB_FILE [-t THREADS] [-m RAM] [-c CHR] [-s START] [-e END] [-n GENE_NAME]"
    echo "  -i: Input directory containing paired-end fastq files"
    echo "  -o: Output directory for results"
    echo "  -r: Reference genome file"
    echo "  -q: Query file (IS6110.fasta)"
    echo "  -g: GenBank file (NC_000962.3.gb)"
    echo "  -t: Number of threads (default: 4)"
    echo "  -m: RAM in GB for GRIDSS (default: 8)"
    echo "  -c: Chromosome name (default: NC_000962.3)"
    echo "  -s: Gene start position (default: 1000)"
    echo "  -e: Gene end position (default: 2000)"
    echo "  -n: Gene name (default: target_gene)"
    exit 1
}

# Parse command line arguments
while getopts "i:o:r:q:g:t:m:c:s:e:n:" opt; do
    case $opt in
        i) input_dir="$OPTARG";;
        o) Data="$OPTARG";;
        r) Ref_name="$OPTARG";;
        q) Query_file="$OPTARG";;
        g) Gb_file="$OPTARG";;
        t) threads="$OPTARG";;
        m) ram="$OPTARG";;
        c) chr="$OPTARG";;
        s) gene_start="$OPTARG";;
        e) gene_end="$OPTARG";;
        n) gene_name="$OPTARG";;
        ?) usage;;
    esac
done

# Check required parameters
if [ -z "$input_dir" ] || [ -z "$Data" ] || [ -z "$Ref_name" ] || [ -z "$Query_file" ] || [ -z "$Gb_file" ]; then
    usage
fi


# Check for required tools
check_tools

# Check if reference file exists
if [ ! -f "$Ref_name" ]; then
    echo "Error: Reference genome file '$Ref_name' not found."
    exit 1
fi

# Check if make_dash_data.R exists
if [ ! -f "${script_dir}/make_dash_data.R" ]; then
    echo "Error: R script 'make_dash_data.R' not found in '$script_dir'."
    exit 1
fi

# Check if GenBank file exists
if [ ! -f "$Gb_file" ]; then
    echo "Error: GenBank file '$Gb_file' not found."
    exit 1
fi

# Check if Query file exists
if [ ! -f "$Query_file" ]; then
    echo "Error: Query file '$Query_file' not found."
    exit 1
fi

# Check if required Python modules are installed
check_python_modules() {
    if ! python -c "import Bio" &> /dev/null; then
        echo "Error: Biopython is not installed. Please install it using 'pip install biopython'."
        exit 1
    fi
    if ! python -c "import pkg_resources" &> /dev/null; then
        echo "Error: pkg_resources is not installed. Please install it using 'pip install setuptools'."
        exit 1
    fi
}

# Call the new function to check for Python modules
check_python_modules

# Create output directory if it doesn't exist
mkdir -p "$Data"
echo "Output directory created or already exists: $Data"

# Change to output directory
cd "$Data" || { echo "Failed to change directory to $Data"; exit 1; }

# Replace the existing show_progress function with a simple Bash progress bar
show_progress() {
    local total=$1
    local current=0
    echo "Processing files..."
    while [ $current -lt $total ]; do
        current=$((current + 1))
        echo -ne "\rProgress: $((current * 100 / total))%"
        sleep 1  # Simulate work being done
    done
    echo -e "\nProcessing complete."
}

# Update the output of coverage results to be formatted as a table
Split_gene() {
    local chr="$1"
    local gene_start="$2"
    local gene_end="$3"
    local sample="$4"
    local outdir="$5"
    local tempdir="$6"
    local bam="$7"
    local gene_name="$8"

    echo "Running Split_gene for sample: $sample"
    samtools depth -a -r "$chr":${gene_start}-${gene_end} "${tempdir}/${sample}_sample.splitters.unsorted.bam" | \
        awk -v tmp="${sample}_split" '{print tmp"\t"$0}' | \
        column -t -s $'\t' > "${outdir}/${sample}_${gene_name}_splits.coverage.tsv"
    
    echo "Calculating median split reads on gene..."
    awk ' { a[i++]=$3; } END { print a[int(i/2)]; }' "${outdir}/${sample}_${gene_name}_splits.coverage.tsv" >> \
        "${outdir}/${sample}_${gene_name}_splits.coverage.info.txt"
    
    echo "Counting positions with ≥1 split read..."
    awk '{ if ($3!=0) print $0}' "${outdir}/${sample}_${gene_name}_splits.coverage.tsv" | \
        wc -l >> "${outdir}/${sample}_${gene_name}_splits.coverage.info.txt"
    
    samtools depth -a -r "$chr":${gene_start}-${gene_end} "$bam" | \
        awk -v tmp="${sample}_all" '{print tmp"\t"$0}' | \
        column -t -s $'\t' >> "${outdir}/${sample}_${gene_name}_splits.coverage.tsv"
}

# BWA mapping function
bwa_map() {
    local ref="$1"
    local r1="$2"
    local r2="$3"
    local sample="$4"
    local outdir="$5"
    local threads="$6"
    local tempdir="$7"

    echo "Mapping sample: $sample with BWA..."
    # Index reference if needed
    if [[ ! -e ${ref}.bwt ]]; then
        echo "Indexing reference genome..."
        bwa index "$ref"
    fi

    if [[ ! -e "${ref}.fai" ]]; then
        echo "Creating FASTA index..."
        samtools faidx "$ref"
    fi

    echo "Processing sample: $sample"
    
    # Perform alignment
    bwa mem -M -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina" \
        -t "$threads" "$ref" "$r1" "$r2" > "${tempdir}/${sample}_svs.sam"
    
    # Convert to BAM, sort, and index
    samtools view -b -h -@ "$threads" --reference "$ref" "${tempdir}/${sample}_svs.sam" | \
        samtools sort - -o "${outdir}/${sample}_svs.bam"
    rm "${tempdir}/${sample}_svs.sam"
    samtools index "${outdir}/${sample}_svs.bam"
    samtools flagstat -@ "$threads" "${outdir}/${sample}_svs.bam" > "${outdir}/${sample}_svs_stats"

    # Extract split reads
    samtools view -h "${outdir}/${sample}_svs.bam" | \
        "${script_dir}/extractSplitReads_BwaMem.py" -i stdin | \
        samtools view -Sb - > "${tempdir}/${sample}_sample.splitters.unsorted.bam"


    mkdir -p "${tempdir}/lumpy_tmp"
    samtools index "${tempdir}/${sample}_sample.splitters.unsorted.bam"
}

# Process all paired-end fastq files in the input directory
# First, collect all R1 files and their corresponding R2 files
declare -a r1_files=()
declare -a r2_files=()

for R1 in "$input_dir"/*_R1_001.fastq.gz; do
    if [ -f "$R1" ]; then
        R2="${R1/_R1_/_R2_}"
        if [ -f "$R2" ]; then
            r1_files+=("$R1")
            r2_files+=("$R2")
            echo "Found paired files: $R1 and $R2"  # Debugging output
        else
            echo "Warning: Could not find R2 file for $R1"
        fi
    else
        echo "Warning: No R1 file found matching pattern in $input_dir"
    fi
done

# Check if any paired-end files were found
if [ ${#r1_files[@]} -eq 0 ]; then
    echo "No paired-end fastq files found in $input_dir."
    exit 1
fi

# Show progress for processing files
show_progress ${#r1_files[@]}

# Run ismap for all collected files
echo "Running ismap for all samples..."
ismap_cmd="ismap --reads"

# Add all read pairs to the command
for i in "${!r1_files[@]}"; do
    ismap_cmd+=" \"${r1_files[$i]}\" \"${r2_files[$i]}\""
done

# Add the remaining parameters
ismap_cmd+=" --queries \"$Query_file\""
ismap_cmd+=" --reference \"$Gb_file\""
ismap_cmd+=" --output_dir \"$Data/ismap_results\""

# Execute ismap command
echo "Executing ismap command: $ismap_cmd"
eval "$ismap_cmd"

# Process each sample with BWA, GRIDSS, and split read analysis
for i in "${!r1_files[@]}"; do
    R1="${r1_files[$i]}"
    R2="${r2_files[$i]}"
    
    # Extract sample name
    sample=$(basename "$R1" | sed 's/_R1_001.fastq.gz//')
    
    # Create temporary directory
    Temp="${Data}/${sample}"
    mkdir -p "$Temp"
    
    # Run BWA mapping (which now includes split read extraction)
    bwa_map "$Ref_name" "$R1" "$R2" "$sample" "$Data" "$threads" "$Temp"
    
    # Run split read analysis
    bam="${Data}/${sample}_svs.bam"
    Split_gene "$chr" "$gene_start" "$gene_end" "$sample" "$Data" "$Temp" "$bam" "$gene_name"
    
    # Run GRIDSS
    echo "Running GRIDSS for sample: $sample"
    gridss \
        -r "$Ref_name" \
        -j "${script_dir}/gridss.jar" \
        -o "${Data}/${sample}_gridss.vcf" \
        -a "${Data}/${sample}_gridss.bam" \
        -t $threads \
        "$bam"

    # Generate samplot visualization
    echo "Generating samplot visualization for sample: $sample"
    samplot plot -n "$sample" -b "$bam" -o "${Data}/${sample}_region.png" \
        -c "$chr" -s $((gene_start - 1000)) -e $((gene_end + 1000)) -t "ROI"

    echo "Generating dash data for sample: $sample"
    # Ensure HOME environment variable is set
    if [ -z "$HOME" ]; then
        export HOME=$(eval echo ~$USER)
    fi
    Rscript "${script_dir}/make_dash_data.R" "$Data" $gene_start $gene_end $gene_name "$script_dir"
    rm "$R1" "$R2"
done

echo "Script execution completed."
