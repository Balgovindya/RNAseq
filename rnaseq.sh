#!/bin/bash

# directories
FASTQ_DIR="/mnt/d/rnaseq/bal"
INDEX_PATH="/mnt/d/rnaseq/index"
OUTPUT_DIR="/mnt/d/rnaseq/bal/aligned"
TRIMMED_DIR="/mnt/d/rnaseq/bal/trimmed"
FEATURES_PATH="/mnt/d/rnaseq/annotation/Physcomitrium.GTF"
ADAPTERS_FILE="/mnt/d/rnaseq/TruSeq3-SE.fa" 
INDEX_PREFIX="$INDEX_PATH/Physcomitrium-genome"

# Create output directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $TRIMMED_DIR

# Step 1: FastQC
#echo "Running FastQC on all FASTQ files..."
#fastqc $FASTQ_DIR/*.fastq.gz -o $FASTQ_DIR/fastqc_reports

# Step 2: Trimmomatic
#echo "Running Trimmomatic for quality trimming..."
# Adjust Trimmomatic options for paired-end reads
#TRIM_OPTS="ILLUMINACLIP:$ADAPTERS_FILE:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"

#for FASTQ_FILE in $FASTQ_DIR/*_R1.fastq.gz; do
    BASENAME=$(basename "$FASTQ_FILE" _R1.fastq.gz)
    R1_FILE="$FASTQ_DIR/${BASENAME}_R1.fastq.gz"
    R2_FILE="$FASTQ_DIR/${BASENAME}_R2.fastq.gz"
    TRIMMED_R1="$TRIMMED_DIR/${BASENAME}_R1_trimmed.fastq.gz"
    TRIMMED_R2="$TRIMMED_DIR/${BASENAME}_R2_trimmed.fastq.gz"
    
    echo "Trimming $R1_FILE and $R2_FILE..."
    trimmomatic PE -phred33 "$R1_FILE" "$R2_FILE" "$TRIMMED_R1" "$TRIMMED_R1".unpaired "$TRIMMED_R2" "$TRIMMED_R2".unpaired "$TRIM_OPTS"

    if [ $? -ne 0 ]; then
        echo "Error trimming $R1_FILE and $R2_FILE"
    else
        echo "Successfully trimmed $R1_FILE and $R2_FILE"
    fi
#done

# Step 3: HISAT2
#echo "Running HISAT2 for alignment..."
#for TRIMMED_R1 in $TRIMMED_DIR/*_R1_trimmed.fastq.gz; do
    BASENAME=$(basename "$TRIMMED_R1" _R1_trimmed.fastq.gz)
    TRIMMED_R2="$TRIMMED_DIR/${BASENAME}_R2_trimmed.fastq.gz"
    OUTPUT_BAM="$OUTPUT_DIR/${BASENAME}.bam"
    
    echo "Aligning $TRIMMED_R1 and $TRIMMED_R2..."
    hisat2 -q --rna-strandness R -x $INDEX_PREFIX -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" | samtools sort -o "$OUTPUT_BAM"

    if [ $? -ne 0 ]; then
        echo "Error aligning $TRIMMED_R1 and $TRIMMED_R2"
    else
        echo "Successfully aligned $TRIMMED_R1 and $TRIMMED_R2 to $OUTPUT_BAM"
    fi
#done

# Step 4: featureCounts
echo "Counting features with featureCounts..."
featureCounts -a $FEATURES_PATH -o $OUTPUT_DIR/feature_counts.txt -p $OUTPUT_DIR/*.bam

# Check if featureCounts ran successfully
if [ $? -ne 0 ]; then
    echo "Error running featureCounts"
else
    echo "Successfully ran featureCounts"
fi

echo "RNA-seq analysis completed."
