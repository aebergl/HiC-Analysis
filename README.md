# HiC-Analysis
## Hi-C tools

Taken from [NGS Learning Hub](https://ngs101.com/how-to-analyze-hi-c-data-for-absolute-beginners-from-raw-reads-to-3d-genome-organization-with-juicer/)

### Installing Juicer
```properties
cd SOFTWARE/
mkdir -p hic_tools
cd hic_tools/
git clone https://github.com/aidenlab/juicer.git
cd juicer
mkdir -p scripts/common
cp CPU/*.* scripts/common
cp CPU/common/* scripts/common
wget https://github.com/aidenlab/Juicebox/releases/download/v2.17.00/juicer_tools_2.17.00.jar
 
mv juicer_tools_2.17.00.jar scripts/common/juicer_tools.jar
```


## Story first example
### Topologically Associating Domains (TAD) Detection

```
PROJ_DIR="/Users/berglund.anders/Documents/USR/STORY/HiC/"
HIC_TOOLS="/Users/berglund.anders/Documents/SOFTWARE/hic_tools/juicer/scripts/common/"


# Run Arrowhead algorithm to detect TADs genome-wide
java -Xmx32g -jar /Users/berglund.anders/Documents/SOFTWARE/hic_tools/juicer/scripts/common/juicer_tools.jar arrowhead \
    -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
    -m 5000 \                           # Minimum resolution for TAD detection
    -r 25000 \                          # Resolution for TAD boundary detection
    -k KR \                             # Normalization method (Knight-Ruiz)
    /Users/berglund.anders/Documents/USR/STORY/HiC/C/inter_30.hic \
    /Users/berglund.anders/Documents/USR/STORY/HiC/C/tads_genome_wide \
# Convert BEDPE TAD output to BED format for Juicebox visualization
awk 'NR>1 {print $1"\t"$2"\t"$6"\tTAD_"NR-1"\t1000\t."}' \
    ~/hic_project/aligned/tads_genome_wide/25000_blocks.bedpe | tail -n +2 > \
    ~/hic_project/aligned/tad_domains_25kb.bed


# C sample  
java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar arrowhead -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX -m 5000 -r 25000 -k KR $PROJ_DIR/C/C_inter_30.hic $PROJ_DIR/C/tads_genome_wide

java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar arrowhead -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX -m 2000 -r 5000 -k KR $PROJ_DIR/C/C_inter_30.hic $PROJ_DIR/C/tads_genome_wide



mv $PROJ_DIR/C/tads_genome_wide/25000_blocks.bedpe $PROJ_DIR/C/tads_genome_wide/C_25000_blocks.bedpe
mv $PROJ_DIR/C/tads_genome_wide/5000_blocks.bedpe $PROJ_DIR/C/tads_genome_wide/C_5000_blocks.bedpe

awk 'NR>1 {print $1"\t"$2"\t"$6"\tTAD_"NR-1"\t1000\t."}' $PROJ_DIR/C/tads_genome_wide/C_25000_blocks.bedpe | tail -n +2 > $PROJ_DIR/C/tads_genome_wide/C_tad_domains_25kb.bed

awk 'NR>1 {print $1"\t"$2"\t"$6"\tTAD_"NR-1"\t1000\t."}' $PROJ_DIR/C/tads_genome_wide/C_5000_blocks.bedpe | tail -n +2 > $PROJ_DIR/C/tads_genome_wide/C_tad_domains_5kb.bed


# H sample
java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar arrowhead -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX -m 5000 -r 25000 -k KR $PROJ_DIR/H/H_inter_30.hic $PROJ_DIR/H/tads_genome_wide

java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar arrowhead -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX -m 2000 -r 5000 -k KR $PROJ_DIR/H/H_inter_30.hic $PROJ_DIR/H/tads_genome_wide

mv $PROJ_DIR/H/tads_genome_wide/25000_blocks.bedpe $PROJ_DIR/H/tads_genome_wide/H_25000_blocks.bedpe
mv $PROJ_DIR/H/tads_genome_wide/5000_blocks.bedpe $PROJ_DIR/H/tads_genome_wide/H_5000_blocks.bedpe

awk 'NR>1 {print $1"\t"$2"\t"$6"\tTAD_"NR-1"\t1000\t."}' $PROJ_DIR/H/tads_genome_wide/H_25000_blocks.bedpe | tail -n +2 > $PROJ_DIR/H/tads_genome_wide/H_tad_domains_25kb.bed

awk 'NR>1 {print $1"\t"$2"\t"$6"\tTAD_"NR-1"\t1000\t."}' $PROJ_DIR/H/tads_genome_wide/H_5000_blocks.bedpe | tail -n +2 > $PROJ_DIR/H/tads_genome_wide/H_tad_domains_5kb.bed

```

### Contact Matrix Extraction

```
#-----------------------------------------------
# STEP 9: Extract contact matrices and perform compartment analysis
#-----------------------------------------------
 
#=============================================
# Extract Contact Matrices for Detailed Analysis
#=============================================

PROJ_DIR="/Users/berglund.anders/Documents/USR/STORY/HiC/"
HIC_TOOLS="/Users/berglund.anders/Documents/SOFTWARE/hic_tools/juicer/scripts/common/"


# Extract intra-chromosomal contact matrix for chromosome 1 at 25kb resolution
java -Xmx32g -jar ~/hic_tools/juicer/scripts/common/juicer_tools.jar dump \
    observed KR \                            # Use KR normalization (recommended)
    ~/hic_project/aligned/inter_30.hic \     # Input Hi-C file
    chr1 chr1 \                              # Extract chr1 vs chr1 interactions
    BP 25000 \                               # Resolution: 25kb bins
    ~/hic_project/aligned/chr1_contacts_25kb.txt  # Output contact matrix
 
# Extract chromosome 1 contact matrix (alternative format for full chromosome)
java -Xmx32g -jar ~/hic_tools/juicer/scripts/common/juicer_tools.jar dump \
    observed KR \                            # Use KR normalization
    ~/hic_project/aligned/inter_30.hic \     # Input Hi-C file
    chr1 \                                   # Extract all chr1 interactions
    BP 25000 \                               # Resolution: 25kb bins
    ~/hic_project/aligned/chr1_25kb.txt      # Output contact matrix

java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar dump observed KR $PROJ_DIR/inter_30.hic 1 1 BP 25000 $PROJ_DIR/Chrs/chr1_contacts_25kb.txt

# Not sure if this is correct
java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar dump norm KR $PROJ_DIR/inter_30.hic 1 BP 25000 $PROJ_DIR/Chrs/chr1_25kb.txt  


```
### Compartment Analysis
```
#=============================================
# A/B Compartment Analysis Using Eigenvector Decomposition
#=============================================
 
# Perform eigenvector analysis to identify A/B compartments
java -Xmx32g -jar ~/hic_tools/juicer/scripts/common/juicer_tools.jar eigenvector \
    KR \                                     # Normalization method
    ~/hic_project/aligned/inter_30.hic \     # Input Hi-C file
    chr1 \                                   # Target chromosome
    BP 100000 \                              # Resolution: 100kb bins (optimal for compartments)
    ~/hic_project/aligned/chr1_eigenvector_100kb.txt  # Output eigenvector values
 
java -Xmx128g -jar $HIC_TOOLS/juicer_tools.jar eigenvector -p KR $PROJ_DIR/inter_30.hic 1 BP 250000 $PROJ_DIR/compartments/chr1_eigenvector_250kb.txt 

# Note: Positive eigenvector values = A compartment (active chromatin)
#       Negative eigenvector values = B compartment (inactive chromatin)
 
#=============================================
# Pearson Correlation Analysis for Compartment Visualization
#=============================================
 
# Calculate Pearson correlation matrix for compartment analysis
java -Xmx32g -jar ~/hic_tools/juicer/scripts/common/juicer_tools.jar pearsons \
    KR \                                     # Normalization method
    ~/hic_project/aligned/inter_30.hic \     # Input Hi-C file
    chr1 \                                   # Target chromosome
    BP 100000 \                              # Resolution: 100kb bins
    ~/hic_project/aligned/chr1_pearsons_100kb.txt  # Output correlation matrix

java -Xmx32g -jar $HIC_TOOLS/juicer_tools.jar/juicer_tools.jar pearsons KR $PROJ_DIR/inter_30.hic 1 BP 250000 $PROJ_DIR/compartments/chr1_pearsons_250kb.txt 
 
# Note: Pearson correlation reveals compartmental organization
#       - Strong positive correlations indicate same compartment type
#       - Negative correlations indicate different compartment types
#       - Creates characteristic plaid pattern in heatmaps
 
#=============================================
# Prepare Files for Juicebox Visualization
#=============================================

```

```
# Detect chromatin loops using HiCCUPS algorithm
java -Xmx32g -jar ~/hic_tools/juicer/scripts/common/juicer_tools.jar hiccups \
    -m 512 \                                # Memory allocation
    -r 5000,10000 \                         # Resolutions for loop detection
    -f 0.1,0.1 \                           # FDR thresholds
    -p 4,2 \                               # Peak calling parameters
    -i 7,5 \                               # Iteration parameters
    -d 20000,20000 \                       # Distance thresholds
    ~/hic_project/aligned/inter_30.hic \
    ~/hic_project/aligned/loops
    
java -Xmx128g -jar $HIC_TOOLS/juicer_tools.jar hiccups -m 512 -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 $PROJ_DIR/inter_30.hic  $PROJ_DIR//loops
```

