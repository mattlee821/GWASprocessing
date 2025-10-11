# based on https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
# based on https://dougspeed.com/reference-panel/

# 1. Setup directories
DIRECTORY_RAW=/data/GWAS_data/files/references/1000genomes/phase3/raw/
DIRECTORY_PROCESSED=/data/GWAS_data/files/references/1000genomes/phase3/processed/
DIRECTORY_TOOLS=~/tools/

# Create directories if they don't exist
mkdir -p ${DIRECTORY_RAW}
mkdir -p ${DIRECTORY_PROCESSED}
mkdir -p ${DIRECTORY_TOOLS}

# 2. Install LDAK (for stats calculation)
cd ${DIRECTORY_TOOLS}
# Download standalone Linux binary of LDAK
wget https://github.com/dougspeed/LDAK/raw/main/ldak6.1.linux
# Make it executable
chmod +x ldak6.1.linux

# 3. Download raw 1000 Genomes files
cd ${DIRECTORY_RAW}
wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst
wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst
wget https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam
wget https://www.dropbox.com/s/0omyj2tyu7jmmw9/deg1_phase3.king.cutoff.out.id?dl=1
wget https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1
wget https://genetics.ghpc.au.dk/doug/genetic_map_b37.zip

# Unzip genetic map (used later to assign genetic distances)
unzip genetic_map_b37.zip

# Rename files to remove URL artifacts and standardize names
mv phase3_corrected.psam all_phase3.psam
mv deg1_phase3.king.cutoff.out.id?dl=1 deg1_phase3.king.cutoff.out.id
mv deg2_phase3.king.cutoff.out.id?dl=1 deg2_phase3.king.cutoff.out.id

# 4. Decompress pgen and pvar files using PLINK2
module add apps/plink2/2.00a68LM

plink2 --zst-decompress ${DIRECTORY_RAW}all_phase3.pgen.zst > ${DIRECTORY_PROCESSED}all_phase3.pgen
plink2 --zst-decompress ${DIRECTORY_RAW}all_phase3.pvar.zst > ${DIRECTORY_PROCESSED}all_phase3.pvar
cp ${DIRECTORY_RAW}all_phase3.psam ${DIRECTORY_PROCESSED}all_phase3.psam

# 5. Remove related individuals
cd ${DIRECTORY_PROCESSED}
# Remove 2nd-degree relatives (1st-degree can be kept if desired)
plink2 --pfile all_phase3 \
       --remove ${DIRECTORY_RAW}deg2_phase3.king.cutoff.out.id \
       --make-pgen

# 6. Identify super-population sample IDs
# Check unique populations in the psam file (column 5)
awk '{print $5}' ${DIRECTORY_PROCESSED}all_phase3.psam | sort | uniq

# Create a "keep" file for each super-population (column 1 = FID, column 2 = IID)
awk '($5=="AFR"){print 0, $1}' ${DIRECTORY_PROCESSED}all_phase3.psam > AFR.keep
awk '($5=="AMR"){print 0, $1}' ${DIRECTORY_PROCESSED}all_phase3.psam > AMR.keep
awk '($5=="EAS"){print 0, $1}' ${DIRECTORY_PROCESSED}all_phase3.psam > EAS.keep
awk '($5=="EUR"){print 0, $1}' ${DIRECTORY_PROCESSED}all_phase3.psam > EUR.keep
awk '($5=="SAS"){print 0, $1}' ${DIRECTORY_PROCESSED}all_phase3.psam > SAS.keep

# 7. Convert to binary PLINK format per super-population
# This step restricts to autosomal SNPs, removes duplicates and rare variants
SUPER_POPULATION=("AFR" "AMR" "EAS" "EUR" "SAS")

for POPULATION in "${SUPER_POPULATION[@]}"; do
    mkdir -p "${POPULATION}"
    echo "." > exclude.snps  # placeholder for any SNPs to exclude
    
    plink2 --make-bed \
        --out "${POPULATION}/${POPULATION}" \
        --pgen all_phase3.pgen \
        --pvar all_phase3.pvar \
        --psam all_phase3.psam \
        --maf 0.01 \
        --autosome \
        --snps-only just-acgt \
        --max-alleles 2 \
        --rm-dup exclude-all \
        --exclude exclude.snps \
        --keep "${POPULATION}.keep"
done

# Handle the "ALL" population separately (no one is excluded)
POPULATION=ALL
mkdir ${POPULATION}
echo "." > exclude.snps
plink2 --make-bed \
    --out "${POPULATION}/${POPULATION}" \
    --pgen all_phase3.pgen \
    --pvar all_phase3.pvar \
    --psam all_phase3.psam \
    --maf 0.01 \
    --autosome \
    --snps-only just-acgt \
    --max-alleles 2 \
    --rm-dup exclude-all \
    --exclude exclude.snps

# 8. Clean FAM files and generate generic SNP names
SUPER_POPULATION=("AFR" "ALL" "AMR" "EAS" "EUR" "SAS")

for POPULATION in "${SUPER_POPULATION[@]}"; do
    # Add population and sex info to FAM files; replace original IDs with generic names
    awk '(NR==FNR){arr[$1]=$5"_"$6;ars[$1]=$4;next}{$1=$2;$2=arr[$1];$5=ars[$1];print $0}' \
        all_phase3.psam "${POPULATION}/${POPULATION}.fam" \
        > "${POPULATION}/${POPULATION}_clean.fam"
    
    mv "${POPULATION}/${POPULATION}_clean.fam" "${POPULATION}/${POPULATION}.fam"

    # Generate generic SNP names: Chr:BP
    awk '{print $1":"$4, $2}' "${POPULATION}/${POPULATION}.bim" > "${POPULATION}/${POPULATION}.names"
done

# 9. Insert genetic distances using PLINK1.9
module add apps/plink1.9/1.90-b77

for POPULATION in "${SUPER_POPULATION[@]}"; do
    plink1.9 --bfile "${POPULATION}/${POPULATION}" \
        --cm-map "${DIRECTORY_RAW}genetic_map_b37/genetic_map_chr@_combined_b37.txt" \
        --make-bed \
        --out "${POPULATION}/${POPULATION}_ref"
    
    # Replace original bed/bim/fam files with the ones with genetic distances
    mv "${POPULATION}/${POPULATION}_ref.bed" "${POPULATION}/${POPULATION}.bed"
    mv "${POPULATION}/${POPULATION}_ref.bim" "${POPULATION}/${POPULATION}.bim"
    mv "${POPULATION}/${POPULATION}_ref.fam" "${POPULATION}/${POPULATION}.fam"
done

# 10. Calculate minor allele frequencies using LDAK
for POPULATION in "${SUPER_POPULATION[@]}"; do
    ${DIRECTORY_TOOLS}ldak6.1.linux --calc-stats "${POPULATION}/stats" \
                                    --bfile "${POPULATION}/${POPULATION}"
done

# 11. Format final data for downstream use
FILES=/data/GWAS_data/files/
WORK=/data/GWAS_data/work/
PROCESSED=/processed/
mkdir ${WORK}references/1000genomes/
mkdir ${WORK}references/1000genomes/phase3/
DIRECTORY=references/1000genomes/phase3/
SUPER_POPULATION=("AFR" "ALL" "AMR" "EAS" "EUR" "SAS")
for POPULATION in "${SUPER_POPULATION[@]}"; do
    mkdir -p ${WORK}${DIRECTORY}${POPULATION}   # create population directory
    chmod go-rwxs ${WORK}${DIRECTORY}           # restrict access to workspace directory

    # Copy only relevant files (.bed/.bim/.fam and stats) while preserving folder structure
    rsync -av --include="*/" \
              --include="*.bed" --include="*.bim" --include="*.fam" --include="stats*" \
              --exclude="*" \
              "${FILES}${DIRECTORY}${PROCESSED}${POPULATION}/" \
              "${WORK}${DIRECTORY}${POPULATION}/"

    # Set read/execute permissions for group
    chmod g+rx ${WORK}${DIRECTORY}
    chmod g+rx ${WORK}${DIRECTORY}/*
done

