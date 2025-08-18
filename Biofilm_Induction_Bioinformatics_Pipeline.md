# SARGASSUM-DERIVED BIOFILM INDUCTION ANALYSIS
# OVERVIEW
1. Adapter trimming and quality filtering - <i>BBDuk</i>
2. Co-assemble QC reads into contig - <i>MEGAHIT</i>
3. Mapping reads to co-assembled contigs for binning - <i>Bowtie2</i>
4. Bin bMAGs - <i>CONCOCT, MaxBin 2, MetaBAT 2</i>
5. Refine, Re-assmeble, Quantify, Classify bMAGs - <i>MetaWRAP, GTDB-tk</i>
6. bMAG metabolic annotation and % ANI tree- <i>Anvi'o</i>
7. Identify and classify viral sequences and prophages - <i>geNomad</i>
8. Extend viral sequences - <i>COBRA</i>
9. Re-identify and classify extended viral sequences - <i>geNomad</i>
10. Dereplicate viral sequences - <i>Virathon</i>
11. Quality assessment of viral sequences and additional prophage identification - <i>CheckV</i>
12. Prophage induction analysis - <i>smalt, samtools, PropagAtE</i>
13. Linking bMAGs to prophages - <i>geNomad, blastn</i>
14. Assign viral infection strategies - <i>DeepPL</i>
15. Viral abundances - <i>smalt, samtools</i>
16. Viral host prediction - <i>iPHoP</i>
17. Viral gene annotations - <i>MetaCerberus</i>
18. AphA protein analysis 

# 1. Adapter trimming and quality filtering
```bash
BBDuk Documentation: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
```
```bash
# QC Paired Reads - [BBDuk]
for f in *_R1_001.fastq.gz; do 
  name=$(echo ${f} | sed 's/_R1_001.fastq.gz//g') 
  bbduk.sh -Xmx512m -da \
  in1=${name}_R1_001.fastq.gz in2=${name}_R2_001.fastq.gz \
  out1=${name}_R1_001_out1.fastq out2=${name}_R2_001_out1.fastq \
  ktrim=rl k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=30 \
  ref= /path/to/adapters.fa 
  bbduk.sh -Xmx512m -da \
  in1=${name}_R1_001_out1.fastq in2=${name}_R2_001_out1.fastq \
  out1=${name}_R1_001_out2.fastq out2=${name}_R2_001_out2.fastq maq=30 
  bbduk.sh -Xmx512m -da \
  in1=${name}_R1_001_out2.fastq in2=${name}_R2_001_out2.fastq \
  out1=${name}_R1_001_noent.fastq out2=${name}_R2_001_noent.fastq \
  ref=/path/to/phix174_ill.ref.fa.gz k=31 hdist=1 
done 
```
```bash
FastQC Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
```
```bash
# Visualize QC Reads - [FastQC]
for f in *_R1_001_noent.fastq; do
  name=$(echo ${f} | sed 's/_R1_001_noent.fastq//g') 
  fastqc ${name}_R1_001_noent.fastq
done
```

# 2. Co-assemble QC reads into contig
```bash
MEGAHIT Documentation: https://github.com/voutcn/megahit
```
```bash
# Co-assemble reads by the four different sample types - [MEGAHIT]
# MC-PreInd (n=5)
megahit --presets meta-large \
-1 MC-PreInd-BF1_R1_001_noent.fastq,MC-PreInd-BF2_R1_001_noent.fastq,MC-PreInd-BF3_R1_001_noent.fastq,MC-PreInd-BF4_R1_001_noent.fastq,MC-PreInd-BF5_R1_001_noent.fastq \
-2 MC-PreInd-BF1_R2_001_noent.fastq,MC-PreInd-BF2_R2_001_noent.fastq,MC-PreInd-BF3_R2_001_noent.fastq,MC-PreInd-BF4_R2_001_noent.fastq,MC-PreInd-BF5_R2_001_noent.fastq \
--min-contig-len 1000 -m 0.85 -o /path/to/output -t 64
# rename contig file
cp final_contigs.fa MC-PreInd_megahit_contigs.fa

# CT-PreInd (n=5)
megahit --presets meta-large \
-1 CT-PreInd-BF1_R1_001_noent.fastq,CT-PreInd-BF2_R1_001_noent.fastq,CT-PreInd-BF3_R1_001_noent.fastq,CT-PreInd-BF4_R1_001_noent.fastq,CT-PreInd-BF5_R1_001_noent.fastq \
-2 CT-PreInd-BF1_R2_001_noent.fastq,CT-PreInd-BF2_R2_001_noent.fastq,CT-PreInd-BF3_R2_001_noent.fastq,CT-PreInd-BF4_R2_001_noent.fastq,CT-PreInd-BF5_R2_001_noent.fastq \
--min-contig-len 1000 -m 0.85 -o /path/to/output  -t 64
# rename contig file 
cp final_contigs.fa CT-PreInd_megahit_contigs.fa

# MC-12PostInd (n=5)
megahit --presets meta-large \
-1 MC-12PostInd-BF1_R1_001_noent.fastq,MC-12PostInd-BF2_R1_001_noent.fastq,MC-12PostInd-BF3_R1_001_noent.fastq,MC-12PostInd-BF4_R1_001_noent.fastq,MC-12PostInd-BF5_R1_001_noent.fastq \
-2 MC-12PostInd-BF1_R2_001_noent.fastq,MC-12PostInd-BF2_R2_001_noent.fastq,MC-12PostInd-BF3_R2_001_noent.fastq,MC-12PostInd-BF4_R2_001_noent.fastq,MC-12PostInd-BF5_R2_001_noent.fastq \
--min-contig-len 1000 -m 0.85 -o /path/to/output -t 64
# rename contig file
cp final_contigs.fa MC-12PostInd_megahit_contigs.fa

# CT-12PostInd (n=5)
megahit --presets meta-large \
-1 CT-12PostInd-BF1_R1_001_noent.fastq,CT-12PostInd-BF2_R1_001_noent.fastq,CT-12PostInd-BF3_R1_001_noent.fastq,CT-12PostInd-BF4_R1_001_noent.fastq,CT-12PostInd-BF5_R1_001_noent.fastq \
-2 CT-12PostInd-BF1_R2_001_noent.fastq,CT-12PostInd-BF2_R2_001_noent.fastq,CT-12PostInd-BF3_R2_001_noent.fastq,CT-12PostInd-BF4_R2_001_noent.fastq,CT-12PostInd-BF5_R2_001_noent.fastq \
--min-contig-len 1000 -m 0.85 -o /path/to/output -t 25
# rename contig file 
cp final_contigs.fa CT-12PostInd_megahit_contigs.fa
```

# 3. Mapping reads to contigs for binning
```bash
Bowtie2 Documentation: https://bowtie-bio.sourceforge.net/bowtie2/index.shtml
```
```bash
# Make a Bowtie2 database for each co-assembled contig file - [Bowtie2]
# (MC-PreInd, CT-PreInd, MC-12PostInd, CT-12PostInd)
for f in *_megahit_contigs.fa; do
  name=$(echo ${f} | sed 's/_megahit_contigs.fa//g')
  bowtie2-build ${name}_megahit_contigs.fa ${name}_bowtieDB
done

# Map reads to co-assembled contigs to get coverage data by sample type - [Bowtie2]
# Map MC-PreInd reads (n=5) to MC-PreInd bowtie database
for f in *_R1_001_noent.fastq; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  bowtie2 --threads 20 -x MC-PreInd_bowtieDB \
  -1 ${sample}_R1_001_noent.fastq -2 ${sample}_R2_001_noent.fastq -S ${sample}.sam
  samtools view -S -b ${sample}.sam | samtools sort > ${sample}.sorted.bam
  samtools index ${sample}.sorted.bam
  rm ${sample}.sam
done

# Map CT-PreInd reads (n=5) to CT-PreInd bowtie database
for f in *_R1_001_noent.fastq; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  bowtie2 --threads 20 -x CT-PreInd_bowtieDB \
  -1 ${sample}_R1_001_noent.fastq -2 ${sample}_R2_001_noent.fastq -S ${sample}.sam
  samtools view -S -b ${sample}.sam | samtools sort > ${sample}.sorted.bam
  samtools index ${sample}.sorted.bam
  rm ${sample}.sam
done

# Map MC-12PostInd reads (n=5) to CT-12PostInd bowtie database
for f in *_R1_001_noent.fastq; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  bowtie2 --threads 20 -x MC-12PostInd_bowtieDB \
  -1 ${sample}_R1_001_noent.fastq -2 ${sample}_R2_001_noent.fastq -S ${sample}.sam
  samtools view -S -b ${sample}.sam | samtools sort > ${sample}.sorted.bam
  samtools index ${sample}.sorted.bam
  rm ${sample}.sam
done

# Map CT-12PostInd reads (n=5) to CT-12PostInd bowtie database
for f in *_R1_001_noent.fastq; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  bowtie2 --threads 20 -x CT-12PostInd_bowtieDB \
  -1 ${sample}_R1_001_noent.fastq -2 ${sample}_R2_001_noent.fastq -S ${sample}.sam
  samtools view -S -b ${sample}.sam | samtools sort > ${sample}.sorted.bam
  samtools index ${sample}.sorted.bam
  rm ${sample}.sam
done
```

# 4. Bin bMAGs
```bash
CONCOCT Documentation: https://github.com/BinPro/CONCOCT
MaxBin2 Documentation: https://sourceforge.net/projects/maxbin/
MetaBAT2 Documentation: https://bitbucket.org/berkeleylab/metabat/src/master/
```
```bash
# Generate depth and coverage files 
for f in *.bam; do
  sample=$(echo ${f} | sed 's/.bam//g')
  jgi_summarize_bam_contig_depths --outputDepth ${sample}.alignment.sorted.depth.txt ${sample}.sorted.bam
  pileup.sh in=${sample}.sorted.bam out=${sample}.alignment.cov.txt overwrite=true
  awk '{print $1"\t"$5}' ${sample}.alignment.cov.txt | grep -v '^#' > ${sample}.alignment.abundance.txt
done
# Move all bam files and generated coverage files into their own directories for 
# CT-PreInd, MC-PreInd, CT-12PostInd, and MC-12PostInd

# Binning - [CONCOCT]
# MC-PreInd
cut_up_fasta.py MC-PreInd_megahit_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
# generate the samp-names file by listing out every bam file in the MC_PreInd directory
# remove the path and the .sorted.bam to get samples names 
ls MC-PreInd/*.bam | sed 's/.sorted.bam//g' | sed 's/MC-PreInd\///g' > samp-names.tsv
BAMS=`ls MC-PreInd/*.bam | python -c 'import sys; print(" ".join([x.strip() for x in sys.stdin.readlines()]))'`
concoct_coverage_table.py contigs_10K.bed --samplenames samp-names.tsv $BAMS > coverage_table.csv
# run concoct
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.csv -b /path/to/output
# Merge subcontig clustering into original contig clustering
merge_cutup_clustering.py clustering_gt1000.csv > clustering_merged.csv
# Extract bins as individual FASTA
extract_fasta_bins.py MC-PreInd_megahit_contigs.fa clustering_merged.csv \
          --output_path /path/to/output
# CT-PreInd
cut_up_fasta.py CT-PreInd_megahit_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
ls CT-PreInd/*.bam | sed 's/.sorted.bam//g' | sed 's/CT-PreInd\///g' > samp-names.tsv
BAMS=`ls MC-PreInd/*.bam | python -c 'import sys; print(" ".join([x.strip() for x in sys.stdin.readlines()]))'`
concoct_coverage_table.py contigs_10K.bed --samplenames samp-names.tsv $BAMS > coverage_table.csv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.csv -b /path/to/output
merge_cutup_clustering.py clustering_gt1000.csv > clustering_merged.csv
extract_fasta_bins.py CT-PreInd_megahit_contigs.fa clustering_merged.csv \
          --output_path /path/to/output
# MC-12PostInd
cut_up_fasta.py MC-12PostInd_megahit_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
ls MC-12PostInd/*.bam | sed 's/.sorted.bam//g' | sed 's/CT-PreInd\///g' > samp-names.tsv
BAMS=`ls MC-12PostInd/*.bam | python -c 'import sys; print(" ".join([x.strip() for x in sys.stdin.readlines()]))'`
concoct_coverage_table.py contigs_10K.bed --samplenames samp-names.tsv $BAMS > coverage_table.csv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.csv -b /path/to/output
merge_cutup_clustering.py clustering_gt1000.csv > clustering_merged.csv
extract_fasta_bins.py MC-12PostInd_megahit_contigs.fa clustering_merged.csv \
          --output_path /path/to/output
# CT-12PostInd
cut_up_fasta.py CT-12PostInd_megahit_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
ls CT-12PostInd/*.bam | sed 's/.sorted.bam//g' | sed 's/CT-PreInd\///g' > samp-names.tsv
BAMS=`ls CT-12PostInd/*.bam | python -c 'import sys; print(" ".join([x.strip() for x in sys.stdin.readlines()]))'`
concoct_coverage_table.py contigs_10K.bed --samplenames samp-names.tsv $BAMS > coverage_table.csv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.csv -b /path/to/output
merge_cutup_clustering.py clustering_gt1000.csv > clustering_merged.csv
extract_fasta_bins.py CT-12PostInd_megahit_contigs.fa clustering_merged.csv \
          --output_path /path/to/output
# 
```
```bash
# Binning - [MaxBin2]
# MC-PreInd
ls MC-PreInd/*.alignment.abundance.txt > MC-PreInd/abund.list.txt
run_MaxBin.pl -thread 4 -contig MC-PreInd_megahit_contigs.fa -out /path/to/output -abund_list MC-PreInd/abund.list.txt
# CT-PreInd
ls CT-PreInd/*.alignment.abundance.txt > CT-PreInd/abund.list.txt
run_MaxBin.pl -thread 4 -contig CT-PreInd_megahit_contigs.fa -out /path/to/output -abund_list CT-PreInd/abund.list.txtt
# MC-12PostInd
ls MC-12PostInd/*.alignment.abundance.txt > MC-12PostInd/abund.list.txt
run_MaxBin.pl -thread 4 -contig MC-12PostInd_megahit_contigs.fa -out /path/to/output -abund_list MC-12PostInd/abund.list.txt
# CT-12PostInd
ls CT-12PostInd/*.alignment.abundance.txt > CT-12PostInd/abund.list.txt
run_MaxBin.pl -thread 4 -contig CT-12PostInd_megahit_contigs.fa -out /path/to/output -abund_list CT-12PostInd/abund.list.txt
```
```bash
# Binning - [MetaBAT2]
# MC-PreInd
jgi_summarize_bam_contig_depths --outputDepth MC-PreInd/depth.txt MC-PreInd/*.bam
metabat2 -m 1500 -t 24 -i MC-PreInd_megahit_contigs.fa -a MC-PreInd/depth.txt -o /path/to/output
# CT-PreInd
jgi_summarize_bam_contig_depths --outputDepth CT-PreInd/depth.txt CT-PreInd/*.bam
metabat2 -m 1500 -t 24 -i CT-PreInd_megahit_contigs.fa -a CT-PreInd/depth.txt -o /path/to/output
# MC-12PostInd
jgi_summarize_bam_contig_depths --outputDepth MC-12PostInd/depth.txt MC-12PostInd/*.bam
metabat2 -m 1500 -t 24 -i MC-12PostInd_megahit_contigs.fa -a MC-12PostInd/depth.txt -o /path/to/output
# CT-12PostInd
jgi_summarize_bam_contig_depths --outputDepth CT-12PostInd/depth.txt CT-12PostInd/*.bam
metabat2 -m 1500 -t 24 -i CT-12PostInd_megahit_contigs.fa -a CT-12PostInd/depth.txt -o /path/to/output
```

# 5. Refine, re-assmeble, quantify, classify, quality bMAGs
```bash
MetaWRAP Documentation: https://github.com/bxlab/metaWRAP
GTDB-tk Documentation: https://github.com/Ecogenomics/GTDBTk
CheckM2 Documentation: https://github.com/chklovski/CheckM2
```
```bash
# Bin refinement - [MetaWRAP]
# run bin_refinement for each sample type (CT-PreInd, MC-PreInd, CT-12PostInd, MC-12PostInd)
metawrap bin_refinement -o /path/to/output -c 45 -x 10 -t 16 \
                        -A /path/to/metabat2_bins/ \
                        -B /path/to/maxbin2_bins/ \
                        -C /path/to/concoct_bins/

```
```bash
# Bin re-assembly - [MetaWRAP]
# run bin reassembly for each sample type (CT-PreInd, MC-PreInd, CT-12PostInd, MC-12PostInd)
metawrap reassemble_bins -o /path/to/output \
-1 /path/to/reads_1.fastq \
-2 /path/to/reads_2.fastq \
-t 64 -m 800 -c 50 -x 10  -b /path/to/metwrap_bins
```
```bash
# Bin quantification - [MetaWRAP]
# Use non re-assembled bins for bin_quant module
# Run bin quantification on for al bins with all read files 
metawrap quant_bins -b /path/to/all_metwrap_bins -o /path/to/output \
-a megahit_contigs.fa \
/path/to/all/reads/*.fastq
```
```bash
# Bin classification - [GTDB-tk]
# Use re-assembled bins
BINS="reassembled_bins_50_10/"
OUT="gtdbtk"
gtdbtk classify_wf --genome_dir $BINS \
--out_dir $OUT --cpus 62 --extension fa --skip_ani_screen
```
```bash
# Bin quality - [CheckM2]
checkm2 predict --threads 30 --input reassembled_bins_50_10/ -x .fa \
--output-directory /path/to/output
```

# 6. bMAG metabolic annotation and % ANI tree
```bash
Anvio Documentation: https://anvio.org/
```
```bash
# Use re-assembled bins
# Contig database - [Anvi'o]
for f in *.fa; do
            bin=$(echo ${f} | sed 's/.fa//g')
            anvi-gen-contigs-database -T 16 -f ${bin}.fa \
            -o ${bin}_contig.db
done
```
```bash
# Set-up KEGG and add kofams to contig database  - [Anvi'o]
anvi-setup-kegg-data --kegg-data-dir /path/to/output
for f in *_contig.db ;do
  bin=$(echo ${f} | sed 's/_contig.db//g')
  anvi-run-kegg-kofams -c ${bin}_contig.db --kegg-data-dir /path/to/kegg_db
done
```
```bash
# bMAG metabolic annotations - [Anvi'o]
for f in *_contig.db ;do
  bin=$(echo ${f} | sed 's/_contig.db//g')
  anvi-estimate-metabolism -c ${bin}_contig.db \
  -O ${bin} \
  --kegg-data-dir /path/to/kegg_db
done
```
```bash
# Dereplicate genomes at 97% ANI for ANI tree - [Anvi'o]
# Create external-genomes file
# get bin names
ls | sed -e 's/.fa//g' > bin_names.txt
sed -e 's/bin_names.txt//g' | sed '/^$/d' > bin_names1.txt
# get contig db paths
realpath *.db > paths.txt
# create column header
echo 'name' > name.txt
echo 'contigs_db_path' > contig_path.txt
# paste everything together
paste name.txt contig_path.txt > name_contigs_path.txt
paste bin_names.txt paths.txt > bin_names_paths.txt
cat  name_contigs_path.txt bin_names_paths.txt > external-genome.txt

# Run anvi-dereplicate
anvi-dereplicate-genomes -e external-genome.txt \
  --program fastANI \
  -o Biofilm_Induction_97ANI \
  --min-full-percent-identity 0.97 \
  --similarity-threshold 0.97
```

# 7. Identify and classify viral sequences and prophages
```bash
geNomad Documentation: https://github.com/apcamargo/genomad
```
```bash
# Viral and prophage identification and classification - [geNomad]
# concatenate all co-assemblies into one file, where the file name
# is added to the contig id 
for f in *_megahit_contigs.fa ;do
  sample=$(echo ${f} | sed 's/_megahit_contigs//g')
  sed -i 's/>/>'${sample}'/g' ${sample}_megahit_contigs.fa
done
cat *_megahit_contigs.fa > All_megahit_contigs.fa
# run geNomad to identify viral sequences
DATABASE="/path/to/genomad_db"
OUTPUT="/path/to/output"
genomad end-to-end All_megahit_contigs.fa $OUTPUT $DATABASE
```

# 8. Extend viral sequences 
```bash
COBRA Documentation: https://github.com/linxingchen/cobra
```
```bash
# Map reads to viral contigs and get coverage data - {SMALT, samtools]
# create smalt index
smalt index Megahit_contig_virusDB All_megahit_contigs_virus.fna 
# map all sample reads to viral contigs 
for f in *_R1_001_noent.fastq; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  smalt map -n 40 -y 0.95 -o ${sample}.sam \
  Megahit_contig_virusDB \
  ${sample}_R1_001_noent.fastq.fastq ${sample}_R2_001_noent.fastq
  samtools sort -o ${sample}.sort.bam ${sample}.sam
  rm ${sample}.sam
done
```
```bash
# Get coverage data and format it for COBRA - [COBRA]
for f in *.bam; do
    sample=$(echo ${f} | sed 's/.bam//g')
    jgi_summarize_bam_contig_depths --outputDepth ${sample}.sorted.depth.txt ${sample}.sort.bam
    python /path/to/cobra/coverage.transfer.py -i ${sample}.sorted.depth.txt -o ${sample}.coverage.txt
done
```
```bash
# Run COBRA to extend viral sequences - [COBRA]
for f in *.coverage.txt; do
  sample=$(echo ${f} | sed 's/.aln.sort.coverage.txt//g')
  cobra-meta -f All_megahit_contigs_050224.fa \
  -q All_megahit_contigs_virus.fna  \
  -c ${sample}.coverage.txt \
  -m ${sample}.sort.bam \
  -a megahit -mink 27 -maxk 127 \
  -o ${sample}_COBRA -t 64
done
```

# 9. Re-identify and classify extended COBRA contigs 
```bash
# Viral and prophage identification and classification - [geNomad]
DATABASE="/projectnb/viralecology/databases/genomad/genomad_db"
for f in *_contigs.new.fa ;do
        sample=$(echo ${f} | sed 's/_contigs.new.fa//g')
        genomad end-to-end ${sample}_contigs.new.fa ${sample}_GENOMAD $DATABASE
done
# Rename viral contigs that cobra identified with which sample the contig was
# extended with
for f in *_virus.fna ;do
  sample=$(echo ${f} | sed 's/_virus.fna//g')
  sed -i 's/>/>'${sample}'/g' ${sample}_virus.fna
done
# concatenate all viral sequence files together
cat *_virus.fna > All_cobra_extend_viruses.fna
```

# 10. Dereplicate viral sequences
```bash
Virathon Documentation: https://github.com/felipehcoutinho/virathon
```
```bash
# Dereplicate viral sequences - [Virathon]
python3 Virathon.py --genome_files All_cobra_extend_viruses.fna --make_pops True --threads 24
grep 'True' Seq_Info.tsv | awk '{print $1}' > Population_repIDs.txt
seqkit grep -n -f  Population_repIDs.txt All_cobra_extend_viruses.fna > Biofilm_Induction_viruses_derep.fa
```

# 11. Quality assessment of viral sequences and additional prophage identification
```bash
CheckV Documentation: https://bitbucket.org/berkeleylab/checkv/src/master/
```
```bash
# Quality assessment of viral sequences - [CheckV]
checkv end_to_end Biofilm_Induction_viruses_derep.fa Biofilm_Induction_viruses_derep_CheckV \
-d /path/to/checkv-db-v1.4 -t 16

# Get mediumd and high quality viruses
grep 'Medium-quality\|High-quality' quality_summary.tsv | awk '{print $1}' > Quality_virusIDs.txt
seqkit grep -n -f Quality_virusIDs.txt Biofilm_Induction_viruses_derep.fa > Biofilm_Induction_quality_viruses_derep.fa

# add prophages that CheckV identified that had greater than 2,000bp
seqkit grep -n -f CheckV_prophage_2000bp_IDs.txt prophage.fna > CheckV_prophage_2000bp.fa
cat Biofilm_Induction_quality_viruses_derep.fa CheckV_prophage_2000bp.fa > Biofilm_Induction_quality_viruses_prophage.fa
# add geNomad prophages that did not meat quality thresholds
seqkit grep -n -f geNomad_prophages_to_add_IDs.txt  Biofilm_Induction_viruses_derep.fa > geNomad_viruses_to_add.fa
cat  Biofilm_Induction_quality_viruses_prophage.fa  geNomad_viruses_to_add.fa > Biofilm_Induction_quality_viruses_all_prophages_derep.fa
```

# 12. Prophage induction analysis
```bash
PropagAtE Documentation: https://github.com/AnantharamanLab/PropagAtE
```
```bash
# Mapping of reads to original contigs with prophage sequences - [SMALT, samtools]
# get original contigs that prophages were identified in 
seqkit grep -n -f Prophage_original_contigIDs.txt  All_cobra_contigs.new.fa > Prophage_original_contigs.fa
# make smalt index 
smalt index Prophage_og_contigDB Prophage_original_contigs.fa
# Map reads for all samples to smalt index at 99% ID
for f in *_R1_001_noent.fastq ; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  smalt map -n 40 -y 0.99 -o ${sample}.sam \
  Prophage_og_contigDB \
  ${sample}_R1_001_noent.fastq ${sample}_R2_001_noent.fastq
  samtools sort -o ${sample}.sort.bam ${sample}.sam
  rm ${sample}.sam
done
```
```bash
# Coverage based VHR calculations - [PropagAtE]
for f in *.sort.bam; do
  name=$(echo ${f} | sed 's/.sort.bam//g')
  Propagate -f Prophage_original_contigs.fa \
  # from geNomad and CheckV 
  -v All_prophage_coordinates.tsv \
  -b ${name}.sort.bam \
  -o PropagAtE_${name} -c 1.5
done
```

# 13. Linking bMAGs to prophages
```bash
# Run geNomad on bMAGs - [geNomad]
DATABASE="/projectnb/viralecology/databases/genomad/genomad_db"
for f in *.fa ;do
        bin=$(echo ${f} | sed 's/.fa//g')
        genomad end-to-end ${bin}.fa ${bin}_GENOMAD $DATABASE
done
```
```bash
# blastn for Prophage.11 (MC-12PostInd-BF1_megahit_k127_567814_flag_0_multi_29.0000_len_13848_extended_partial|provirus_3_16519) - [blastn]
makeblastdb -in MC-12PostInd-BF1_megahit_k127_567814_provirus.fa -dbtype nucl -out k127_567814_prophageDB
for f in *.fa; do
  bin=$(echo ${f} | sed 's/.fa//g')
  blastn -db k127_567814_prophageDB \
  -query ${bin}.fa \
  -out /path/to/output \
  -outfmt "6 qseqid sseqid pident qcovs qlen slen length mismatch qstart qend sstart send evalue"
done 
```

# 14. Assign viral infection strategies
```bash
DeepPL Documentation: https://github.com/Wu-Microbiology/DeepPL
```
```bash
# get all viruses in their own fasta file
awk '/^>/{close(s); s=$0".fasta"; sub(/^>/,"",s)} {print > s}' Biofilm_Induction_quality_viruses_all_prophages_derep.fa
# make the fasta files single line
for f in *.fasta ;do
  name=$(echo ${f} | sed 's/.fasta//g')
  awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${f} > ${name}_singleline.fasta
done
# Run DeepPL
python mutilpredict_lyso_vs_lytic.py \
--model_path /path/to/deeppl_ckpt-340000 \
--fasta_folder /path/to/fasta/files \
--output_csv /path/to/output
```

# 15. Viral abundances 
```bash
# make smalt index - [smalt]
smalt index Biofilm_Induction_virusDB Biofilm_Induction_quality_viruses_all_prophages_derep.fa
# Map reads for all samples to smalt index at 95% ID - [smalt, samtools]
for f in *_R1_001_noent.fastq ; do
  sample=$(echo ${f} | sed 's/_R1_001_noent.fastq//g')
  smalt map -n 40 -y 0.95 -o ${sample}.sam \
  Biofilm_Induction_virusDB \
  ${sample}_R1_001_noent.fastq ${sample}_R2_001_noent.fastq
  samtools sort -o ${sample}.sort.bam ${sample}.sam
  rm ${sample}.sam
  samtools idxstats ${name}.sort.bam > ${name}.idxstats
done
```

# 16. Viral host prediction
```bash
iPHoP Documentation: https://bitbucket.org/srouxjgi/iphop/src/main/
```
```bash
# Predict viral-host pairs - [iPHoP]
iphop predict --fa_file Biofilm_Induction_quality_viruses_all_prophages_derep.fa  \
--db_dir /path/to/iphop_db/Aug_2023_pub_rw \
--out_dir /path/to/output
```

# 17. Viral gene annotations
```bash
MetaCerberus Documentation: https://github.com/raw-lab/MetaCerberus
```
```bash
# Annotate viral sequences using prodigal and all possible hmm provided by
# MetaCerberus- [MetaCerberus]
metacerberus.py --prodigal Biofilm_Induction_quality_viruses_all_prophages_derep.fa --hmm ALL --dir_out /path/to/output
```
