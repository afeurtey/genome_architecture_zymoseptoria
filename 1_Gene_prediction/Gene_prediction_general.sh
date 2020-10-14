#!/bin/bash
#
#  example submit script for a serial job
#  submit by  sbatch serial-job.sh
#
#  specify the job name
#SBATCH --job-name=annZt05
#  how many cpus are requested
#SBATCH --ntasks=10
#  run on one node, important if you have more than 1 ntasks
#SBATCH --nodes=1
#  maximum walltime
#SBATCH --time=100:00:00
#  maximum requested memory
#SBATCH --mem=100G
#  write std out and std error to these files
#SBATCH --error=annZt05.sterr
#SBATCH --output=annZt05.stout
#  send a mail for job start, end, fail, etc.
#SBATCH --mail-type=NONE
#SBATCH --mail-user=feurtey@evolbio.mpg.de
#  which partition?
#  there are global,testing,highmem,standard,fast
#SBATCH --partition=global


sample="Zt05"

# Directories
# -----------
general_files=/home/feurtey/Assemblies_and_annotations/2_Annotations/General_files/
work_dir=/home/feurtey/Assemblies_and_annotations/2_Annotations/${sample}/
raw_data_dir="0-Needed_data/"
filtered_dir="1-Trimmed_reads/"
aligned_dir="2-Alignments/"
predict_dir="3-Gene_prediction/"
abinitio_dir="0-Ab_initio/"
mapping_dir="1-Mapping_based/"
denovo_dir="2-Transcript_assembly/"
consensus_dir="3-Weighted_consensus/"

source_genome_dir=/groups/mpistaff/Stukenbrock_lab/Assemblies_and_annotations/Pacbio_Assemblies_2018/
unmasked_genome_dir="2_Final_unmasked_assemblies/"
masked_genome_dir="3_Final_repeat_softmasked_assemblies/"



# Software
# --------
PASAHOME="/data/biosoftware/PASA/PASApipeline-v2.3.3/"
TRINITYPATH="/data/biosoftware/trinity/trinityrnaseq-Trinity-v2.8.4/"
AUGUSTUSPATH="/data/biosoftware/augustus/augustus-3.3/config/"
GENEMARKPATH="/data/biosoftware/GeneMark/gm/"
EVMPATH="/data/biosoftware/EVidenceModeler/EVidenceModeler-1.1.1/"
MAKERPATH="/data/biosoftware/maker/maker3/maker/bin/"
BAMTOOLSPATH="/data/biosoftware/bamtools/bamtools/"
SAMTOOLSPATH="/data/biosoftware/samtools/samtools-1.4/samtools"


# Files
# ------
adapters_file=${general_files}adapters.txt
fobamfiles=${work_dir}${predict_dir}${denovo_dir}Fastq_files_description.txt
PASA_contig_file=${PASAHOME}sample_data/sqlite.confs/alignAssembly.config
EVM_weight_files=${work_dir}${predict_dir}${consensus_dir}weights.txt
header_file=${general_files}header_and_necessary_info.sh

source_unmasked_genome_file=${source_genome_dir}${unmasked_genome_dir}${sample}_pacbio_assembly_2019_filtered_on_cov_2X_unmasked.fa
unmasked_genome_file=${work_dir}${raw_data_dir}${sample}_pacbio_assembly_2019_filtered_on_coverage_2X_unmasked.fa
source_masked_genome_file=${source_genome_dir}${masked_genome_dir}${sample}_assembly_2019_final_repeats_softmasked.fa
masked_genome_file=${work_dir}${raw_data_dir}${sample}_assembly_2019_final_repeats_softmasked.fa

gff_all_pred=${work_dir}${predict_dir}${consensus_dir}${sample}_all_predictions.gff
gff_transcripts=${work_dir}${predict_dir}${consensus_dir}${sample}_all_transcripts.gff





#    <<>><<>><<>><<>><<>><<>><<>><<>><<>>
# |   Directories, read and genome prep   |
#    <<>><<>><<>><<>><<>><<>><<>><<>><<>>

cd ${work_dir}

# To convert from my old file name to the new one
#if [ -d "0-Untrimmed_reads" ]; then
#mv 0-Untrimmed_reads ${raw_data_dir%/}
#fi

# Create directories
#mkdir ${work_dir}${raw_data_dir}
#mkdir ${work_dir}${filtered_dir}
#mkdir ${work_dir}${predict_dir}
#mkdir ${work_dir}${predict_dir}${abinitio_dir}
#mkdir ${work_dir}${predict_dir}${mapping_dir}
#mkdir ${work_dir}${predict_dir}${denovo_dir}
#mkdir ${work_dir}${predict_dir}${consensus_dir}
#mkdir ${work_dir}${aligned_dir}By_strand/

# Get reference genomes
#cp ${source_unmasked_genome_file} ${unmasked_genome_file}
#cp ${source_masked_genome_file} ${masked_genome_file}


# Trimming reads
#cd ${work_dir}${raw_data_dir}
#for fq_file in *fastq* ; do
#    echo ${fq_file}
#    echo ${work_dir}${raw_data_dir}${fq_file%fastq*}adapt.qual30.minlen30.fastq.gz
#    fastqc ${fq_file}
#    java -jar /data/biosoftware/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar \
#     SE -threads 10 \
#     ${fq_file} \
#     ${work_dir}${filtered_dir}${fq_file%fastq*}adapt.qual30.minlen30.fastq.gz \
#     ILLUMINACLIP:${adapters_file}:2:30:10 \
#     LEADING:30 TRAILING:30 MINLEN:30

#    fastqc ${work_dir}${filtered_dir}${fq_file%fastq*}adapt.qual30.minlen30.fastq.gz
#done



#    <<>><<>><<>><<>>
# |  Read alignments  |
#    <<>><<>><<>><<>>


#cd ${work_dir}${filtered_dir}
#hisat2-build ${unmasked_genome_file} ${work_dir}${aligned_dir}${sample}

#for fq_file in *adapt.qual30.minlen30.fastq.gz; do
#    aln_file=${fq_file%.adapt.qual30.minlen30.fastq.gz}
#    hisat2 \
#     -x ${work_dir}${aligned_dir}${sample} \
#     -U $fq_file \
#     -S ${work_dir}${aligned_dir}${aln_file}.sam \
#     --un ${work_dir}${aligned_dir}${aln_file}_unaligned_reads.fastq \
#     --no-unal --met-stderr --phred33 \
#     --rna-strandness R \
#     --min-intronlen 20 --max-intronlen 15000 \
#     --summary-file ${work_dir}${aligned_dir}${aln_file}_summary.txt \
#     --threads 10

#    samtools flagstat ${work_dir}${aligned_dir}${aln_file}.sam \
#      > ${work_dir}${aligned_dir}${aln_file}_flagstat.txt
#    samtools view -b ${work_dir}${aligned_dir}${aln_file}.sam > ${work_dir}${aligned_dir}${aln_file}.bam

    # Separate the reads by strand
#    samtools view  -F 2308 -f 16 -b ${work_dir}${aligned_dir}${aln_file}.bam -o ${work_dir}${aligned_dir}By_strand/${aln_file}.reverse.bam
#    samtools view  -F 2324 -b ${work_dir}${aligned_dir}${aln_file}.bam -o ${work_dir}${aligned_dir}By_strand/${aln_file}.forward.bam

#    samtools sort ${work_dir}${aligned_dir}By_strand/${aln_file}.forward.bam --output-fmt BAM -o ${work_dir}${aligned_dir}By_strand/${aln_file}.forward.sorted.bam
#    samtools index ${work_dir}${aligned_dir}By_strand/${aln_file}.forward.sorted.bam

#    samtools sort ${work_dir}${aligned_dir}By_strand/${aln_file}.reverse.bam --output-fmt BAM -o ${work_dir}${aligned_dir}By_strand/${aln_file}.reverse.sorted.bam
#    samtools index ${work_dir}${aligned_dir}By_strand/${aln_file}.reverse.sorted.bam

#done




#   <<>><<>><<>><<>><<>>
# |   Gene prediction   |
#   <<>><<>><<>><<>><<>>

#   -------------
# |   Ab initio  |
#   -------------

#echo "Starting ab initio gene prediction using genemark-ES"
Download_reads_loop_rerun_02_prefetch.sh
#cd ${work_dir}${predict_dir}${abinitio_dir}
#module load perl/5.24.1

#${GENEMARKPATH}gmes_petap.pl --ES  --fungus \
# --sequence ${masked_genome_file} \
# --cores 10

#${MAKERPATH}genemark_gtf2gff3 ${work_dir}${predict_dir}${abinitio_dir}genemark.gtf \
# > ${work_dir}${predict_dir}${abinitio_dir}${sample}_ab_initio.gff

#echo "Finished ab initio gene prediction using genemark-ES"



#   ----------------
# |   Mapping based  |
#   ----------------

echo "Starting mapping based gene prediction using braker1"

#ls -1 ${work_dir}${aligned_dir}*bam > ${work_dir}${predict_dir}${mapping_dir}List_of_bam_files.txt

#samtools merge -b ${work_dir}${predict_dir}${mapping_dir}List_of_bam_files.txt ${work_dir}${predict_dir}${mapping_dir}${sample}_merged.bam

#samtools sort -n ${work_dir}${predict_dir}${mapping_dir}${sample}_merged.bam > ${work_dir}${predict_dir}${mapping_dir}${sample}_merged.sorted.bam

export AUGUSTUS_CONFIG_PATH=${AUGUSTUSPATH}

module load perl/5.24.1
perl /data/biosoftware/braker1/braker1/braker.pl \
 --genome=${masked_genome_file} \
 --bam=${work_dir}${predict_dir}${mapping_dir}${sample}_merged.sorted.bam \
 --AUGUSTUS_CONFIG_PATH=${AUGUSTUSPATH} \
 --GENEMARK_PATH=${GENEMARKPATH} \
 --BAMTOOLS_PATH=${BAMTOOLSPATH} \
 --SAMTOOLS_PATH=${SAMTOOLSPATH} \
 --workingdir=${work_dir}${predict_dir}${mapping_dir} \
 --species=${sample} \
 --fungus \
 --cores 10


${MAKERPATH}genemark_gtf2gff3 ${work_dir}${predict_dir}${mapping_dir}braker/${sample}/augustus.gtf \
 > ${work_dir}${predict_dir}${mapping_dir}${sample}_mapping_based.gff

echo "Finished mapping based gene prediction using braker1"


#   ----------------
# |   Trinity/PASA   |
#   ----------------

# Genome guided trinity assemblies
echo "Starting genome-guided transcripts assembly using Trinity"

#samtools sort ${work_dir}${predict_dir}${mapping_dir}${sample}_merged.bam > ${work_dir}${predict_dir}${mapping_dir}${sample}_merged.sorted2.bam

#module load perl/5.26.1
#module load python/python-anaconda-5.0.1
#${TRINITYPATH}Trinity \
# --genome_guided_bam ${work_dir}${predict_dir}${mapping_dir}${sample}_merged.sorted2.bam \
# --SS_lib_type R \
# --max_memory 100G --CPU 10 \
# --genome_guided_max_intron 10000 \
# --output ${work_dir}${predict_dir}${denovo_dir}trinity_outputs

#gzip ${work_dir}${filtered_dir}*fastq
assembled_transcripts=${work_dir}${predict_dir}${denovo_dir}trinity_outputs/Trinity-GG.fasta
#echo "Finished genome-guided transcripts assembly using	Trinity"

# Alignment with PASA
echo "Starting alignment with PASA"

module load perl/5.26.1
#cd ${work_dir}${predict_dir}${denovo_dir}trinity_outputs/

#${PASAHOME}bin/seqclean  ${assembled_transcripts} -c 10
#${PASAHOME}Launch_PASA_pipeline.pl \
#           --config ${PASA_contig_file} \Download_reads_loop_rerun_02_prefetch.sh
#           --create \
#           --run \
#           --genome ${unmasked_genome_file} \
#           --transcripts ${assembled_transcripts}.clean \
#           -T -u ${assembled_transcripts} \
#           --ALIGNERS blat,gmap --CPU 5 \
#           --stringent_alignment_overlap 30.0 \
#           --transcribed_is_aligned_orient
#echo "Finished alignment with PASA"


# Preparing weights file
echo "ABINITIO_PREDICTION   GeneMarkES   2" > ${EVM_weight_files}
echo "ABINITIO_PREDICTION   BRAKER1   5" >> ${EVM_weight_files}
echo "TRANSCRIPT      PASA_transcript_assemblies      10" >> ${EVM_weight_files}

#Preparing gff file

echo "##gff-version 3" > ${gff_all_pred}
grep -v "#" ${work_dir}${predict_dir}${abinitio_dir}${sample}_ab_initio.gff \
   | awk '  BEGIN {OFS = "\t"}  {{$2 = "GeneMarkES"}; print } ' > ${gff_all_pred}
grep -v "#" ${work_dir}${predict_dir}${mapping_dir}${sample}_mapping_based.gff \
   | awk '  BEGIN {OFS = "\t"}  {{$2 = "BRAKER1"}; print } ' >> ${gff_all_pred}

cat ${work_dir}${predict_dir}${denovo_dir}trinity_outputs/sample_mydb_pasa.sqlite.pasa_assemblies.gff3 \
   | awk '  BEGIN {OFS = "\t"}  {{$2 = "PASA_transcript_assemblies"}; print } '  > ${gff_transcripts}




#Running the consensus analysis

cd ${work_dir}${predict_dir}${consensus_dir}
${EVMPATH}EvmUtils/partition_EVM_inputs.pl \
 --genome ${masked_genome_file} \
 --gene_predictions ${gff_all_pred} \
 --transcript_alignments ${gff_transcripts} \
 --segmentSize 100000 --overlapSize 10000 \
 --partition_listing ${work_dir}${predict_dir}${consensus_dir}${sample}_partitions_list.out

${EVMPATH}EvmUtils/write_EVM_commands.pl \
 --genome ${masked_genome_file} \
 --weights ${EVM_weight_files} \
 --gene_predictions ${gff_all_pred} \
 --transcript_alignments ${gff_transcripts} \
 --output_file_name ${sample}_evm.out \
 --partitions ${work_dir}${predict_dir}${consensus_dir}${sample}_partitions_list.out \
 >  ${work_dir}${predict_dir}${consensus_dir}${sample}_commands.list

contigs="$(grep ">" ${unmasked_genome_file} | sed 's/>//' )"

for i in $contigs;
do
echo ${i}
cat $header_file > EVM_$i.sh
sed -i "s/TK_to_change/EVM_${i}/" EVM_$i.sh
sed -i "s/TK_sample_name/${sample}/" EVM_$i.sh
grep -w "$i" ${work_dir}${predict_dir}${consensus_dir}${sample}_commands.list >> EVM_$i.sh
sbatch EVM_$i.sh
done
