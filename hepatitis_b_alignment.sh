# Preprocess (trimming and host filtering) using nf-core taxprofiler
 nextflow run nf-core/taxprofiler -r dev -profile singularity -w ${results}/work --max_cpus 4 \
--input samplesheet.csv \
--outdir $results \
--databases database.csv \
--perform_shortread_qc \
--perform_shortread_hostremoval \
--hostremoval_reference $human_genome \
--shortread_hostremoval_index $human_index_dir_bowtie \
--save_hostremoval_unmapped


# Run alignments
for sample in ${samples[@]}; do       
       
       # Build genome index
       bowtie2-build ${genome_dir}/${species}.fasta ${genome_dir}/${species}

       # ALign taxprofiler preprocessed reads to reference genome using bowtie
       (bowtie2 -x ${genome_dir}/${species} --very-sensitive \
       -1 ${results}/bowtie2/align/${sample}_${run}.unmapped_1.fastq.gz \
       -2 ${results}/bowtie2/align/${sample}_${run}.unmapped_2.fastq.gz \
       -S ${results}/bowtie/${run}/${sample}_${run}_${species}.sam) \
       2> ${results}/bowtie/${run}/${sample}_${run}_${species}_log.txt

       # Extract aligned reads and convert to bam
       samtools view -bF 4 -h ${results}/bowtie/${run}/${sample}_${run}_${species}.sam |
       samtools sort - \
       > ${results}/bowtie/${run}/${sample}_${run}_${species}.bam

       # Remove PCR duplicated
       samtools collate -@ 4 -O -u ${results}/bowtie/${run}/${sample}_${run}_${species}.bam |
       samtools fixmate -@ 4 -m -u - - |
       samtools sort -@ 4 -u - |
       samtools markdup -@ 4 -r -d 2500 - ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup.bam \
       -f ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup_stats.txt 

       # Index bam file
       samtools index ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup.bam

       # Calculate depth
       samtools depth -a ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup.bam \
       > ${results}/bowtie/${run}/${sample}_${run}_${species}_depth.txt

       # Create mpileup file for variant calling
       samtools mpileup ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup.bam \
       --fasta-ref ${genome_dir}/${species}.fasta -a \
       > ${results}/bowtie/${run}/${sample}_${run}_${species}_mpileup.txt

       #  Create consensus sequence
       samtools consensus -f fasta -a \
       ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup.bam \
       > ${results}/bowtie/${run}/${sample}_${run}_${species}_consensus.fasta

       # Calculate genome coverage
       samtools coverage ${results}/bowtie/${run}/${sample}_${run}_${species}_dedup.bam \
       > ${results}/bowtie/${run}/${sample}_${run}_${species}_coverage.txt

done
