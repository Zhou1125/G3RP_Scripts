#clean 
sample=$1
fastp -i ${sample}_r1.fq.gz -o ${sample}_1.clean.fq.gz -I ${sample}_r2.fq.gz -O ${sample}_2.clean.fq.gz -f 3 -t 3 -l 36 -r -W 4 -M 15 --detect_adapter_for_pe --html "${sample}.html
 --json "${sample}.json --thread 10 
fastp -i ${sample}.fastq.gz -o ${sample}.clean.fq.gz -l 35 -5 -W 3 -M 15 -q 20 --html ${sample}.html --json ${sample}.json --thread 10 

#align
sample=$1
bwa mem -t 10 -R "@RG\tID:${sample}\tPL:illumina\tLB:${sample}\tSM:${sample}" GCF_001704415.2_ARS1.2.fna ${sample}_1.clean.fq.gz ${sample}_2.clean.fq.gz | samtools view -S -b | samtools sort -@ 10 -m 4G -O bam -o ${sample}.bam

#MarkDuplicates
sample=$1
gatk MarkDuplicates -I ${sample}.bam -M ${sample}.marked_dup_metrics.txt -O ${sample}_rmdup.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true

#gvcf
sample=$1
chr=$2
gatk HaplotypeCaller --emit-ref-confidence GVCF -R GCF_001704415.2_ARS1.2.fna -I ${sample}_rmdup.bam -L ${chr} -O ${sample}.${chr}.g.vcf.gz

#CombineGVCFs
chr=$1
tmp=$2
gatk GenomicsDBImport -R GCF_001704415.2_ARS1.2.fna --genomicsdb-workspace-path database_${chr} -V sample.${chr}.g.vcf.gz -L ${chr} --tmp-dir=${tmp} --reader-threads 10

gatk GenotypeGVCFs -R GCF_001704415.2_ARS1.2.fna -V gendb://database_${chr} -O Chr${chr}.vcf.gz

#SNP and InDel
chr=$1
gatk SelectVariants -R GCF_001704415.2_ARS1.2.fna -V Chr${chr}.vcf.gz -select-type SNP -O Chr${chr}.raw.SNP.vcf.gz
gatk SelectVariants -R GCF_001704415.2_ARS1.2.fna -V Chr${chr}.vcf.gz -select-type INDEL -O Chr${chr}.raw.indel.vcf.gz

#filter1
chr=$1
gatk VariantFiltration -R GCF_001704415.2_ARS1.2.fna -V Chr${chr}.raw.SNP.vcf.gz -filter "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filterName filtered -O Chr${chr}.SNP.mark.vcf.gz
gatk VariantFiltration -V Chr${chr}.raw.indel.vcf.gz -filter "QD < 2.0 || QUAL < 30.0 || FS > 200.0|| ReadPosRankSum < -20.0" --filterName filtered -O Chr${chr}.INDEL.mark.vcf.gz

#filter2
chr=$1
bcftools view  -m 2 -M 2 -v snps -f PASS --threads 6  -i  'F_MISSING<0.1' -O z -o Chr${chr}.SNP.mark.mis0.1.vcf.gz Chr${chr}.SNP.mark.vcf.gz
bcftools view  -m 2 -M 2 -V snps -f PASS --threads 6  -i  'F_MISSING<0.1' -O z -o Chr${chr}.indel.mark.mis0.1.vcf.gz Chr${chr}.indel.mark.vcf.gz

#phase
chr=$1
ne=$2
java -Xmx100g -jar beagle.22Jul22.46e.jar gt=Chr${chr}.SNP.mark.mis0.1.vcf.gz  nthreads=8 ne=${ne} out=Chr${chr}.phased
java -Xmx100g -jar beagle.22Jul22.46e.jar gt=Chr${chr}.indel.mark.mis0.1.vcf.gz  nthreads=8 ne=${ne} out=Chr${chr}.phased
