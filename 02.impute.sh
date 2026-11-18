#step1 panel chunk
chr=$1
window=$2
buffer=$3

GLIMPSE2_chunk_static \
    -I chr${chr}.vcf.gz \
	  --region ${chr} \
	  --window-mb ${window} \
	  --buffer-mb ${buffer} \
	  --sequential \
	  -O ./panelbin/${chr}/chunks_chr${chr}.txt \
	  --threads 2

#step 2 panel split
chr=$1
while IFS="" read -r LINE || [ -n "$LINE" ];
do
	GLIMPSE2_split_reference_static \
		--reference chr${chr}.vcf.gz \
		--input-region ${IRG} \
		--output-region ${ORG} \
		--output ./panelbin/${chr}/chr${chr} \
		--threads 3
done < ./panelbin/${chr}/chunks_chr${chr}.txt

#step3 impute
chr=$1
while IFS="" read -r LINE || [ -n "$LINE" ];
do
	GLIMPSE2_phase_static \
		--bam-list bam.list \
		--fasta GCF_001704415.2_ARS1.2.fna \
		--reference panelbin/${chr}/chr${chr}_${CHR}_${REGS}_${REGE}.bin \
		--output  ${chr}/chr${CHR}_${REGS}_${REGE}.bcf \
		--threads 6 --Kinit 100
	ls "${chr}/chr${CHR}_${REGS}_${REGE}.bcf" >> \
		${chr}.bcflist
done < panelbin/${chr}/chunks_chr${chr}.txt

#step4 bcf_vcf

chr=$1
while IFS="" read -r LINE || [ -n "$LINE" ]; 
do
	GLIMPSE2_ligate_static \
		--input ${chr}.bcflist \
		--output ${chr}.bcf \
		--threads 4
done < panelbin/${chr}/chunks_chr${chr}.txt
	  
bcftools view ${chr}.bcf -Oz -o ${chr}.vcf.gz --threads 4
