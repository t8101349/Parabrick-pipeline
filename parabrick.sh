#!/bin/bash

# === 參數解析 ===
re='^(--help|-h)$'
if [[ $1 =~ $re ]]; then
    Help
else
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -r|--ref) ref="$2"; shift 2;;
            -s|--sample) sample="$2"; shift 2;;
            -k|--knownSites) knownSites="$2"; shift 2;;
            -f|--filter) out="$2"; shift 2;;
            *) echo "unknown option: $1" >&2; exit 1;;
        esac
    done

    # === 檢查必要參數 ===
    if [[ -z "$sample" ]]; then
        echo "❌ 必要參數缺失，請確認 --sample 是否有指定。" >&2
        exit 1
    fi
    if [[ -z "$ref" ]]; then
        ref="Homo_sapiens_assembly38.fasta"
    fi
    if [[ -z "$filter" ]]; then
        filter="3.0103"
    fi
fi


#wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
#wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi

#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

#tar xvf parabricks_sample.tar.gz
#parabricks_sample/ parabricks_sample/Data/ parabricks_sample/Data/markdup_input.bam parabricks_sample/Data/sample_2.fq.gz parabricks_sample/Data/sample_1.fq.gz #parabricks_sample/Data/single_ended.bam parabricks_sample/Ref/ parabricks_sample/Ref/Homo_sapiens_assembly38.fasta parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.pac #parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.ann parabricks_sample/Ref/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi #parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.amb parabricks_sample/Ref/Homo_sapiens_assembly38.dict parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.fai #parabricks_sample/Ref/Homo_sapiens_assembly38.known_indels.vcf.gz parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.bwt #parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.sa


# 合併 L1~L8 的 R1 / R2
#cat /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L*_ds.*/*R1_001.fastq.gz > /home/weber/tmp/${sample}_R1_all.fastq.gz
#cat /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L*_ds.*/*R2_001.fastq.gz > /home/weber/tmp/${sample}_R2_all.fastq.gz

# 合併 R1
(
for LANE in {1..8}; do
    DIR=$(ls -d /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L${LANE}_ds.*)
    ls ${DIR}/*R1_001.fastq.gz | sort
done
) | xargs cat > /home/weber/tmp/${sample}_R1_all.fastq.gz
# 合併 R2
(
for LANE in {1..8}; do
    DIR=$(ls -d /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L${LANE}_ds.*)
    ls ${DIR}/*R2_001.fastq.gz | sort
done
) | xargs cat > /home/weber/tmp/${sample}_R2_all.fastq.gz

mkdir -p /home/weber/QC/${sample}
# 前qc
docker run \
  --rm \
  -v /mnt:/mnt \
  -v /home:/home \
  biocontainers/fastqc:v0.11.9_cv8 \
  fastqc \
    /home/weber/tmp/${sample}_R1_all.fastq.gz \
    /home/weber/tmp/${sample}_R2_all.fastq.gz \
    -o /home/weber/QC/${sample}/


# trim
out_dir=/home/weber/fastp_trim/${sample}
mkdir -p $out_dir

docker run --rm -v /home:/home quay.io/biocontainers/fastp:0.23.2--h79da9fb_0 \
  fastp \
    -i /home/weber/tmp/${sample}_R1_all.fastq.gz \
    -I /home/weber/tmp/${sample}_R2_all.fastq.gz \
    -o ${out_dir}/${sample}_R1_trim.fastq.gz \
    -O ${out_dir}/${sample}_R2_trim.fastq.gz \
    -q 30 -u 50 -l 40 -w 16 -t 4 -T 4 \
    --trim_poly_g \
    -h ${out_dir}/R1_fastp.html \
    -j ${out_dir}/R2_fastp.json 

# post qc
docker run \
  --rm \
  -v /mnt:/mnt \
  -v /home:/home \
  biocontainers/fastqc:v0.11.9_cv8 \
  fastqc \
    /home/weber/fastp_trim/${sample}/${sample}_R1_trim.fastq.gz \
    /home/weber/fastp_trim/${sample}/${sample}_R2_trim.fastq.gz \
    -o /home/weber/QC/${sample}/


mkdir -p /home/weber/outputdir/${sample}
# fq2bam 含oom應對方案
docker run \
  --gpus '"device=1,2,3,4"' \
  --rm \
  --volume /home:/home \
  --volume /mnt:/mnt \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
pbrun fq2bam \
    --ref /home/weber/parabricks_sample/Ref/${ref} \
    --in-fq /home/weber/fastp_trim/${sample}/${sample}_R1_trim.fastq.gz \
            /home/weber/fastp_trim/${sample}/${sample}_R2_trim.fastq.gz \
    --read-group-id-prefix ${sample} --read-group-sm ${sample} --read-group-lb lib1 --read-group-pl ILLUMINA \
    --out-bam /home/weber/outputdir/${sample}/${sample}_fq2bam_output.bam \
    --out-duplicate-metrics /home/weber/outputdir/${sample}/${sample}.dup_metrics.txt \
    --tmp-dir /home/weber/tmp/ \
    --low-memory


# BQSR
docker run \
  --gpus '"device=4,5,6,7"' \
  --rm \
  --volume /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun bqsr  \
    --ref /home/weber/ref/${ref} \
    --in-bam /home/weber/outputdir/${sample}/${sample}_fq2bam_output.bam \
    --knownSites /home/weber/known_sites/dbsnp_146.hg38.vcf.gz \
    --knownSites /home/weber/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --out-recal-file /home/weber/outputdir/${sample}/${sample}.recal_data.table


# Apply BQSR
docker run \
  --gpus '"device=6,7"' \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun applybqsr \
  --ref /home/weber/ref/${ref} \
  --in-bam /home/weber/outputdir/${sample}/${sample}_fq2bam_output.bam \
  --in-recal-file /home/weber/outputdir/${sample}/${sample}.recal_data.table \
  --out-bam /home/weber/outputdir/${sample}/${sample}.recal.bam

# HaplotypeCaller
docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun haplotypecaller \
  --ref /home/weber/ref/${ref} \
  --in-bam /home/weber/outputdir/${sample}/${sample}.recal.bam \
  --logfile /home/weber/outputdir/logs/${sample}.Haplotypecaller.log \
  --out-variants /home/weber/outputdir/${sample}/${sample}.haplotypecaller.vcf.gz


 # gvcf /gvcf不能filter
docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun haplotypecaller \
  --ref /home/weber/ref/${ref} \
  --in-bam /home/weber/outputdir/${sample}/${sample}.recal.bam \
  --logfile /home/weber/outputdir/logs/${sample}.Haplotypecaller.g.log \
  --gvcf \
  --out-variants /home/weber/outputdir/${sample}/${sample}.haplotypecaller.g.vcf.gz


# Deepvariant
docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun deepvariant \
  --ref /home/weber/ref/${ref} \
  --in-bam /home/weber/outputdir/${sample}/${sample}.recal.bam \
  --out-variants /home/weber/outputdir/${sample}/${sample}.deepvariant.vcf.gz \
  --num-gpus 8 \
  --logfile /home/weber/outputdir/logs/${sample}.deepvariant.log


docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun deepvariant \
  --ref /home/weber/ref/${ref} \
  --in-bam /home/weber/outputdir/${sample}/${sample}.recal.bam \
  --out-variants /home/weber/outputdir/${sample}/${sample}.deepvariant.g.vcf.gz \
  --num-gpus 8 \
  --gvcf \
  --logfile /home/weber/outputdir/logs/${sample}.deepvariant.g.log


# filter(vcf)
docker run --rm -v /home/weber:/data broadinstitute/gatk:4.3.0.0 \
  gatk VariantFiltration \
    -R /data/ref/${ref} \
    -V /data/outputdir/${sample}/${sample}.haplotypecaller.vcf.gz \
    -O /data/outputdir/${sample}/${sample}_filtered.haplotypecaller.vcf.gz \
    --filter-name "HardQUAL" --filter-expression "QUAL < 3.0103" \
    --filter-name "LowGQ" --filter-expression "GQ < 0" \
    --filter-name "LowDepth" --filter-expression "DP < 1" \
    --filter-name "LowAF" --filter-expression "AF < 0.2"

docker run --rm -v /home/weber:/data broadinstitute/gatk:4.3.0.0 \
  gatk VariantFiltration \
    -R /data/ref/${ref} \
    -V /data/outputdir/${sample}/${sample}.deepvariant.vcf.gz \
    -O /data/outputdir/${sample}/${sample}_filtered.deepvariant.vcf.gz \
    --filter-name "HardQUAL" --filter-expression "QUAL < 3.0103" \
    --filter-name "LowGQ" --filter-expression "GQ < 0" \
    --filter-name "LowDepth" --filter-expression "DP < 1" \
    --filter-name "LowAF" --filter-expression "AF < 0.2"



# trim版 統計
docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools query -f '%FILTER\n' /data/outputdir/${sample}/${sample}_filtered.haplotypecaller.vcf.gz | sort | uniq -c

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools view -H /data/outputdir/${sample}/${sample}.haplotypecaller.g.vcf.gz | wc -l

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools query -f '%FILTER\n' /data/outputdir/${sample}/${sample}_filtered.deepvariant.vcf.gz | sort | uniq -c

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools view -H /data/outputdir/${sample}/${sample}.deepvariant.g.vcf.gz | wc -l


