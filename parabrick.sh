#!/bin/bash

# === Help ===
Help() {
    echo "
使用方式: $(basename "$0") [參數]

必要參數:
  -s, --sample         Sample name

選用參數:
  -r, --ref            default: Homo_sapiens_assembly38.fasta
  -f, --folder         default: plate1
  -k, --knownSites     knownSites VCF 
  --filter             

其他:
  -h, --help   
  
  example: parallel --joblog log/fastqc.log --bar --eta -j 30 'bash parabrick.sh --folder plate2 --sample {}' :::: sample_list.txt
  note: need docker image 'docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1' first
"
  }


# === 參數解析 ===
re='^(--help|-h)$'
if [[ $1 =~ $re ]]; then
    Help
else
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -r|--ref) ref="$2"; shift 2;;
            -f|--folder) folder="$2" shift 2;;
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
    if [[ -z "$folder" ]]; then
        folder="plate1"
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


# 合併 L1~L8 的 R1 / R2
#cat /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L*_ds.*/*R1_001.fastq.gz > /home/weber/tmp/${sample}_R1_all.fastq.gz
#cat /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L*_ds.*/*R2_001.fastq.gz > /home/weber/tmp/${sample}_R2_all.fastq.gz

tmpdir="/mnt/Results/weber/tmp"
mkdir -p "$tmpdir"

# 判斷是否已合併
if [[ ! -s "/mnt/Results/weber/tmp/${sample}_R1_all.fastq.gz" ]]; then
  # 合併 R1
  (
  for LANE in {1..8}; do
      DIR=$(ls -d /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L${LANE}_ds.*)
      ls ${DIR}/*R1_001.fastq.gz | sort
  done
  ) | xargs cat > /mnt/Results/weber/tmp/${sample}_R1_all.fastq.gz
fi
if [[ ! -s "/mnt/Results/weber/tmp/${sample}_R2_all.fastq.gz" ]]; then
  # 合併 R2
  (
  for LANE in {1..8}; do
      DIR=$(ls -d /mnt/NovaSeqX/GalateaBio/plate1/fastq/${sample}_L${LANE}_ds.*)
      ls ${DIR}/*R2_001.fastq.gz | sort
  done
  ) | xargs cat > /mnt/Results/weber/tmp/${sample}_R2_all.fastq.gz
fi

qc_dir=/home/weber/QC/${sample}
mkdir -p ${qc_dir}

if [[ ! -s "${qc_dir}/${sample}_R1_all_fastqc.html" || ! -s "${qc_dir}/${sample}_R2_all_fastqc.html" ]]; then
  # 前qc
  docker run \
    --rm \
    -v /mnt:/mnt \
    -v /home:/home \
    biocontainers/fastqc:v0.11.9_cv8 \
    fastqc \
      /mnt/Results/weber/tmp/${sample}_R1_all.fastq.gz \
      /mnt/Results/weber/tmp/${sample}_R2_all.fastq.gz \
      -o ${qc_dir}/
fi  
  

# trim
out_dir=/mnt/Results/weber/fastp_trim/${sample}
mkdir -p $out_dir
if [[ ! -s "/mnt/Results/weber/fastp_trim/${sample}/${sample}_R1_trim.fastq.gz" || ! -s "/mnt/Results/weber/fastp_trim/${sample}/${sample}_R2_trim.fastq.gz" ]]; then  
  docker run --rm -v /home:/home -v /mnt:/mnt quay.io/biocontainers/fastp:0.23.2--h79da9fb_0 \
    fastp \
      -i /mnt/Results/weber/tmp/${sample}_R1_all.fastq.gz \
      -I /mnt/Results/weber/tmp/${sample}_R2_all.fastq.gz \
      -o ${out_dir}/${sample}_R1_trim.fastq.gz \
      -O ${out_dir}/${sample}_R2_trim.fastq.gz \
      -q 30 -u 50 -l 40 -w 16 -t 4 -T 4 \
      --trim_poly_g \
      -h ${out_dir}/R1_fastp.html \
      -j ${out_dir}/R2_fastp.json 
fi


if [[ ! -s "${qc_dir}/${sample}_R1_trim_fastqc.html" || ! -s "${qc_dir}/${sample}_R2_trim_fastqc.html" ]]; then
# post qc
docker run \
  --rm \
  -v /mnt:/mnt \
  -v /home:/home \
  biocontainers/fastqc:v0.11.9_cv8 \
  fastqc \
    ${out_dir}/${sample}_R1_trim.fastq.gz \
    ${out_dir}/${sample}_R2_trim.fastq.gz \
    -o ${qc_dir}/
fi

mkdir -p /home/weber/outputdir/${sample}
if [[ ! -s "/home/weber/outputdir/${sample}/${sample}_fq2bam_output.bam" ]]; then
# fq2bam 含oom應對方案
  docker run \
    --gpus all \
      --rm \
    --volume /home:/home \
    --volume /mnt:/mnt \
    nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
    pbrun fq2bam \
      --ref /home/weber/ref/${ref} \
      --in-fq ${out_dir}/${sample}_R1_trim.fastq.gz \
              ${out_dir}/${sample}_R2_trim.fastq.gz \
      --read-group-id-prefix ${sample} --read-group-sm ${sample} --read-group-lb lib1 --read-group-pl ILLUMINA \
      --out-bam /home/weber/outputdir/${sample}/${sample}_fq2bam_output.bam \
      --out-duplicate-metrics /home/weber/outputdir/${sample}/${sample}.dup_metrics.txt \
      --tmp-dir /home/weber/tmp/ \
      --low-memory
fi

if [[ -s "/home/weber/outputdir/${sample}/${sample}_fq2bam_output.bam" ]]; then
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
  
  
   # HaplotypeCaller gvcf /gvcf不能filter
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
  
  
  # Deepvariant gvcf
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
fi

if [[ -s "/home/weber/outputdir/${sample}/${sample}_filtered.haplotypecaller.vcf.gz" || ! -s "/home/weber/outputdir/${sample}/${sample}_filtered.deepvariant.vcf.gz" ]]; then
  # statistics
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
    bcftools query -f '%FILTER\n' /data/outputdir/${sample}/${sample}_filtered.haplotypecaller.vcf.gz | sort | uniq -c
  
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
    bcftools view -H /data/outputdir/${sample}/${sample}.haplotypecaller.g.vcf.gz | wc -l
  
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
    bcftools query -f '%FILTER\n' /data/outputdir/${sample}/${sample}_filtered.deepvariant.vcf.gz | sort | uniq -c
  
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
    bcftools view -H /data/outputdir/${sample}/${sample}.deepvariant.g.vcf.gz | wc -l
fi