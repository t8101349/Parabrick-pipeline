wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

tar xvf parabricks_sample.tar.gz
parabricks_sample/ parabricks_sample/Data/ parabricks_sample/Data/markdup_input.bam parabricks_sample/Data/sample_2.fq.gz parabricks_sample/Data/sample_1.fq.gz parabricks_sample/Data/single_ended.bam parabricks_sample/Ref/ parabricks_sample/Ref/Homo_sapiens_assembly38.fasta parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.pac parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.ann parabricks_sample/Ref/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.amb parabricks_sample/Ref/Homo_sapiens_assembly38.dict parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.fai parabricks_sample/Ref/Homo_sapiens_assembly38.known_indels.vcf.gz parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.bwt parabricks_sample/Ref/Homo_sapiens_assembly38.fasta.sa


# 合併 L1~L8 的 R1 / R2
#cat /mnt/NovaSeqX/GalateaBio/plate1/fastq/00015988_2025-09-12_L*_ds.*/*R1_001.fastq.gz > /home/weber/tmp/00015988_2025-09-12_R1_all.fastq.gz
#cat /mnt/NovaSeqX/GalateaBio/plate1/fastq/00015988_2025-09-12_L*_ds.*/*R2_001.fastq.gz > /home/weber/tmp/00015988_2025-09-12_R2_all.fastq.gz

# 合併 R1
(
for LANE in {1..8}; do
    DIR=$(ls -d /mnt/NovaSeqX/GalateaBio/plate1/fastq/00015988_2025-09-12_L${LANE}_ds.*)
    ls ${DIR}/*R1_001.fastq.gz | sort
done
) | xargs cat > /home/weber/tmp/00015988_2025-09-12_R1_all.fastq.gz
# 合併 R2
(
for LANE in {1..8}; do
    DIR=$(ls -d /mnt/NovaSeqX/GalateaBio/plate1/fastq/00015988_2025-09-12_L${LANE}_ds.*)
    ls ${DIR}/*R2_001.fastq.gz | sort
done
) | xargs cat > /home/weber/tmp/00015988_2025-09-12_R2_all.fastq.gz


# 前qc
docker run \
  --rm \
  -v /mnt:/mnt \
  -v /home:/home \
  biocontainers/fastqc:v0.11.9_cv8 \
  fastqc \
    /home/weber/tmp/00015988_2025-09-12_R1_all.fastq.gz \
    /home/weber/tmp/00015988_2025-09-12_R2_all.fastq.gz \
    -o /home/weber/QC/test_sample/


# trim
sample=00015988_2025-09-12
out_dir=/home/weber/fastp_trim/${sample}
mkdir -p $out_dir

docker run --rm -v /home:/home quay.io/biocontainers/fastp:0.23.2--h79da9fb_0 \
  fastp \
    -i /home/weber/tmp/00015988_2025-09-12_R1_all.fastq.gz \
    -I /home/weber/tmp/00015988_2025-09-12_R2_all.fastq.gz \
    -o ${out_dir}/00015988_2025-09-12_R1_trim.fastq.gz \
    -O ${out_dir}/00015988_2025-09-12_R2_trim.fastq.gz \
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
    /home/weber/fastp_trim/00015988_2025-09-12/00015988_2025-09-12_R1_trim.fastq.gz \
    /home/weber/fastp_trim/00015988_2025-09-12/00015988_2025-09-12_R2_trim.fastq.gz \
    -o /home/weber/QC/test_sample/

"""
# fq2bam
docker run \
  --gpus all \
  --rm \
  --volume /home:/home \
  --volume /mnt:/mnt \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
pbrun fq2bam \
    --ref /home/weber/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta \
    --in-fq /home/weber/fastp_trim/00015988_2025-09-12/00015988_2025-09-12_R1_trim.fastq.gz \
            /home/weber/fastp_trim/00015988_2025-09-12/00015988_2025-09-12_R2_trim.fastq.gz \
    --read-group-id-prefix sample1 --read-group-sm sample1 --read-group-lb lib1 --read-group-pl ILLUMINA \
    --out-bam /home/weber/outputdir/fq2bam_output.bam \
    --out-duplicate-metrics /home/weber/outputdir/sample1.dup_metrics.txt 
"""

# fq2bam 含oom應對方案
docker run \
  --gpus '"device=1,2,3,4"' \
  --rm \
  --volume /home:/home \
  --volume /mnt:/mnt \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
pbrun fq2bam \
    --ref /home/weber/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta \
    --in-fq /home/weber/fastp_trim/00015988_2025-09-12/00015988_2025-09-12_R1_trim.fastq.gz \
            /home/weber/fastp_trim/00015988_2025-09-12/00015988_2025-09-12_R2_trim.fastq.gz \
    --read-group-id-prefix sample1 --read-group-sm sample1 --read-group-lb lib1 --read-group-pl ILLUMINA \
    --out-bam /home/weber/outputdir/fq2bam_output.bam \
    --out-duplicate-metrics /home/weber/outputdir/sample1.dup_metrics.txt \
    --tmp-dir /home/weber/tmp/ \
    --low-memory


# BQSR
docker run \
  --gpus '"device=4,5,6,7"' \
  --rm \
  --volume /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun bqsr  \
    --ref /home/weber/ref/Homo_sapiens_assembly38.fasta \
    --in-bam /home/weber/outputdir/fq2bam_output.bam \
    --knownSites /home/weber/known_sites/dbsnp_146.hg38.vcf.gz \
    --knownSites /home/weber/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --out-recal-file /home/weber/outputdir/test_sample.recal_data.table


# Apply BQSR
docker run \
  --gpus '"device=6,7"' \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun applybqsr \
  --ref /home/weber/ref/Homo_sapiens_assembly38.fasta \
  --in-bam /home/weber/outputdir/fq2bam_output.bam \
  --in-recal-file /home/weber/outputdir/test_sample.recal_data.table \
  --out-bam /home/weber/outputdir/test_sample.recal.bam

# HaplotypeCaller
docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun haplotypecaller \
  --ref /home/weber/ref/Homo_sapiens_assembly38.fasta \
  --in-bam /home/weber/outputdir/test_sample.recal.bam \
  --logfile /home/weber/outputdir/logs/test_sample.Haplotypecaller.log \
  --out-variants /home/weber/outputdir/test_sample_2.vcf.gz


 # gvcf /gvcf不能filter
docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun haplotypecaller \
  --ref /home/weber/ref/Homo_sapiens_assembly38.fasta \
  --in-bam /home/weber/outputdir/test_sample.recal.bam \
  --logfile /home/weber/outputdir/logs/test_sample.Haplotypecaller.g.log \
  --gvcf \
  --out-variants /home/weber/outputdir/test_sample_2.g.vcf.gz


# Deepvariant
docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun deepvariant \
  --ref /home/weber/ref/Homo_sapiens_assembly38.fasta \
  --in-bam /home/weber/outputdir/test_sample.recal.bam \
  --out-variants /home/weber/outputdir/test_sample.deepvariant_2.vcf.gz \
  --num-gpus 8 \
  --logfile /home/weber/outputdir/logs/test_sample.deepvariant.log


docker run \
  --gpus all \
  --rm \
  -v /home:/home \
  nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1 \
  pbrun deepvariant \
  --ref /home/weber/ref/Homo_sapiens_assembly38.fasta \
  --in-bam /home/weber/outputdir/test_sample.recal.bam \
  --out-variants /home/weber/outputdir/test_sample.deepvariant_2.g.vcf.gz \
  --num-gpus 8 \
  --gvcf \
  --logfile /home/weber/outputdir/logs/test_sample.deepvariant.g.log


# filter(vcf)
docker run --rm -v /home/weber:/data broadinstitute/gatk:4.3.0.0 \
  gatk VariantFiltration \
    -R /data/ref/Homo_sapiens_assembly38.fasta \
    -V /data/outputdir/test_sample_2.vcf.gz \
    -O /data/outputdir/test_sample_filtered_2.vcf.gz \
    --filter-name "HardQUAL" --filter-expression "QUAL < 3.0103" \
    --filter-name "LowGQ" --filter-expression "GQ < 0" \
    --filter-name "LowDepth" --filter-expression "DP < 1" \
    --filter-name "LowAF" --filter-expression "AF < 0.2"

docker run --rm -v /home/weber:/data broadinstitute/gatk:4.3.0.0 \
  gatk VariantFiltration \
    -R /data/ref/Homo_sapiens_assembly38.fasta \
    -V /data/outputdir/test_sample.deepvariant_2.vcf.gz \
    -O /data/outputdir/test_sample_filtered.deepvariant_2.vcf.gz \
    --filter-name "HardQUAL" --filter-expression "QUAL < 3.0103" \
    --filter-name "LowGQ" --filter-expression "GQ < 0" \
    --filter-name "LowDepth" --filter-expression "DP < 1" \
    --filter-name "LowAF" --filter-expression "AF < 0.2"

"""
#查看 dragen 格式
bcftools index -o $HOME/00015988_2025-09-12.hard-filtered.gvcf.gz.tbi /mnt/NovaSeqX/GalateaBio/plate1/gvcf/00015988_2025-09-12_ds.21417714ac12427194e62287ee5d0953/00015988_2025-09-12.hard-filtered.gvcf.gz

ln -s /mnt/NovaSeqX/GalateaBio/plate1/gvcf/00015988_2025-09-12_ds.21417714ac12427194e62287ee5d0953/00015988_2025-09-12.hard-filtered.gvcf.gz $HOME/
zcat 00015988_2025-09-12.hard-filtered.gvcf.gz | less -S

# 未trim版統計
docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools query -f '%FILTER\n' /data/outputdir/test_sample_filtered.vcf.gz | sort | uniq -c

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools view -H /data/outputdir/test_sample.g.vcf.gz | wc -l

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools query -f '%FILTER\n' /data/outputdir/test_sample_filtered.deepvariant.vcf.gz | sort | uniq -c

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools view -H /data/outputdir/test_sample.deepvariant.g.vcf.gz | wc -l

docker run --rm -u 1014:1017 -v /mnt/NovaSeqX/GalateaBio/plate1:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools query -f '%FILTER\n' /data/gvcf/00015988_2025-09-12_ds.21417714ac12427194e62287ee5d0953/00015988_2025-09-12.hard-filtered.gvcf.gz | sort | uniq -c


bcftools query -f '%FILTER\n' /mnt/NovaSeqX/GalateaBio/plate1/gvcf/00015988_2025-09-12_ds.21417714ac12427194e62287ee5d0953/00015988_2025-09-12.hard-filtered.gvcf.gz | \
awk '{filter=($0=="PASS")?"PASS":"non-PASS"; count[filter]++} END{for(f in count) print f, count[f]}'

for CHR in {1..22} X Y M; do
  echo " chr$CHR :"
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools query -r chr$CHR -f '%FILTER\n' /data/outputdir/test_sample_filtered.vcf.gz | sort | uniq -c
done

for CHR in {1..22} X Y M; do
  echo  " chr$CHR :"
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools view -H -r chr$CHR /data/outputdir/test_sample.g.vcf.gz | wc -l
done

for CHR in {1..22} X Y M; do
  echo " chr$CHR :"
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools query -r chr$CHR -f '%FILTER\n' /data/outputdir/test_sample_filtered.deepvariant.vcf.gz | sort | uniq -c
done

for CHR in {1..22} X Y M; do
  echo  " chr$CHR :"
  docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools view -H -r chr$CHR /data/outputdir/test_sample.deepvariant.g.vcf.gz | wc -l
done

for CHR in {1..22} X Y M; do
  echo " chr$CHR:"
  bcftools query -r chr$CHR -f '%FILTER\n' $HOME/00015988_2025-09-12.hard-filtered.gvcf.gz | \
  awk '{filter=($0=="PASS")?"PASS":"non-PASS"; count[filter]++} END{for(f in count) print f, count[f]}'
done
"""



# trim版 統計
docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools query -f '%FILTER\n' /data/outputdir/test_sample_filtered_2.vcf.gz | sort | uniq -c

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools view -H /data/outputdir/test_sample_2.g.vcf.gz | wc -l

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools query -f '%FILTER\n' /data/outputdir/test_sample_filtered.deepvariant_2.vcf.gz | sort | uniq -c

docker run --rm -u $(id -u):$(id -g) -v /home/weber:/data biocontainers/bcftools:v1.9-1-deb_cv1 \
bcftools view -H /data/outputdir/test_sample.deepvariant_2.g.vcf.gz | wc -l





