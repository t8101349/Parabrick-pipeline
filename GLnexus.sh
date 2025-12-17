#!/bin/bash

# docker image: docker pull ghcr.io/dnanexus-rnd/glnexus:v1.4.1

# gvcf 整合於一個資料夾
# example:
# dpv/sample1.g.vcf.gz
# dpv/sample2.g.vcf.gz

# 製作清單
ls /home/weber/dpv/*.g.vcf.gz > /home/weber/gvcf.list

# glnexus_cli
docker run --rm \
  -v /home:/home \
  ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
  glnexus_cli --config DeepVariantWGS --list /home/weber/gvcf.list \
  | bcftools view -Oz -o /home/weber/outputdir/merge_cohort.vcf.gz
  
# 有空行
# 解壓移除空行
gunzip -c /home/weber/outputdir/merge_cohort.vcf.gz > /home/weber/outputdir/merge_cohort.vcf

grep -v '^$' /home/weber/outputdir/merge_cohort.vcf > /home/weber/outputdir/merge_cohort.clean.vcf

bgzip /home/weber/outputdir/merge_cohort.clean.vcf

tabix -p vcf /home/weber/outputdir/merge_cohort.clean.vcf.gz

