#!/bin/bash

# gvcf 整合於一個資料夾
# example:
# dpv/sample1.g.vcf.gz
# dpv/sample2.g.vcf.gz

# 製作清單
ls dpv/*.g.vcf.gz > gvcf.list

# glnexus_cli
glnexus_cli --config DeepVariantWGS --list gvcf.list \
  | bgzip -c > merge_corhort.vcf.gz
tabix -p vcf merge_corhort.vcf.gz