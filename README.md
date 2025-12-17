# get parabrick docker first
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
# recommend
Variant calling: DeepVariant \
Joint genotyping: GLnexus \

Note1: GLnexus docker image \
docker pull ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
Note2: fq2bam, HaplotypeCaller 步驟非常容易OOM \
parabrick.sh 可能不適合平行運行 \
ex.bash parabrick.sh --folder plate1 --sample 00021787_2025-09-12 \
bash parabrick.sh --folder plate2 --sample BY31700_A02Wgenom \
fastqc.sh 可以平行運行沒問題(只做完fastqc部分) \
ex. parallel --joblog log/fastqc.log --bar --eta -j 30 'bash fastqc.sh --folder plate2 --sample {}' :::: sample_list.txt \
Note3: 目前包含haplotypecaller, deepvariant(vcf/gvcf) \
自行斟酌要留下哪個
