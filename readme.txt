Note1: 需先下載parabrick 4.6.0-1 docker image
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
Note2: fq2bam, HaplotypeCaller 步驟非常容易OOM
parabrick.sh 可能不適合平行運行
ex.bash parabrick.sh --folder plate2 --sample sample_list.txt
fastqc.sh 可以平行運行沒問題
ex. parallel --joblog log/fastqc.log --bar --eta -j 30 'bash fastqc.sh --sample {}' :::: sample_list.txt
Note3: 目前包含haplotypecaller, deepvariant(vcf/gvcf)
自行斟酌要留下哪個
