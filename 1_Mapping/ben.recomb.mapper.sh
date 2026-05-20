#for i in X201SC19031065-Z01-F003/raw_data/P*;
#do
#newi=${i/X201SC19031065-Z01-F003\/raw_data\//}
#echo $newi $i
#fastp -i $i/*L3_1.fq.gz -I $i/*L3_2.fq.gz -o $newi\_L3_1.fq.gz -O $newi\_L3_2.fq.gz
#fastp -i $i/*L4_1.fq.gz -I $i/*L4_2.fq.gz -o $newi\_L4_1.fq.gz -O $newi\_L4_2.fq.gz
#bwa mem -R '@RG\tID:'$newi'\tSM:'$newi -t 4 genome/Ccin_Oko7.fa $newi\_L3_1.fq.gz $newi\_L3_2.fq.gz | samtools view -bh > $newi.L3.bam
#bwa mem -R '@RG\tID:'$newi'\tSM:'$newi -t 4 genome/Ccin_Oko7.fa $newi\_L4_1.fq.gz $newi\_L4_2.fq.gz | samtools view -bh > $newi.L4.bam
#samtools merge -f -t 2 all_$newi.bam $newi.L*.bam
#rm $newi.L3.bam $newi.L4.bam $newi\_L*fq.gz L4_2.fq.gz fastp*
#done

#samtools merge -@ 8 - all_P*.bam | samtools sort -@ 8 - > bams/all_samples.sorted.bam
#samtools index bams/all_samples.sorted.bam

#rm all_ P*.bam

#freebayes -p 1 -f genome/Ccin_Oko7.fa bams/all_samples.sorted.bam | /mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfkeepgeno - "GT" > vcfs/all_samples.vcf
##I originally wanted to include this, but it is giving problems, I think it is only needed for making the file smaller so not truly necessary
##| /mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfkeepinfo - "AC" "AF" "TYPE"

cat vcfs/all_samples.vcf | bgzip > vcfs/all_samples.vcf.gz
bcftools index vcfs/all_samples.vcf.gz

bcftools isec -n =2 -w 1 -c all vcfs/all_samples.vcf.gz java6/java6.variants.filtered.vcf.gz > vcfs/all_samples.goodcalls.vcf

#grep -n "#" vcfs/all_samples.goodcalls.vcf
tail -n +122 vcfs/all_samples.goodcalls.vcf | awk '{$3=$2+1;print $0}' | cut -d " " --complement -f 4-9 > vcfs/all_samples.goodcalls.bed

#cat vcfs/all_samples.goodcalls.bed | sed "s/ 0/ \-1/g" | sed "s/ \./ -/g" | sed "s/ /\t/g" > vcfs/all_samples.cleaned.bed
cat vcfs/all_samples.goodcalls.bed | sed "s/ 0/ \-1/g" | sed "s/ /\t/g" > vcfs/all_samples.cleaned.bed

head -n 122 vcfs/all_samples.goodcalls.vcf | cut --complement -f 4-9 | tail -n 1 | sed "s/#//" | sed 's/ID/END/' | sed 's/\t/,/g' > vcfs/all_samples.genotypes.csv

cut -f 1,2 genome/Ccin_Oko7.fa.fai > genome/Ccin_Oko7.genome

bedtools makewindows -g genome/Ccin_Oko7.genome -w 50000 > genome/Ccin_Oko7.50kb.windows.bed

bedtools map -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,\
51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,\
101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,\
139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,\
177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,\
215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243 -o mean -a genome/Ccin_Oko7.50kb.windows.bed -b vcfs/all_samples.cleaned.bed | sed "s/\t/,/g" >> vcfs/all_samples.genotypes.csv
