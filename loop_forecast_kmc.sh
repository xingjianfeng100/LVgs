#!/bin/bash
set -e -o pipefail

mkdir -p  intermediate
echo -e "Kmer_len\tGS\tRepeat\tUnique\tModel_fit\tError\tKmer_type\tKmer_total\tKcov" > ${outputpref}.stat.tab
for i in `seq ${bg_len} ${step} ${ed_len}`
do
	mkdir -p intermediate/k${i}
	cd intermediate/k${i}
	mkdir -p fas_tmp
	kmc -k${i} -t${threads} -m${mem} ${kmcoption} ${reads} kmer fas_tmp/ 
	kmc_tools transform kmer histogram k${i}.histo -cx1000000
	genomescope.R -i k${i}.histo -o gsco_k${i} -k ${i}  -p ${pl} ${genomescopeoption}
	rm -r fas_tmp
	gs=`grep "Genome Haploid Length" gsco_k${i}/summary.txt | awk 'FS="\t" {print $6}' | sed 's/,//g'`
	repeat=`grep "Genome Repeat Length" gsco_k${i}/summary.txt | awk 'FS="\t" {print $6}' | sed 's/,//g'`
	Unique=`grep "Genome Unique Length" gsco_k${i}/summary.txt | awk 'FS="\t" {print $6}' | sed 's/,//g'`
	Model=`grep "Model Fit" gsco_k${i}/summary.txt | awk 'FS="\t" {print $4}' | sed 's/%//g'`
	error=`grep "Read Error Rate" gsco_k${i}/summary.txt | awk 'FS="\t" {print $5}' | sed 's/%//g'`
	kmer_type=`awk '{ sum += $2 } END { print sum }' k${i}.his`
	kmer_total=`awk '{ sum += $1*$2 } END { print sum }' k${i}.his`
	kcov=`grep "kmercov " gsco_k${i}/model.txt | awk '{printf"%.2f\n", $2}'`
	cd ../../
	echo -e "$i\t$gs\t$repeat\t$Unique\t$Model\t$error\t$kmer_type\t$kmer_total\t$kcov" >> ${outputpref}.stat.tab	
done

