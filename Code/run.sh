#!/bin/bash
cd /research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/OriginalCodeData/

for k in {1..13};
do
for i in {1..1000};
do
  echo "${k}-${i}"
	qsub -N p"${k}-${i}" -j y -cwd -q 1-day -m abe -l h_vmem=32G -V -o permute"${k}-${i}".out -e permute"${k}-${i}".errnof -b y R CMD BATCH \"--args ${k} ${i}\" /research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Code/analysis.R permute"${k}-${i}".Rout
done
done
