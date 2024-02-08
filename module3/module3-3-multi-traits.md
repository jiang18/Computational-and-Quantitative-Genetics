## Software tool
PLINK
[MPH](https://jiang18.github.io/mph/)

## Data set
The same as in Modules 1 and 2.

## Genetic correlations

```
echo "genome 10" > A.grm.txt

mkdir multi-trait

# Use a single GRM in multi-trait analysis.
mph --minque --grm_list A.grm.txt --phenotype phen.csv --trait milk,fat,fat_percent --num_threads 10 --out ./multi-trait/genome

# Genetic correlation estimates are in ./multi-trait/genome.mq.cor.csv.

# Use multiple GRMs in multi-trait analysis.
mph --minque --save_mem --grm_list chr.grms.txt --phenotype phen.csv --trait milk,fat,fat_percent --num_threads 10 --out ./multi-trait/chromosomes

# Genetic correlation estimates are in ./multi-trait/chromosomes.mq.cor.csv.

```

