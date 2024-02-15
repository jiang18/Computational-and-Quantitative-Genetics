## LSF resources
https://hpc.ncsu.edu/Documents/LSF.php  

## Example script for submitting LSF job

```bash
#!/bin/bash
#BSUB -n 20
#BSUB -W 120
#BSUB -R span[hosts=1]
#BSUB -R "select[model==Gold6226R]"
#BSUB -J fastGWA
#BSUB -o stdout.%J
#BSUB -e stderr.%J


~/bin/gcta64 --grm genome --make-bK-sparse 0.05 --out sp_grm

~/bin/gcta64 --bfile geno --grm-sparse sp_grm --fastGWA-mlm --pheno phen.txt --thread-num 20 --out fast
```

## Select a queue

```bash
bsub -q standard < test.sh
```

