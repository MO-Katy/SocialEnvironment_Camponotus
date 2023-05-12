#!/bin/bash
#SBATCH --job-name=SNGCounting
#SBATCH --output=/scratch/tkay/SNG/code/%x.o.log
#SBATCH --error=/scratch/tkay/SNG/code/%x.e.log
#SBATCH --mail-user=tkay@unil.ch
#SBATCH --mail-type=FAIL,END
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3G
#SBATCH --time=23:00:00

cd /scratch/tkay/SNG/data

mkdir ./count_data

/users/tkay/subread-2.0.2-Linux-x86_64/bin/featureCounts -p -a ./genomes/evmodPasaWccl.gtf -o ./count_data/SNG_brains.txt -Q 10 -s 2 -T 20 ./aligned_transcriptomes/*.sam
