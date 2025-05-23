#!/bin/bash
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=spct
#SBATCH -A zhanghuiminggroup
#SBATCH --qos normal
#SBATCH --mem=128G
#SBATCH --export=TENX_IGNORE_DEPRECATED_OS=1

source ~/miniforge3/bin/activate
source activate ~/miniforge3/envs/SUPER-Tag
ulimit -u 8096

python spatial.py \
    --r1 data/SC20250402XFFPE/SC25040207_L5_1.fq.gz \
    --r2 data/SC20250402XFFPE/SC25040207_L5_2.fq.gz \
    --sample SC25040207 \
    --outdir ~/projects/spatial_ct/SC25040207_test/ \
    --image data/SC20250402XFFPE/SC25032707_registedimage.tif \
    --alignment data/SC20250402XFFPE/SC25032707_result.json \
    --barcode 0_8_38_46_101 --cores 32 \
    --cratac ~/pipelines/cellranger-atac-2.1.0/bin/cellranger-atac \
    --refdata ~/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0

