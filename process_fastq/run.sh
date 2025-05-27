python spatial.py \
    --r1 data/SC20250402/SC25040207_L5_1.fq.gz \
    --r2 data/SC20250402/SC25040207_L5_2.fq.gz \
    --sample SC25040207 \
    --outdir ~/projects/spatial_ct/SC20250402/ \
    --image data/SC20250402/SC20250402_registedimage.tif \
    --alignment data/SC20250402/SC20250402_result.json \
    --barcode 0_8_38_46_101 --cores 32 \
    --cratac ~/pipelines/cellranger-atac-2.1.0/bin/cellranger-atac \
    --refdata ~/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0

