from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen
import os
import argparse
from pathlib import Path
import os
import subprocess as sb
import cv2
import json
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import argparse

def split_fastq(input_file_R1, output_file_R1, output_file_R2, barcode=[0,8,38,46,101]):
    bc1_start, bc1_end, bc2_start, bc2_end, seq_start = barcode
    with gzopen(input_file_R1, "rt") as in_handle_R1, open(output_file_R1.rstrip('.gz'), "w") as out_handle_R1, open(output_file_R2.rstrip('.gz'), "w") as out_handle_R2:
        for title, seq, qual in FastqGeneralIterator(in_handle_R1):
            new_seq_R1 = seq[seq_start:]
            new_qual_R1 = qual[seq_start:]
            barcode = seq[bc1_start:bc1_end] + seq[bc2_start:bc2_end]
            new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]        
            out_handle_R1.write("@%s\n%s\n+\n%s\n" % (title, new_seq_R1, new_qual_R1))
            out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual_R2))
    os.system(f"pigz -f {output_file_R1.rstrip('.gz')}")
    os.system(f"pigz -f {output_file_R2.rstrip('.gz')}")

def spatial_image(img, alignment, out):
    Path(f"{out}/spatial").mkdir(exist_ok=True,parents=True)
    img = cv2.imread(img)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    elems = json.load(open(alignment,encoding="utf-8"))
    try:
        df = pd.DataFrame(elems['spots'])
        points = np.array(elems['points'],np.int32)

        aligned_fiducials = cv2.polylines(img, [points], isClosed=True, color=(255,0,0), thickness=5)
    except:
        df = pd.DataFrame(elems)
        aligned_fiducials = False
        points = False
    spot_diameter = float(df['diameter'].values[0])
    df["pxl_col_in_fullres"] = df['x'].astype(float)
    df["pxl_row_in_fullres"] = df['y'].astype(float)
    df['in_tissue']          = df['tissue'].map({"False":0,"True":1})
    df['x']                  = df['col']
    df['y']                  = df['row']
    df = df.set_index("barcode")
    df = df.sort_values(by=["y",'x'])

    df = df[['in_tissue','y','x','pxl_row_in_fullres','pxl_col_in_fullres']]
    df.index = [i+'-1' for i in df.index]
    df.to_csv(f"{out}/spatial/tissue_positions_list.csv",header=None)
    y,x = img.shape[:2]
    if y >= x:
        d={"tissue_hires_scalef":2000/y,
        "tissue_lowres_scalef":600/y,
        "fiducial_diameter_fullres":spot_diameter+30,
        "spot_diameter_fullres":spot_diameter}
    else:
        d={"tissue_hires_scalef":2000/x,
        "tissue_lowres_scalef":600/x,
        "fiducial_diameter_fullres":spot_diameter+30,
        "spot_diameter_fullres":spot_diameter}

    hires_img = cv2.resize(img,(0,0),fx=d['tissue_hires_scalef'],fy=d['tissue_hires_scalef'])
    lowres_img = cv2.resize(img,(0,0),fx=d['tissue_lowres_scalef'],fy=d['tissue_lowres_scalef'])
    with open(f"{out}/spatial/scalefactors_json.json","w") as outf:
        json.dump(d,outf)

    cv2.imwrite(f"{out}/spatial/tissue_hires_image.png", hires_img)
    cv2.imwrite(f"{out}/spatial/tissue_lowres_image.png", lowres_img)

def spatial_function(R1, R2, sample_name, outdir, image, alignment, cratac, refdata, barcode = [0,8,38,46,101], cores = 24):
    O1 = f'{outdir}/fastq_cr/{sample_name}_S1_L001_R1_001.fastq.gz'
    I1 = f'{outdir}/fastq_cr/{sample_name}_S1_L001_R2_001.fastq.gz'
    O2 = f'{outdir}/fastq_cr/{sample_name}_S1_L001_R3_001.fastq.gz'
    Path(f'{outdir}/fastq_cr/').mkdir(exist_ok = True, parents =True)
    print(f"processing {sample_name}...")
    
    def run_cr():
        print("running cellrang atac...")
        command_str = f"{cratac} count --id={sample_name} \
        --reference={refdata} \
        --fastqs={outdir}/fastq_cr/ \
        --sample={sample_name} --localcores={cores} --localmem=64 --force-cells 2500 "
        sb.call(command_str, shell=True)

    os.chdir(outdir)

    split_fastq(R1, O1.rstrip('.gz'), I1.rstrip('.gz'), barcode)
    if not os.path.exists(O2):
        os.system(f"cp {R2} {O2}")

    run_cr()
    spatial_image(image, alignment, f"{outdir}/")


if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--r1', required=True, help="Read1")
    ap.add_argument('--r2', required=True, help="Read2")
    ap.add_argument('-s', '--sample', required=True, help="Read2")
    ap.add_argument('-o', '--outdir', required=True, help="Read2")
    ap.add_argument('--image', required=True, help="Image file for spatial function")
    ap.add_argument('--alignment', required=True, help="JSON file for spatial function")
    ap.add_argument('--barcode', required=False, help="barcode loci", default = '0_8_38_46_95')
    ap.add_argument('--cores', required=False, help="CPU cores to use", default = 24)
    ap.add_argument('--cratac', required=True, help="path to cratac")
    ap.add_argument('--refdata', required=True, help="path to 10x reference")

    args = vars(ap.parse_args())
    barcode = [int(i) for i in args['barcode'].split('_')]
    spatial_function(args['r1'], args['r2'], args['sample'], args['outdir'], args['image'], args['alignment'], args['cratac'], args['refdata'], barcode, args['cores'])