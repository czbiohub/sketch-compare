#!/usr/bin/env python3

import argparse
import os
import logging
import subprocess
import itertools
import gzip
import boto3
import sys
from datetime import date

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
logger = logging.getLogger(__name__)

def run(cmd, **kwargs):
    logger.info(f"Runnung command: {cmd}")
    subprocess.run(cmd, check=True, stderr=subprocess.STDOUT, **kwargs)

    
def truncate_fastqs(fastq, size):
    '''
    truncate fastq files from gzipped files
    
    '''
    new_file = fastq.replace(".fastq.gz", f"_{size}.fastq")
    logger.info(f"creating new file {new_file}")
    with gzip.open(fastq, 'rb') as f_in:
        lines = list(itertools.islice(f_in, int(size)))
        with open(new_file, 'wb') as f_out:
            for j in lines:
                f_out.write(j)
    return new_file


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("bucket_name")
    parser.add_argument("key1")
    parser.add_argument("key2")
    parser.add_argument("size")
    parser.add_argument("outbucket")
    parser.add_argument("--workdir", default=None)

    args = parser.parse_args()

    s3 = boto3.client('s3')

    try:
        
        fq_1 = os.path.basename(args.key1)
        fq_2 = os.path.basename(args.key2)

        logger.info(f"downloading from {args.bucket_name} files {args.key1}, {args.key2}")
        
        s3.download_file(args.bucket_name, args.key1, fq_1)
        s3.download_file(args.bucket_name, args.key2, fq_2)
        
        sketch_name = "_".join([fq_1.replace("_R1", ""), args.size]).replace(".fastq.gz", "")

        logger.info(f"truncating files to the first {args.size} lines")
        trunc_fq_1 = truncate_fastqs(fq_1, args.size)
        trunc_fq_2 = truncate_fastqs(fq_2, args.size)

        logger.info(f"running ska and creating sketch {sketch_name}")

        sketch_file = sketch_name + ".skf"
        sketch_command = [
            "ska",
            "fastq",
            "-o",
            sketch_name,
            trunc_fq_1,
            trunc_fq_2
        ]

        run(sketch_command)

        
        logger.info(f"uploading sketch {sketch_file} to bucket {args.outbucket}")
        s3_upload = boto3.resource('s3')
        s3_upload.meta.client.upload_file(sketch_file, args.outbucket, sketch_file)

    except Exception as e:
        logger.info(e)
        raise
