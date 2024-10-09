#!/usr/bin/env python3

import os
import json
import subprocess
import glob
import argparse
from shutil import copyfile
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import re
import shutil

def run_command(command):
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def copy_txt_files(src_dir, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    for txt_file in glob.glob(os.path.join(src_dir, "*.txt")):
        dest_file = os.path.join(dest_dir, os.path.basename(txt_file))
        copyfile(txt_file, dest_file)
        print(f"Copied {txt_file} to {dest_file}")

def run_gutsmash(config):
    input_dir = config['genome_file']
    output_base_dir = config['gutsmash_output']
    
    for file_path in glob.glob(os.path.join(input_dir, "*.fna")) + glob.glob(os.path.join(input_dir, "*.fa")) + glob.glob(os.path.join(input_dir, "*.fasta")):
        file_name = os.path.basename(file_path)
        file_base_name = 'results_' + os.path.splitext(file_name)[0]
        output_dir = os.path.join(output_base_dir, file_base_name)
        os.makedirs(output_dir, exist_ok=True)
        command = f"python {config['gutsmash_path']} --cb-knownclusters --enable-genefunctions {file_path} --output-dir {output_dir} --genefinding-tool prodigal -c {config['cpu']}"
        run_command(command)

def run_gut_extractor(config):
    command = f"python3 script/step1.py -n {config['gene_name']} -i {config['gutsmash_output']} -o {config['output_dir']}"
    run_command(command)

def translate_cds_to_protein(config, protein_dir):
    os.makedirs(protein_dir, exist_ok=True)
    for cds_file in glob.glob(f"{config['cds_dir']}/*.fna"):
        protein_file = os.path.join(protein_dir, os.path.basename(cds_file).replace(".fna", ".faa"))
        command = f"{config['seqkit_path']} translate {cds_file} -o {protein_file}"
        run_command(command)

def run_hmmsearch(config, protein_dir, hmm_output_dir):
    os.makedirs(hmm_output_dir, exist_ok=True)
    combined_output_file = os.path.join(hmm_output_dir, "HMM.txt")
    with open(combined_output_file, 'w') as outfile:
        for hmm_file in config['hmm_files']:
            outfile.write(hmm_file)
            outfile.write("\n")  

    for protein_file in glob.glob(f"{protein_dir}/*.faa"):
        for hmm_file in config['hmm_files']:
            hmm_path = os.path.join(config['hmm_dir'], hmm_file)
            output_file = os.path.join(hmm_output_dir, f"{os.path.basename(protein_file).replace('.faa', '')}_{hmm_file.replace('.hmm', '')}.tab")
            command = f"{config['hmmsearch_path']} --domtblout {output_file} --acc --noali --notextw -E {config['e_value']} --domT {config['domT']} --cpu {config['cpu']} {hmm_path} {protein_file}"
            run_command(command)

def concatenate_tab_files(hmm_output_dir, output_file):
    with open(output_file, 'w') as outfile:
        for tab_file in glob.glob(f"{hmm_output_dir}/*.tab"):
            with open(tab_file) as infile:
                outfile.write(infile.read())

def run_hmm_info_extract(config, hmm_output_dir):
    command = f"python3 script/step2_1.py HMM {hmm_output_dir} {hmm_output_dir} {hmm_output_dir} {config['score']}"
    run_command(command)

def run_hmm_gene_extract(config, hmm_tab_txt_file,species_dir):
    command = f"python3 script/step2_2.py -i {config['cds_dir']} -f {hmm_tab_txt_file} -s {species_dir} -o {config['output_dir']} -g {config['gap']} --seqkit {config['seqkit_path']}"
    run_command(command)

def main(config_file):
    with open(config_file) as f:
        config = json.load(f)
    if os.path.exists(config['gutsmash_output']):
        shutil.rmtree(config['gutsmash_output'])
        print(f"Deleted folder: {config['gutsmash_output']}")
    if os.path.exists(config['output_dir']):
        shutil.rmtree(config['output_dir'])
        print(f"Deleted folder: {config['output_dir']}")    
    
    # Step 0: Run GutSMASH
    run_gutsmash(config)
    
    # 检查 GutSMASH 输出目录是否存在
    if not os.path.exists(config['gutsmash_output']):
        print("Error: GutSMASH output directory does not exist.")
        return
    
    # Step 1: Copy gene_lab2 .txt files to ./gutsmash_output/gene_lab2/
    gene_lab2_dest_dir = os.path.join(config['gutsmash_output'], "gene_lable2")
    copy_txt_files(config['gene_lable2'], gene_lab2_dest_dir)
    species_dir = os.path.join(config['gutsmash_output'], "species")
    copy_txt_files(config['species_path'], species_dir)
    
    # Step 2: 创建 Gut-extractor 输出目录
    os.makedirs(config['output_dir'], exist_ok=True)
    run_gut_extractor(config)
    
    # Step 3: Translate CDS to protein sequences
    # 检查蛋白序列输出目录是否存在
    protein_dir=config['protein_dir']
    if not os.path.exists(protein_dir):
        print(f"Error: Protein output directory {protein_dir} does not exist.")
        return

    # Step 4: Run HMMsearch
    hmm_output_dir = os.path.join(config['output_dir'], "hmm_output")
    run_hmmsearch(config, protein_dir, hmm_output_dir)
    
    # 检查 HMMsearch 输出目录是否存在
    if not os.path.exists(hmm_output_dir):
        print(f"Error: HMMsearch output directory {hmm_output_dir} does not exist.")
        return

    # Step 5: Concatenate all .tab files into HMM.tab
    hmm_tab_file = os.path.join(hmm_output_dir, "HMM.tab")
    concatenate_tab_files(hmm_output_dir, hmm_tab_file)
    
    # 检查 HMM.tab 文件是否存在
    if not os.path.exists(hmm_tab_file):
        print(f"Error: HMM.tab file {hmm_tab_file} does not exist.")
        return

    # Step 6: Run HMM_info_extract
    run_hmm_info_extract(config, hmm_output_dir)
    
    # 检查 HMM_info_extract 输出文件是否存在
    hmm_info_output_file = os.path.join(hmm_output_dir, "HMM.tab.txt")
    if not os.path.exists(hmm_info_output_file):
        print(f"Error: HMM_info_extract output file {hmm_info_output_file} does not exist.")
        return

    # Step 7: Run HMM_gene_extract
    run_hmm_gene_extract(config, hmm_info_output_file,species_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the complete pipeline")
    parser.add_argument("config", help="Path to the configuration file")
    args = parser.parse_args()
    main(args.config)
