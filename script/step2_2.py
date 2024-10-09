#!/bin/python3

import os
import re
import sys
import glob
import argparse
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import subprocess

def getSeq(fafile, genelist):
    rank = 0
    out_seq = ''
    seq_id = ''
    fa_file = open(fafile)
    for line in fa_file:
        if line.startswith('>') and rank == 0:
            seq_id = line.strip()
            if seq_id.split(' ')[0][1:] == genelist or (re.search('protein_id=.*?\]', seq_id) != None and re.search('protein_id=(.*?)\]', seq_id).group(1) == genelist):
                rank = 1
        elif line.startswith('>') and rank == 1:
            break
        elif rank == 1:
            out_seq += line.strip()
        else:
            continue
    return seq_id, out_seq

def speciesExtract(species, species_info):
    print(f"Processing species file: {species}")
    try:
        try:
            species_file = pd.read_table(species, usecols=[0, 7, 8], header=None, engine='python')
        except:
            species_file = pd.read_table(species, encoding='unicode_escape', usecols=[0, 7, 8], header=None, engine='python')
    except Exception as e:
        print(f"Error reading species file {species}: {e}")
        return species_info
    for i in range(species_file.shape[0]):
        if str(species_file.iloc[i, 2]) == 'nan':
            species_info[species_file.iloc[i, 0]] = species_file.iloc[i, 1]
        else:
            species_info[species_file.iloc[i, 0]] = ' '.join([species_file.iloc[i, 1], str(species_file.iloc[i, 2])])
    return species_info

def extract_contig_info(fafile):
    contig_info = {}
    for record in SeqIO.parse(fafile, "fasta"):
        header = record.description
        contig_id = re.search(r'lcl\|([^|]+)_cds', header).group(1)
        contig_info[record.id] = contig_id
    return contig_info

def create_output_structure(outdir, genome_id, gene_id, seq_id, seq, species_real, species_info, seqkit_path, cluster_id=None):
    out_path = f'{outdir}/{genome_id}'
    if cluster_id is not None:
        out_path = f'{out_path}/cluster{cluster_id}'
    os.makedirs(out_path, exist_ok=True)
    outfile_path = f'{out_path}/{gene_id}.fasta'
    print(f"Writing to file: {outfile_path}")
    with open(outfile_path, 'w') as outfile:
        species_info_str = species_info.get(species_real, "Unknown species")
        outfile.write(f'{seq_id} [{species_info_str}]\n{seq}\n')

    protein_outfile_path = outfile_path.replace('.fasta', '.faa')
    command = f"{seqkit_path} translate {outfile_path} -o {protein_outfile_path}"
    print(f"Translating to protein: {command}")
    subprocess.run(command, shell=True, check=True)

def process_file(file_path, refdir, outdir, gap, species_info, seqkit_path):
    gene_positions = defaultdict(list)
    with open(file_path) as f:
        for line in f:
            line_info = re.split('\s{1,}', line.strip())
            genome_id = line_info[0][:15]  # 使用第一列的前15个字符
            gene_id = line_info[1]  # 修改为第二列
            ref_file_pattern = f'{refdir}/{genome_id}*.fna'
            ref_files = glob.glob(ref_file_pattern)
            if not ref_files:
                ref_file_pattern = f'{refdir}/{genome_id}*.fa'
                ref_files = glob.glob(ref_file_pattern)
            if not ref_files:
                ref_file_pattern = f'{refdir}/{genome_id}*.fasta'
                ref_files = glob.glob(ref_file_pattern)
            if not ref_files:
                print(f"No reference files found for {genome_id}")
                continue
            ref_file = ref_files[0]
            seq_id, seq = getSeq(ref_file, gene_id)
            if seq_id == '' or seq == '':
                print(f"No sequence found for {gene_id} in {ref_file}")
                continue
            contig_id = re.search(r'lcl\|([^|]+)_cds', seq_id).group(1)
            gene_positions[genome_id].append((contig_id, seq_id, seq, gene_id))

    for genome_id, genes in gene_positions.items():
        genes.sort(key=lambda x: (x[0], x[3]))  # 按照 contig_id 和 gene_id 排序
        cluster_id = 0
        previous_gene = None
        species_real = genome_id
        for gene in genes:
            contig_id, seq_id, seq, gene_id = gene
            if previous_gene:
                prev_contig_id, prev_seq_id, prev_seq, prev_gene_id = previous_gene
                if contig_id == prev_contig_id and abs(genes.index(gene) - genes.index(previous_gene)) <= gap:
                    create_output_structure(outdir, genome_id, gene_id, seq_id, seq, species_real, species_info, seqkit_path, cluster_id)
                else:
                    cluster_id += 1
                    create_output_structure(outdir, genome_id, gene_id, seq_id, seq, species_real, species_info, seqkit_path, cluster_id)
            else:
                create_output_structure(outdir, genome_id, gene_id, seq_id, seq, species_real, species_info, seqkit_path, cluster_id)
            previous_gene = gene

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract seq')
    parser.add_argument("-i", '--refdir', required=True, help='reference dir')
    parser.add_argument("-f", '--result', required=True, help="analysis result")
    parser.add_argument("-s", '--speciesdir', required=True, help='species info dir')
    parser.add_argument("-o", '--outdir', required=True, help='output dir')
    parser.add_argument("-g", '--gap', required=False, type=int, default=5, help='gap between genes')
    parser.add_argument("--seqkit", required=True, help='path to seqkit executable')
    args = parser.parse_args()

    species_info = {}
    for file in glob.glob(f'{args.speciesdir}/*txt'):
        species_info = speciesExtract(file, species_info)

    process_file(args.result, args.refdir, args.outdir, args.gap, species_info, args.seqkit)
