#!/env/python3 
import json, os, re
from Bio import SeqIO
from collections import defaultdict
import sys
import glob
import argparse
import pandas as pd

def getRegion(html, product):
    region_list = []
    html = re.sub('\s', '', html)
    datas = re.findall('<divclass="regbutton(.*?)</div>', html)
    for data in datas:
        if product.lower() in data.lower():
            region = re.findall('<ahref=.*?>(.*?)<', data)
            region_id = re.findall('<ahref="#(.*?)"', data)
            region_list.append(region_id[0])
    return region_list

def getGbk(html):
    gbk_dic = {}
    html = re.sub('\s', '', html)
    html = re.sub('"', "'", html)
    ids = re.findall("<divclass='page'id='(.*?)'style", html)[1:]
    gbks = re.findall("<divclass='region-download'><ahref=(.*?)>", html)
    for id, gbk in zip(ids, gbks):
        gbk_dic[id] = gbk
    return gbk_dic


def getGene(product):
    gene_list = []
    datas = open(f'{tag_folder}/gene_lable2/{product}.txt', 'r', encoding='utf-8').readlines()
    for data in datas:
        item = data.replace('\n', '').split('/')
        gene_list.append(item) if ((not re.search('^$', data.strip())) and (item not in gene_list)) else ''
    return gene_list

def getData(path):
    with open(path, 'r', encoding='unicode_escape') as f:
        while True:
            logLine = f.readline()
            if not logLine:
                break
            yield logLine

def delSequence(gen_rs):
    location_left = []
    location_right = []
    ready_del = defaultdict(list)
    already_del = []
    del_list = []
    retain_list = []
    new_gen_rs = defaultdict(list)
    for key, value in gen_rs.items():
        for va in value:
            location_left.append(va[0][0])
            location_right.append(va[0][1])
            if len(value) == 1:
                retain_list.append(va[0])
            else:
                ready_del[key].append(va[0])
    if ready_del:
        for key in ready_del:
            for r_d in ready_del[key]:
                if r_d not in retain_list and key not in already_del:
                    if r_d[0] == min(location_left):
                        del_list.append(r_d)
                        already_del.append(key)
                    elif r_d[1] == max(location_right):
                        del_list.append(r_d)
                        already_del.append(key)
    for key, value in gen_rs.items():
        if len(value) > 1:
            for va in value:
                if va[0] not in del_list:
                    new_gen_rs[key].append(va)
        else:
            for va in value:
                new_gen_rs[key].append(va)
    return new_gen_rs

def speciesExtract(species, species_info):
    try:
        species_file = pd.read_table(species, usecols=[0, 7, 8], low_memory=False)
    except:
        species_file = pd.read_table(species, encoding='unicode_escape', usecols=[0, 7, 8], low_memory=False)
    for i in range(species_file.shape[0]):
        if str(species_file.iloc[i, 2]) == 'nan':
            species_info[species_file.iloc[i, 0]] = species_file.iloc[i, 1]
        else:
            species_info[species_file.iloc[i, 0]] = ' '.join([species_file.iloc[i, 1], str(species_file.iloc[i, 2])])
    return species_info
def getSequence(file, special_genes):
    gen_sequence = defaultdict(list)
    gbk = f'{ori_folder}/{file}'
    gbk_file = SeqIO.parse(gbk, 'genbank')
    for record in gbk_file:
        for feature in record.features:
            gene_functions = feature.qualifiers.get("gene_functions")
            location = feature.location
            is_special = False
            if gene_functions:
                if any(special_gene.lower() in gene_function.lower() for special_gene in special_genes for gene_function in gene_functions):
                    is_special = True
                    gene_cutC = feature.qualifiers.get('gene')
                    print(gene_cutC)
                    try:
                        gene_functions = gene_functions + gene_cutC
                        print(gene_functions)
                    except:
                        pass
                gene_functions = [gene.replace('\n', '').lower() for gene in gene_functions]
                gene_functions = [re.sub('biosynthetic.*?rule-based-clusters\) ', '', gene) for gene in gene_functions]
                gene_functions = sorted(gene_functions)
                if feature.type == 'CDS':
                    location = list(str(location).split(']')[0][1:].split(':'))
                    location = [int(lo.replace('>', '').replace('<', '')) for lo in location]
                    rs = feature.qualifiers.get("locus_tag")[0]
                    translation = feature.qualifiers.get("translation")[0]
                    cds = feature.extract(record.seq)
                    gen_sequence['/'.join(gene_functions).lower()].append([list(location), translation, cds, rs, is_special])
    return gen_sequence

def outputResults(new_gen_rs, result_folder, folder_name, product, region, right_name):
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    aa_txt = open(f'{result_folder}/{folder_name}_{product}_{region}_cds.fasta', 'w')
    protein_txt = open(f'{result_folder}/{folder_name}_{product}_{region}_protein.fasta', 'w')
    diff_sequence = []
    for key, value in new_gen_rs.items():
        for val in value:
            if val not in diff_sequence:
                diff_sequence.append(val)
                va = val[3]
                is_special = val[4]
                special_marker = "[SPECIAL]" if is_special else ""
                aa_txt.write(f'>{va} {key} {str([right_name])} {special_marker}\n')
                protein_txt.write(f'>{va} {key} {str([right_name])} {special_marker}\n')
                aa_txt.write(str(val[2]) + '\n')
                protein_txt.write(str(val[1]) + '\n')
    aa_txt.close()
    protein_txt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get gene info')
    parser.add_argument("-i", '--inpath', required=True, help='input path, original data path')
    parser.add_argument("-n", '--type', required=True, help="gene name")
    parser.add_argument("-o", '--out', required=True, help='output dir')
    parser.add_argument("-s", '--specialg', required=False, default='TMA,bai_operon', help='special gene type')
    args = parser.parse_args()

    tag_folder = args.inpath
    product = args.type
    res_folder = args.out
    special_genes = args.specialg.split(',')

    species_name_dic = {}
    for file in glob.glob(f'{tag_folder}/species/*txt'):
        species_name_dic = speciesExtract(file, species_name_dic)

    for folder in os.listdir(tag_folder):
        try:
            if 'results_' in folder:
                print(f'正在匹配文件夹--{folder}...')
                folder_name = folder.replace('results_', '')
                ori_folder = f'{tag_folder}/{folder}'
                ori_files = os.listdir(ori_folder)
                if 'ERR' in ori_folder:
                    species_name = folder.split('_')[-1]
                else:
                    s = folder.split('_')
                    species_name = f'{s[1]}_{s[2]}'
                try:
                    right_name = species_name_dic[species_name]
                except:
                    right_name = ''
                for file in ori_files:
                    if 'json' in file:
                        json_data = json.load(open(f"{ori_folder}/{file}", 'r'))['records'][0]['features']
                        break
                html = open(f'{ori_folder}/index.html').read()
                region_list = getRegion(html, product)
                gbk_dic = getGbk(html)
                gene_list = getGene(product)
                for region in region_list:
                    try:
                        file = gbk_dic[region]
                        gen_sequence = getSequence(file, special_genes)
                        gen_rs = defaultdict(list)
                        for gens in gene_list:
                            gens = sorted(gens)
                            j = k = 0
                            if len(gens) == 2:
                                for gen_key, gen_value in gen_sequence.items():
                                    if gens[0].lower() in gen_key.split('/') and gens[-1].lower() in gen_key.split('/'):
                                        j = 1
                                        for gen in gen_value:
                                            gen_rs['/'.join(gens)].append(gen)
                                if j == 0:
                                    for gen_key, gen_value in gen_sequence.items():
                                        if gens[0].lower() in gen_key.split('/') or gens[-1].lower() in gen_key.split('/'):
                                            k = 1
                                            for gen in gen_value:
                                                gen_rs['/'.join(gens)].append(gen)
                            else:
                                for gen_key, gen_value in gen_sequence.items():
                                    if gens[0].lower() in gen_key.split('/'):
                                        for gen in gen_value:
                                            gen_rs['/'.join(gens)].append(gen)
                        result_folder = f'{res_folder}/Core_gene_{product}/{folder_name}'
                        if len(gen_rs.keys()) == len(gene_list):
                            p = 0
                            for key, value in gen_rs.items():
                                if len(value) == 2:
                                    p += 1
                            if p != len(gen_rs.keys()):
                                new_gen_rs = delSequence(gen_rs)
                            else:
                                new_gen_rs = gen_rs
                            same_sequence = []
                            for key, value in new_gen_rs.items():
                                if value not in same_sequence:
                                    same_sequence.append(value)
                            if len(same_sequence) == len(new_gen_rs.keys()):
                                outputResults(new_gen_rs, result_folder, folder_name, product, region, right_name)
                    except Exception as e:
                        print(e)
                        pass
        except Exception as e:
            print(e)
            pass
