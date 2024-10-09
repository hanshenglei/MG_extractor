#!/bin/python3

import re
import os
import sys
import glob
import pandas as pd
import argparse
from collections import defaultdict

def getLable(lable):
    '''
    get gene name in lable
    '''
    lable_file = open(lable)
    lables = []
    for line in lable_file:
        lables.append(line.strip())
    lable_file.close()
    return set(lables)

def getHmmIndex(stext, target=' '):
    target_idnex = []
    i = 0
    start_r = 0
    for s in stext:
        if s == '-' and start_r == 0:
            start_r = 1
        if s == target and start_r != 0:
            target_idnex.append(i)
        i += 1
    return target_idnex

def getData(file_path, lable, score, out):
    infile = open(file_path)
    headers = {}
    score_left = 0
    score_right = 0
    for line in infile:
        header_line = line.strip()
        if line.startswith('# target name'):
            total = 0
            file_index = 0
            score_max = 0
            already = []
            sep_line = next(infile)
            sep_index = getHmmIndex(sep_line)
            if len(headers) == 0:
                start = 0
                for idx in sep_index:
                    key = line[start:idx].strip()
                    for i in range(1, 10):
                        if key in headers:
                            key += str(i)
                    headers[key] = []
                    start = idx
                headers[line[start:].strip()] = []
                headers['file'] = []
            idx_key = {}
            headers_tmp = {}
            start = 0
            for idx in sep_index:
                key = line[start:idx].strip()
                for i in range(1, 10):
                    if key in headers_tmp:
                        key += str(i)
                headers_tmp[key] = ''
                idx_key[idx] = key
                if key == 'score':
                    score_left = start
                    score_right = idx
                start = idx
            while True:
                line = next(infile).strip()
                if not line.startswith('#'):
                    score_tmp = float(line[score_left:score_right].strip())
                    if score_tmp >= score:
                        total += 1
                        start = 0
                        file_index = 1
                        already.append(line[0:sep_index[0]].strip())
                        for idx in sep_index:
                            key = idx_key[idx]
                            headers[key].append(line[start:idx].strip())
                            start = idx + 1
                        headers[header_line[start:]].append(line[start:].strip())

                elif line.startswith('# Program') and total == 0:
                    break
                elif line.startswith('# Target file') and file_index == 1:
                    file_name = os.path.basename(line.split(' ')[-1]).strip()
                    for i in range(total):
                        headers['file'].append(file_name)
                    file_index = 0
                    total = 0
                elif line.startswith('# [ok]'):
                    break
    infile.close()

    out_data = pd.DataFrame(headers)
    cols = ['GenomesID']
    cols.extend(list(out_data.columns))
    out_data['GenomesID'] = out_data['file'].apply(lambda x: '.'.join(x.split('.')[:-1]))
    out_data = out_data.loc[:, cols]

    out_data.to_csv(out, sep='\t', index=False, header=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='HMM_info_extract')
    parser.add_argument("label", help='label name')
    parser.add_argument("label_path", help="label path")
    parser.add_argument("data_path", help='data path')
    parser.add_argument("output", help='output path')
    parser.add_argument("score", type=int, help='score')
    args = parser.parse_args()

    labels = getLable(f'{args.label_path}/{args.label}.txt')

    if not os.path.exists(args.data_path):
        sys.exit(1)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    for file_path in glob.glob(f'{args.data_path}/*.tab'):
        out_name = os.path.basename(file_path)
        if re.search(args.label, file_path):
            getData(file_path, labels, args.score, f'{args.output}/{out_name}.txt') 
