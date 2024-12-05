configfile: 'config.yaml'

import re
import math

import numpy as np
import pandas as pd

run_list = pd.read_csv(config['run-list'], header = 0)


group_info = pd.read_csv(config['group.file'],
                        sep = '\t', header = None,
                        names = ['group', 'chrom', 'pos', 'ref', 'alt', 'weight'])


rule create_chunk_groups:
    output:
        'output/group-info/chr{chrom}/chunk_{num}'
    run:
        chunk_size = config['groups.chunk.size']
        index = int(wildcards.num)
        chrom_groups = group_info.loc[group_info['chrom'].astype(str) == str(wildcards.chrom), 'group'].drop_duplicates()
        split_groups = np.array_split(chrom_groups, range(chunk_size, len(chrom_groups), chunk_size))
        split_groups[index].to_csv(output[0], header=False, index=False)


def create_all_chunk_groups_input(w):
    chrom_groups = group_info.loc[group_info['chrom'].astype(str) == str(w.chrom), 'group'].drop_duplicates()
    chunk_size = config['groups.chunk.size']
    chunk_count = math.ceil(len(chrom_groups) / chunk_size)
    return expand(rules.create_chunk_groups.output[0], chrom=w.chrom, num=range(0, chunk_count))

checkpoint create_all_chunk_groups:
    input:
        create_all_chunk_groups_input
    output:
        'output/group-info/chr{chrom}/all-chunk-count'
    run:
        with open(output[0], 'w') as f:
            f.write(str(len(input)))
