#!/usr/bin/env python

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import re
import sys
import time
import ete3
import numpy
import scipy

aas = 'ARNDCQEGHILKMFPSTWYV'

in_dir = sys.argv[1]
out_dir = sys.argv[2]
input_scheme = sys.argv[3]

mangrove_sp = ['ama', 'rap', 'sal']
nonmangrove_sp =['sin', 'ptr', 'egr']
ref_sp = ['osa']

def group_scheme(scheme):
    if scheme == 'GS1':
        group_divide = ['PAGS','CV','T','QNED','LIFWMY','HRK','X']
    elif scheme == 'GS2':
        group_divide = ['AGV','ILFP','YMTS','HNQW','RK','DE','C','X']
    elif scheme == 'GS3':
        group_divide = ['C','AGPST','NDQE','RHK','ILMV','FWY','X']
    elif scheme == 'GS4':
        group_divide = ['G','P','F','AILV','ST','DE','C','NQ','RH','K','WY','M','X']
    return group_divide
group_divide = group_scheme(input_scheme)

def main():
    print(f'[{time.asctime(time.localtime())}] Start ...')
    dir_list = sorted([x for x in os.listdir(in_dir)
                       if x.endswith('fasta')])
    print(f'{len(dir_list)} genes to process')
    with open(f'{out_dir}/{input_scheme}.ouput.tsv', 'w') as f:
        for pref in dir_list:
            fsts = open(f'{in_dir}/{pref}').read().strip().split()
            sp_names = [x[1:].split('|')[0] for x in fsts[::2]]
            seqs = fsts[1::2]
            sp_seq_dict = dict(zip(sp_names, seqs))

            mangrove_seq, nonmangrove_seq, ref_seq = [], [], []
            for m in mangrove_sp:
                mangrove_seq.append(sp_seq_dict[m])
            for n in nonmangrove_sp:
                nonmangrove_seq.append(sp_seq_dict[n])
            for r in ref_sp:
                ref_seq.append(sp_seq_dict[r])
            m_cnt, n_cnt, m_list , n_list = count_conv(mangrove_seq, nonmangrove_seq, ref_seq)
            print((f'{pref}' + '\t' + f'{m_cnt:d}\t{n_cnt:d}' 
                    + '\t(' + ','.join(m_list) + ')\t(' 
                    + ','.join(n_list) + ')'), file=f)
            print(f'[count obs.]{pref}\t{m_cnt}\t{n_cnt}')
    print(f'[{time.asctime(time.localtime())}] Finish.')

def count_conv(mangrove, nonmangrove, ref):
    m_cnt, n_cnt = 0, 0
    m_list, n_list = [], []
    m_seqs, n_seqs, ref_seqs = mangrove, nonmangrove, ref
    m_seqs_group, n_seqs_group = aas_classify(mangrove), aas_classify(nonmangrove)
    for i in range(len(m_seqs[0])):
        m_single_site = [x[i] for x in m_seqs_group]
        m_most_occur = max(m_single_site, key=m_single_site.count)
        m_unique = list(set(m_single_site))
        n_single_site = [x[i] for x in n_seqs_group]
        n_most_occur = max(n_single_site, key=n_single_site.count)
        n_unique = list(set(n_single_site))

        if (len(m_unique) <= 2 and
                all([x[i] == n_seqs[0][i] for x in n_seqs[1:]]) and
                m_most_occur != n_seqs_group[0][i] and
                ref_seqs[0][i] == n_seqs[0][i]):
            m_cnt += 1
            m_sites = ''.join([x[i] for x in m_seqs])
            n_sites = ''.join([x[i] for x in n_seqs])
            r_sites = ''.join([x[i] for x in ref_seqs])
            m_list.append(f'{i}_{m_sites}{n_sites}{r_sites}')

        if (all([x[i] == m_seqs[0][i] for x in m_seqs[1:]]) and
                len(n_unique) <= 2 and
                m_seqs_group[0][i] != n_most_occur and
                ref_seqs[0][i] == m_seqs[0][i]):
            n_cnt += 1
            m_sites = ''.join([x[i] for x in m_seqs])
            n_sites = ''.join([x[i] for x in n_seqs])
            r_sites = ''.join([x[i] for x in ref_seqs])
            n_list.append(f'{i}_{m_sites}{n_sites}{r_sites}')

    return m_cnt, n_cnt, m_list, n_list

def aas_classify(seq):
    seq_1 = []
    for s in range(len(seq)):
        seq_among = ''
        for i in range(len(seq[s])):
            for g in range(len(group_divide)):
                if seq[s][i] in group_divide[g]:
                    seq_among = seq_among + chr(g+97)
        seq_1.append(seq_among)
    return seq_1

if __name__ == "__main__":
    main()
