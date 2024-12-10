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
focal_sp_str = sys.argv[4]
sister_sp_str = sys.argv[5]
ref_sp_str = sys.argv[6]

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

            focal_sp_seq, sister_sp_seq, ref_sp_seq = [], [], []
            focal_sp_list = focal_sp_str.split(',')
            sister_sp_list = sister_sp_str.split(',')
            ref_sp_list = ref_sp_str.split(',')
            for m in focal_sp_list:
                focal_sp_seq.append(sp_seq_dict[m])
            for n in sister_sp_list:
                sister_sp_seq.append(sp_seq_dict[n])
            for r in ref_sp_list:
                ref_sp_seq.append(sp_seq_dict[r])
            m_cnt, n_cnt, m_list , n_list = count_conv(focal_sp_seq, sister_sp_seq, ref_sp_seq)
            print((f'{pref}' + '\t' + f'{m_cnt:d}\t{n_cnt:d}' 
                    + '\t(' + ','.join(m_list) + ')\t(' 
                    + ','.join(n_list) + ')'), file=f)
            print(f'[count obs.]{pref}\t{m_cnt}\t{n_cnt}')
    print(f'[{time.asctime(time.localtime())}] Finish.')

def count_conv(focal_sp, sister_sp, ref_sp):
    f_cnt, s_cnt = 0, 0
    f_list, s_list = [], []
    f_seqs, s_seqs, ref_seqs = focal_sp, sister_sp, ref_sp
    f_seqs_group, s_seqs_group = aas_classify(focal_sp), aas_classify(sister_sp)
    for i in range(len(f_seqs[0])):
        f_single_site = [x[i] for x in f_seqs_group]
        f_most_occur = max(f_single_site, key=f_single_site.count)
        f_unique = list(set(f_single_site))
        s_single_site = [x[i] for x in s_seqs_group]
        s_most_occur = max(s_single_site, key=s_single_site.count)
        s_unique = list(set(s_single_site))

        if (len(f_unique) <= 2 and
                all([x[i] == s_seqs[0][i] for x in s_seqs[1:]]) and
                f_most_occur != s_seqs_group[0][i] and
                ref_seqs[0][i] == s_seqs[0][i]):
            f_cnt += 1
            f_sites = ''.join([x[i] for x in f_seqs])
            s_sites = ''.join([x[i] for x in s_seqs])
            r_sites = ''.join([x[i] for x in ref_seqs])
            f_list.append(f'{i}_{f_sites}{s_sites}{r_sites}')

        if (all([x[i] == f_seqs[0][i] for x in f_seqs[1:]]) and
                len(s_unique) <= 2 and
                f_seqs_group[0][i] != s_most_occur and
                ref_seqs[0][i] == f_seqs[0][i]):
            s_cnt += 1
            f_sites = ''.join([x[i] for x in f_seqs])
            s_sites = ''.join([x[i] for x in s_seqs])
            r_sites = ''.join([x[i] for x in ref_seqs])
            s_list.append(f'{i}_{f_sites}{s_sites}{r_sites}')

    return f_cnt, s_cnt, f_list, s_list

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
