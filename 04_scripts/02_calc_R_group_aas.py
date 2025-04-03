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
sps1_str = sys.argv[4]
sps2_str = sys.argv[5]
mat_file = '/data1/chenshanshan/pipeline_data/00_raw_data/lg.dat'

def group_scheme(scheme):
    if scheme == 'GS1':
        group_divide = ['PAGS','CV','T','QNED','LIFWMY','HRK']
    elif scheme == 'GS2':
        group_divide = ['AGV','ILFP','YMTS','HNQW','RK','DE','C']
    elif scheme == 'GS3':
        group_divide = ['C','AGPST','NDQE','RHK','ILMV','FWY']
    elif scheme == 'GS4':
        group_divide = ['G','P','F','AILV','ST','DE','C','NQ','RH','K','WY','M']
    return group_divide

group_divide = group_scheme(input_scheme)

def main():
    print(f'[{time.asctime(time.localtime())}] Start ...')
    dir_list = sorted([x for x in os.listdir(in_dir)
                       if os.path.isdir(f'{in_dir}/{x}')])
    print(f'{len(dir_list)} genes to process')

    for pref in dir_list:
        gene_in_dir = f'{in_dir}/{pref}'
        tree = read_ancestral_tree(f'{gene_in_dir}/rst')
        freqs = get_freqs(tree)
        rates = read_rate(f'{gene_in_dir}/rates')
        smat = read_mat(mat_file)
        sps1_list = sps1_str.split(',')
        sps2_list = sps2_str.split(',')
        group_list = get_mrca_node_number(sps1_list, sps2_list, tree)
        print(group_list)
        print(f'{len(group_list)} groups to check')
        with open(f'{in_dir}/{pref}/branch_group.list', 'w') as f:
            for item in group_list:
                temp = ' '.join(item)
                f.write(temp + '\n')
        with open(f'{out_dir}/{pref}.{input_scheme}.obsconv.tsv', 'w') as f:
            for group in group_list:
                print(group)
                p_cnt, c_cnt, p_list , c_list = count_conv(group, tree)
                print((','.join(group) + '\t' + f'{p_cnt:d}\t{c_cnt:d}' 
                       + '\t(' + ','.join(p_list) + ')\t(' 
                       + ','.join(c_list) + ')'), file=f)
                print(f'[count obs.]{pref}\t{",".join(group)}\t{p_cnt}\t{c_cnt}')
        with open(f'{out_dir}/{pref}.{input_scheme}.expconv.tsv', 'w') as f:
            for group in group_list:
                p_prob, c_prob, p_list, c_list = calc_conv(group, tree, rates, freqs, smat)
                print((','.join(group) + '\t' + f'{p_prob:.8f}\t{c_prob:.8f}'), file=f)
                print(f'[calculate exp.]{pref}\t{",".join(group)}\t{p_prob:.2e}\t{c_prob:.2e}')
    print(f'[{time.asctime(time.localtime())}] Finish.')

def read_ancestral_tree(rst_file_name):
    rst_file = open(rst_file_name)
    flag0 = False
    flag1 = False
    flag2 = True
    species_list = []
    for line in rst_file:
        if (flag2 == True) and line.startswith('('):
            length_tree = ete3.Tree(line.strip())
            flag2 = False
        if flag0 == True:
            species_tree = ete3.PhyloTree(line.strip(), format=8)
            re_root = re.search(r'\)\s+([_\-\.\w]+)\s+;', line)
            if re_root:
                species_tree.name = re_root.group(1)
            for node in species_tree.traverse():
                if node.is_leaf():
                    node.name = '_'.join(node.name.split('_')[1:])
                    species_list.append(node.name)
            line_set = set(species_list + ['node',])
            flag0 = False
            flag1 = True
        if (flag1 == True) and (len(line) > 1) and (line.split()[0] in line_set):
            cols = line.strip().split()
            if cols[0] in species_list:
                (species_tree & cols[0]).sequence = ''.join(cols[1:])
            else:
                (species_tree & cols[1][1:]).sequence = ''.join(cols[2:])
        if line.startswith("tree with node labels for Rod Page's TreeView"):
            flag0 = True
    for node in species_tree.traverse('preorder'):
        leaves = set(node.get_leaf_names())
        for length_node in length_tree.traverse('preorder'):
            if set(length_node.get_leaf_names()) == leaves:
                node.dist = length_node.dist
    return species_tree

def read_rate(file_name):
    rates = []
    flag = False
    with open(file_name) as f:
        for line in f:
            cols = line.strip().split()
            if (flag == True) and (len(cols) == 5):
                rates.append(float(cols[3]))
            if 'posterior' in line:
                flag = True
    return rates

def read_mat(mat_file):
    mat = numpy.zeros((len(aas), len(aas)))
    i = 0
    for line in open(mat_file):
        if i >= 20:
            break
        cols = line.strip().split()
        if len(cols) != 0:
            mat[i+1, :len(cols)] = [float(x) for x in cols]
        i += 1
    mat += mat.T
    return mat

def get_freqs(tree):
    seqs = [x.sequence for x in tree.get_leaves()]
    # gene model
    freq = []
    seq = ''.join(seqs)
    for aa in aas:
        freq.append(seq.count(aa))
    freq = numpy.array(freq)
    return freq / freq.sum()

def count_conv(group, tree):
    p_cnt, c_cnt = 0, 0
    p_list , c_list = [], []
    a_seqs, b_seqs = [], []
    for branch in group:
        a_seqs.append((tree & branch).up.sequence)
        b_seqs.append((tree & branch).sequence)
        a_seqs_group = aas_classify(a_seqs)
        b_seqs_group = aas_classify(b_seqs)
    for i in range(len(a_seqs[0])):
        if (all([x[i] == a_seqs_group[0][i] for x in a_seqs_group[1:]]) and
            all([x[i] == b_seqs_group[0][i] for x in b_seqs_group[1:]]) and
            a_seqs_group[0][i] != b_seqs_group[0][i]):
            p_cnt += 1
            a_sites = ''.join([x[i] for x in a_seqs])
            b_sites = ''.join([x[i] for x in b_seqs])
            p_list.append(f'{i}_{a_sites}{b_sites}')
        elif (all([x[i] != y[i] for x, y in zip(a_seqs_group, b_seqs_group)]) and
              all([x[i] == b_seqs_group[0][i] for x in b_seqs_group[1:]])):
            c_cnt += 1
            a_sites = ''.join([x[i] for x in a_seqs])
            b_sites = ''.join([x[i] for x in b_seqs])
            c_list.append(f'{i}_{a_sites}{b_sites}')
    return p_cnt, c_cnt, p_list , c_list

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

def calc_conv(group, tree, rates, freqs, smat, group_divide):
    p_list, c_list = [], []
    seq_len = len(rates)
    if len(freqs.shape) == 1:
        freqs = numpy.repeat(freqs, seq_len, axis=0).reshape(-1, seq_len).T
    freqs = freqs / freqs.sum(axis=1, keepdims=True) 
    ancnode = tree.get_common_ancestor(group)

    for idx in range(seq_len):
        qmat = smat.dot(numpy.diag(freqs[idx]))
        for i in range(len(aas)):
            qmat[i, i] = 0
            qmat[i, i] = - numpy.sum(qmat[i, :])
        scale_f = numpy.sum(freqs[idx] * numpy.diag(qmat))
        flag = True
        if numpy.abs(scale_f) > 0.0:
            qmat = qmat / numpy.abs(scale_f)
        else:
            flag = False
        anc = numpy.zeros((len(group_divide),), dtype = int)

        anc = numpy.array([1.0 if x == ancnode.sequence[idx] else 0. for x in aas])   
        totalmat_list = []
        for b_name in group:
            a_root_dist = (tree & b_name).up.get_distance(ancnode) * rates[idx]
            ab_dist = (tree & b_name).dist * rates[idx]           
            if flag:  
                a_prob = evo(anc, a_root_dist, qmat)
                imat = numpy.identity(smat.shape[0])
                b_condmat = evo(imat, ab_dist, qmat)
            else:
                a_prob = anc
                b_condmat = imat
            ab_totalmat = numpy.multiply(a_prob.reshape(-1, 1), b_condmat)
            totalmat_list.append(ab_totalmat)
        pp, pc = sum_prob(totalmat_list, group_divide)
        p_list.append(pp)
        c_list.append(pc)
    return sum(p_list), sum(c_list), p_list, c_list

def evo(n0, t, mat):
    n = numpy.dot(n0, scipy.linalg.expm(mat * t))
    return n

def sum_prob(tmat_list, group_divide):
    ppmat = numpy.ones_like(tmat_list[0])
    size = ppmat.shape[0]
    pcvec = numpy.ones(size)
    group_map = {}
    for g_idx, g_str in enumerate(group_divide):
        for aa in g_str:
            group_map[aa] = g_idx
    
    mat_process = []
    for tmat in tmat_list:
        tmat1 = tmat.copy()
        for i in range(size):
            tmat1[i, i] = 0.0
        for i in range(size):
            for j in range(size):
                if tmat1[i,j] > 0:
                    if group_map[aas[i]] == group_map[aas[j]]:
                        # Substitutions within same group can not be counted 
                        tmat1[i,j] = 0.0
        mat_process.append(tmat1)
        tvec1 = tmat1.sum(axis=0)
        pcvec *= tvec1

    mat1 = mat_process[0]
    mat2 = mat_process[1]
    aa_index = {aa: i for i, aa in enumerate(aas)}
    mat_sum = 0.0
    for g in group_divide:
        g_indices = [aa_index[aa] for aa in g]
        for idx in g_indices:
            mat1_sum = mat1[:, idx].sum()
            filtered_idx = [x for x in g_indices if x != idx]
            mat2_sum = mat2[:, filtered_idx].sum()
            mat_sum += mat1_sum * mat2_sum

    pp = mat_sum
    pc = pcvec.sum()
    return pp, pc

def get_mrca_node_number(focal_sp_1, focal_sp_2, tree):
    print(focal_sp_1)
    print(focal_sp_2)
    group_list =[]
    leaf_focal_1 = []
    for leaf in tree.iter_leaves():
        if leaf.name in focal_sp_1:
            leaf_focal_1.append(leaf.name)
    focal_1_ancnode_num = tree.get_common_ancestor(leaf_focal_1).name

    leaf_focal_2 = []
    for leaf in tree.iter_leaves():
        if leaf.name in focal_sp_2:
            leaf_focal_2.append(leaf.name)
    print(leaf_focal_2)
    focal_2_ancnode_num = tree.get_common_ancestor(leaf_focal_2).name
    group_list.append([focal_1_ancnode_num, focal_2_ancnode_num])

    return group_list

if __name__ == "__main__":
    main()
