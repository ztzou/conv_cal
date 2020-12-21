#!/usr/bin/env python

import os
import re
import sys
import time

import ete3
import numpy
import scipy.linalg

aas = 'ARNDCQEGHILKMFPSTWYV'

in_dir = sys.argv[1]
out_dir = sys.argv[2]
group_file = sys.argv[3]
mat_file = sys.argv[4]
if len(sys.argv) > 5:
    freq_mode = sys.argv[5]
else:
    freq_mode = 'gene'


def main():
    print(f'[{time.asctime(time.localtime())}] Start ...')
    dir_list = sorted([x for x in os.listdir(in_dir) 
                       if os.path.isdir(f'{in_dir}/{x}')])
    print(f'{len(dir_list)} genes to process')
    with open(group_file) as f:
        group_list = [x.split() for x in f.read().strip().split('\n')]
    print(f'{len(group_list)} groups to check')
    for pref in dir_list:
        gene_in_dir = f'{in_dir}/{pref}'
        tree = read_ancestral_tree(f'{gene_in_dir}/rst')
        freqs = get_freqs(tree, freq_mode)
        rates = read_rate(f'{gene_in_dir}/rates')
        smat = read_mat(mat_file)
        with open(f'{out_dir}/{pref}.obsconv.tsv', 'w') as f:
            for group in group_list:
                p_cnt, c_cnt, p_list , c_list = count_conv(group, tree)
                print((','.join(group) + '\t' + f'{p_cnt:d}\t{c_cnt:d}' 
                       + '\t(' + ','.join(p_list) + ')\t(' 
                       + ','.join(c_list) + ')'), file=f)
                print(f'[count obs.]{pref}\t{",".join(group)}\t{p_cnt}\t{c_cnt}')
        with open(f'{out_dir}/{pref}.expconv.tsv', 'w') as f:
            for group in group_list:
                p_prob, c_prob, p_list, c_list = calc_conv(group, tree, rates, freqs, smat)
                p_liststr = ','.join([f'{x:.5e}' for x in p_list])
                c_liststr = ','.join([f'{x:.5e}' for x in c_list])
                print((','.join(group) + '\t' + f'{p_prob:.8f}\t{c_prob:.8f}' 
                       + f'\t({p_liststr})\t({c_liststr})'), file=f)
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
        if i > 20:
            break
        cols = line.strip().split()
        if len(cols) != 0:
            mat[i, :len(cols)] = [float(x) for x in cols]
        i += 1
    mat += mat.T
    return mat


def get_freqs(tree, freq_mode='gene'):
    seqs = [x.sequence for x in tree.get_leaves()]
    if freq_mode == 'gene':
        freq = []
        seq = ''.join(seqs)
        for aa in aas:
            freq.append(seq.count(aa))
        freq = numpy.array(freq)
        return freq / freq.sum()
    elif freq_mode == 'site':
        d = dict(zip(aas,list(range(len(aas)))))
        freq = numpy.zeros((len(seqs[0]), len(aas)))
        for seq in seqs:
            for i in range(len(seq)):
                freq[i, d[seq[i]]] += 1.
        return freq / freq.sum(axis=1, keepdims=True)
    else:
        return numpy.loadtxt(freq_mode)


def count_conv(group, tree):
    p_cnt, c_cnt = 0, 0
    p_list , c_list = [], []
    a_seqs, b_seqs = [], []
    for branch in group:
        a_seqs.append((tree & branch).up.sequence)
        b_seqs.append((tree & branch).sequence)
    for i in range(len(a_seqs[0])):
        if (all([x[i] == a_seqs[0][i] for x in a_seqs[1:]]) and 
            all([x[i] == b_seqs[0][i] for x in b_seqs[1:]]) and 
            a_seqs[0][i] != b_seqs[0][i]):
            p_cnt += 1
            a_sites = ''.join([x[i] for x in a_seqs])
            b_sites = ''.join([x[i] for x in b_seqs])
            p_list.append(f'{i}_{a_sites}{b_sites}')
        elif (all([x[i] != y[i] for x, y in zip(a_seqs, b_seqs)]) and 
              all([x[i] == b_seqs[0][i] for x in b_seqs[1:]])):
            c_cnt += 1
            a_sites = ''.join([x[i] for x in a_seqs])
            b_sites = ''.join([x[i] for x in b_seqs])
            c_list.append(f'{i}_{a_sites}{b_sites}')
    return p_cnt, c_cnt, p_list , c_list


def calc_conv(group, tree, rates, freqs, smat):
    p_list, c_list = [], []
    seq_len = len(rates)
    d = 0.01
    mat_scale = 100
    t_scale = mat_scale / d
    if len(freqs.shape) == 1:
        freqs = numpy.repeat(freqs, seq_len, axis=0).reshape(-1, seq_len).T
    freqs = freqs / freqs.sum(axis=1, keepdims=True) 
    for idx in range(seq_len):
        qmat = smat.dot(numpy.diag(freqs[idx]))
        for i in range(len(aas)):
            qmat[i, i] = 0
            qmat[i, i] = - numpy.sum(qmat[i, :])
        scale_f = numpy.sum(freqs[idx] * numpy.diag(qmat))
        if numpy.abs(scale_f) > 0.0:
            qmat = qmat / numpy.abs(scale_f)
            pmat = scipy.linalg.expm(qmat * d) / mat_scale
        else:
            pmat = numpy.identity(smat.shape[0])
        for i in range(pmat.shape[0]):
            pmat[i, i] = 0
            pmat[i, i] = 1 - numpy.sum(pmat[i, :])
        anc = numpy.array([1.0 if x == tree.sequence[0] else 0. for x in aas])
        totalmat_list = []
        for b_name in group:
            a_root_dist = (tree & b_name).up.get_distance(tree)
            ab_dist = (tree & b_name).dist
            a_prob = evo(anc, a_root_dist * t_scale, pmat)
            imat = numpy.identity(smat.shape[0])
            b_condmat = evo(imat, ab_dist * t_scale, pmat)
            ab_totalmat = numpy.multiply(a_prob.reshape(-1, 1),b_condmat)
            totalmat_list.append(ab_totalmat)
        pp, pc = sum_prob(totalmat_list)
        p_list.append(pp)
        c_list.append(pc)
    return sum(p_list), sum(c_list), p_list, c_list


def evo(n0,t,mat):
    t = int(round(t))
    n = numpy.dot(n0, numpy.linalg.matrix_power(mat,t))
    return n

def sum_prob(tmat_list):
    ppmat = numpy.ones_like(tmat_list[0])
    pcvec = numpy.ones(ppmat.shape[0])
    for tmat in tmat_list:
        tmat1 = tmat.copy()
        for i in range(ppmat.shape[0]):
            tmat1[i, i] = 0.0
        ppmat *= tmat1
        tvec1 = tmat1.sum(axis=0)
        pcvec *= tvec1
    pp = ppmat.sum()
    pc = pcvec.sum()
    return pp, pc


if __name__ == "__main__":
    main()