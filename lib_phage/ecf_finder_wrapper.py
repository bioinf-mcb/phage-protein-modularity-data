### API for domain/ECF finder algorithm ###

import pandas as pd
import pickle

def load_and_filter_data(work_dir, filters_used, filters_params):

    ### load data
    output_dirpath   = work_dir + 'output/prot-families/all-by-all/'
    hhr_table        = pd.read_csv(output_dirpath + 'table-hhr.txt', sep=',')
    annotation_table = pd.read_csv(output_dirpath + 'repr-annot.txt', sep=',')

    ### get dataset size
    dataset_size = len(list(hhr_table['qname'].unique())) # determine number of proteins in dataset
    id_width     = len(str(dataset_size)) # determine lenght of protein identifiers

    ### filter data according to previously set parameters
    # by probability
    if 'prob' in filters_used:
        hhr_table = hhr_table[hhr_table['prob'] >= filters_params['prob_threshold']]

    # by e-val
    if 'eval' in filters_used:
        hhr_table = hhr_table[hhr_table['eval'] <= filters_params['eval_threshold']]

    # by coverage
    if 'cov' in filters_used:
        hhr_table = hhr_table.assign(cov = lambda x: (x.qend - x.qstart + 1) / x.qlength)
        hhr_table = hhr_table[hhr_table['cov'] >= filters_params['coverage_cutoff']]

    # by-directional coverage, qcov and scov
    if 'cov-both' in filters_used:
        hhr_table = hhr_table.assign(qcov = lambda x: (x.qend - x.qstart + 1) / x.qlength)
        hhr_table = hhr_table.assign(scov = lambda x: (x.send - x.sstart + 1) / x.slength)
        hhr_table = hhr_table[hhr_table['qcov'] <= filters_params['coverage_cutoff']]
        #hhr_table = hhr_table[hhr_table['scov'] <= filters_params['coverage_cutoff']]

    if 'hit-len' in filters_used:
        hhr_table = hhr_table.assign(hlen = lambda x: (x.qend - x.qstart + 1))
        hhr_table = hhr_table[hhr_table['hlen'] >= filters_params['min_hit_len']]

    # backup self-hits (for singleton nodes creation when needed)
    hhr_self = hhr_table[hhr_table['sname'] == hhr_table['qname']]

    # exclude self-hits
    if 'self' in filters_used:
        hhr_table = hhr_table[hhr_table['sname'] != hhr_table['qname']]

    # exclude self-hits but only if they encompass full protein
    if 'self-full' in filters_used:
        hhr_self_not_full = hhr_self[(hhr_self['qend'] - hhr_self['qstart'] + 1) != hhr_self['qlength']]
        hhr_table         = hhr_table[hhr_table['sname'] != hhr_table['qname']]
        hhr_table         = pd.concat([hhr_self_not_full, hhr_table])

    return hhr_table, annotation_table, dataset_size, id_width

def store_scan_results(work_dir, results, mode='txt'):

    ecf_out_dirpath = work_dir + 'output/ecf-search/'

    if mode == 'txt':
        fecf = open(ecf_out_dirpath + 'ecfs_results', 'w')
        fecf.write('qname,ecf_start,ecf_stop\n')

        for prot_id, ecfs in results.items():
            for ecf in ecfs:
                fecf.write(','.join([prot_id, str(ecf[0]), str(ecf[1])]) + '\n')
        fecf.close()

    elif mode == 'pickle':
        fecf = open(ecf_out_dirpath + 'ecfs_results', 'wb')
        pickle.dump(results, fecf)
        fecf.close()

    print('Results of ECF scan stored.')
