import pandas as pd
import sys

if __name__ == '__main__':
    work_dir = sys.argv[1]
    run_mode = sys.argv[2]
    pair_table_id = sys.argv[3]
    prob_threshold = sys.argv[4]

    sys.path.append(sys.argv[5])
    from lib_phage.repr_hits_pairwise import get_prob_cov

    ## define paths
    output_dir = work_dir + 'output/'
    inter_dir = work_dir + 'intermediate/'
    all_by_all_output_dir = output_dir + 'prot-families/all-by-all/' + run_mode + '/'
    families_output_dir = output_dir + 'prot-families/families/'

    # load pair table from path & set output name
    parts_dir = work_dir + 'tmp/prot-families/pair_table_chunks/'
    pair_table = pd.read_csv(parts_dir + 'pair-table-' + str(pair_table_id) + '.csv', sep=',')
    output_filename = 'repr-hits-pairwise-prob' + str(prob_threshold) + '-' + str(pair_table_id) + '.txt'

    # load hhr_table
    table_hhr_filename = all_by_all_output_dir + 'table-hhr.txt'
    table_hhr = pd.read_csv(table_hhr_filename, sep=',')

    ## for each pair, calculate prob and cov
    pair_results = {}

    for i, pair in pair_table.iterrows():
        try:
            this_pair_qname = pair['qname']
            this_pair_sname = pair['sname']

            prob_ab=0
            prob_ba=0
            scov_ab=0
            scov_ba=0
            qcov_ab=0
            qcov_ba=0
            pident_ab=0
            pident_ba=0

            this_pair_hits_ab = table_hhr[(table_hhr['qname'] == this_pair_qname) & (table_hhr['sname'] == this_pair_sname)]

            this_pair_prob_cov_ab = get_prob_cov(this_pair_hits_ab, prob_threshold = prob_threshold)
            if this_pair_prob_cov_ab:
                prob_ab = this_pair_prob_cov_ab[0]
                scov_ab = this_pair_prob_cov_ab[1]
                qcov_ab = this_pair_prob_cov_ab[2]
                pident_ab = this_pair_prob_cov_ab[3]


            this_pair_hits_ba = table_hhr[(table_hhr['sname'] == this_pair_qname) & (table_hhr['qname'] == this_pair_sname)]

            this_pair_prob_cov_ba = get_prob_cov(this_pair_hits_ba, prob_threshold = prob_threshold)
            if this_pair_prob_cov_ba:
                prob_ba = this_pair_prob_cov_ba[0]
                scov_ba = this_pair_prob_cov_ba[1]
                qcov_ba = this_pair_prob_cov_ba[2]
                pident_ba = this_pair_prob_cov_ba[3]


            min_prob = 0
            if (prob_ab and prob_ba):
                min_prob = min(prob_ab, prob_ba)
            this_pair_prob = min_prob
            this_pair_scov = min(scov_ab, scov_ba)
            this_pair_qcov = min(qcov_ab, qcov_ba)
            min_pident = 0
            if (pident_ab and pident_ba):
                min_pident = min(pident_ab, pident_ba)
            this_pair_pident = min_pident

            pair_results[i] = [this_pair_qname, this_pair_sname, this_pair_qcov,
                               this_pair_scov, this_pair_prob, this_pair_pident]
        except:
            print('error in pair', pair['qname'], pair['sname'], 'omitted')

    ### now to dataframe and drop to file
    pair_results_df = pd.DataFrame.from_dict(pair_results, orient='index',
                                             columns=['qname', 'sname', 'qcov', 'scov', 'prob', 'pident'])
    pair_results_df.to_csv(families_output_dir + output_filename, index=False, float_format="%.4f")
