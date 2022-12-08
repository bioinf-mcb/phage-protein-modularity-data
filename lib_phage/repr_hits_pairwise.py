import pandas as pd
import numpy as np

def get_prob_cov(pair_hits, prob_threshold=50, ecf=False):
    """FIXME: add docstring"""
    no_hits = len(pair_hits)
    out = []
    if no_hits == 0:
        out = [float("nan"), float("nan"), float("nan"), float("nan"), no_hits]
    elif no_hits > 0:
        subject_len = pair_hits['slength'].iloc[0]
        query_len = pair_hits['qlength'].iloc[0]
        if no_hits == 1:
            hit_prob = pair_hits['prob'].iloc[0] / 100
            hit_scov = abs(pair_hits['send'] - pair_hits['sstart'] + 1) / subject_len
            hit_qcov = abs(pair_hits['qend'] - pair_hits['qstart'] + 1) / query_len
            hit_scov = hit_scov.iloc[0]
            hit_qcov = hit_qcov.iloc[0]
            hit_pident = pair_hits['pident'].iloc[0] / 100
        elif ecf:
            # take only the longest hit for ECF qcov/scov calculations
            pair_hits = pair_hits.assign(qcov = lambda x: (x.qend - x.qstart + 1) / x.qlength)
            longest_hit = pair_hits[pair_hits['qcov'] == pair_hits['qcov'].max()]

            hit_prob = longest_hit['prob'].iloc[0] / 100
            hit_scov = abs(longest_hit['send'] - longest_hit['sstart'] + 1) / subject_len
            hit_qcov = abs(longest_hit['qend'] - longest_hit['qstart'] + 1) / query_len
            hit_scov = hit_scov.iloc[0]
            hit_qcov = hit_qcov.iloc[0]
            hit_pident = longest_hit['pident'].iloc[0] / 100
        else:
            # subject
            prob_per_pos = pd.DataFrame()
            for hid, hit in pair_hits.iterrows():
                this_hit_pos = [hit['prob'] if pos in range(hit['sstart']-1, hit['send']) else 0.0 for pos in range(subject_len)]
                prob_per_pos[str(hid)] = this_hit_pos

            pident_per_pos = pd.DataFrame()
            for hid, hit in pair_hits.iterrows():
                this_hit_pos = [hit['pident'] if pos in range(hit['sstart']-1, hit['send']) else 0.0 for pos in range(subject_len)]
                pident_per_pos[str(hid)] = this_hit_pos


            prob_per_pos['max'] = prob_per_pos.max(axis=1)
            pident_per_pos['max'] = pident_per_pos.max(axis=1)
            pos_pass = (prob_per_pos['max'] >= prob_threshold)

            hit_scov = sum(pos_pass) / len(prob_per_pos)
            hit_sprob = np.mean(prob_per_pos['max'][pos_pass]) / 100
            hit_spident = np.mean(pident_per_pos['max'][pos_pass]) / 100

            # query
            prob_per_pos = pd.DataFrame()
            for hid, hit in pair_hits.iterrows():
                this_hit_pos = [hit['prob'] if pos in range(hit['qstart']-1, hit['qend']) else 0.0 for pos in range(query_len)]
                prob_per_pos[str(hid)] = this_hit_pos

            pident_per_pos = pd.DataFrame()
            for hid, hit in pair_hits.iterrows():
                this_hit_pos = [hit['pident'] if pos in range(hit['qstart']-1, hit['qend']) else 0.0 for pos in range(query_len)]
                pident_per_pos[str(hid)] = this_hit_pos


            prob_per_pos['max'] = prob_per_pos.max(axis=1)
            pident_per_pos['max'] = pident_per_pos.max(axis=1)
            pos_pass = (prob_per_pos['max'] >= prob_threshold)

            hit_qcov = sum(pos_pass) / len(prob_per_pos)
            hit_qprob = np.mean(prob_per_pos['max'][pos_pass]) / 100
            hit_qpident = np.mean(pident_per_pos['max'][pos_pass]) / 100

            hit_prob = min(hit_sprob, hit_qprob)
            hit_pident = min(hit_spident, hit_qpident)

        out = [hit_prob, hit_scov, hit_qcov, hit_pident, no_hits]

    return out
