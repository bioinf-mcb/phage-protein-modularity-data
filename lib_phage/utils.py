### various pipeline functions ###

import os
import glob
import pandas as pd
from subprocess import call
from Bio import SeqIO
from csb.bio.io import HHOutputParser

def setup_dir_tree(work_dir):

    """Setup all directories used in the pipeline."""

    ### function should be commented with dir tree description ###

    ### tmp dirs
    for d in ['tmp', 'tmp/repr-proteins', 'tmp/mmseqs', 'tmp/all-by-all',
              'tmp/all-by-all/individual-seqs']:
        try:
            os.mkdir(work_dir + d)
        except FileExistsError:
            print('/' + d + '/ directory already set up')

    ### outputs
    for d in ['output', 'output/prot-families', 'output/prot-families/representative',
              'output/prot-families/all-by-all']:
        try:
            os.mkdir(work_dir + d)
        except FileExistsError:
            print('/' + d + '/ directory already set up')

    ### intermediates
    for d in ['intermediate', 'intermediate/prot-families', 'intermediate/prot-families/profiles',
              'intermediate/prot-families/all-by-all', 'intermediate/prot-families/db']:
        try:
            os.mkdir(work_dir + d)
        except FileExistsError:
            print('/' + d + '/ directory already set up')

    ### logs
    try:
        os.mkdir(work_dir + 'log')
    except FileExistsError:
        print('/log/ directory already set up')

def fetch_and_rename_protein_ids(work_dir, clustering_filepath, cds_all_filepath):

    """Rename fasta files from NCBI ids to project internal ids."""

    clustering_results = pd.read_csv(clustering_filepath, sep='\t', header=None)

    # get unique protein ids
    repr_prot_names = clustering_results[0].unique()

    # !!! minimize number of proteins for quick code check
    repr_prot_names = repr_prot_names[:16]
    print('DEVEL: Restricting input to 16 proteins for fast calculation')

    repr_prot_names_filepath = '{}/repr-names'.format(work_dir + 'tmp/repr-proteins')

    f = open(repr_prot_names_filepath, 'w')
    for prot in repr_prot_names:
        f.write(prot + '\n')
    f.close()

    no_repr_prot = len(repr_prot_names)


    repr_prot_seqs_filepath = '{}/repr-seqs.fa'.format(work_dir + 'tmp/repr-proteins')
    extract_repr_cmd = 'seqtk subseq {} {} > {}'.format(cds_all_filepath,
                                                        repr_prot_names_filepath,
                                                        repr_prot_seqs_filepath)
    call(extract_repr_cmd, shell=True)

    # create table mapping reprseq-like id to ncbi id
    name_table_filepath = '{}/name-table.txt'.format(work_dir + 'output/prot-families/representative')

    f = open(name_table_filepath, 'w')
    f.write('repr.name,prot.name\n')
    protid2reprid = {}
    for repr, prot in enumerate(repr_prot_names):
        reprid = 'reprseq{:0>{id_width}}'.format(repr+1, id_width=len(str(no_repr_prot)))
        f.write(reprid + ',' + prot + '\n')
        protid2reprid[prot] = reprid
    f.close()

    # rename proteins in fasta file and store file with seqs lengths
    fasta                         = SeqIO.parse(repr_prot_seqs_filepath, "fasta")
    repr_fasta                    = []
    repr_prot_seqs_filepath_final = '{}/repr-seqs.fa'.format(work_dir + 'output/prot-families/representative')
    repr_prot_lens_filepath_final = '{}/repr-seqs-lengths.txt'.format(work_dir + 'output/prot-families/representative')
    repr_prot_lens = pd.DataFrame({'name':[], 'length':[]})

    for record in fasta:
        record.id          = protid2reprid[record.id]
        record.name        = ''
        record.description = ''
        repr_prot_lens     = repr_prot_lens.append({'name':record.id, 'length':len(record)}, ignore_index=True)
        repr_fasta.append(record)
    SeqIO.write(repr_fasta, repr_prot_seqs_filepath_final, "fasta")

    repr_prot_lens['length'] = repr_prot_lens['length'].astype('int')
    repr_prot_lens           = repr_prot_lens.sort_values('name')
    repr_prot_lens.to_csv(repr_prot_lens_filepath_final, index=False)

    return no_repr_prot, name_table_filepath

def build_hhr_table(work_dir):

    """Build a table of results from hhr files."""

    output_hhblits_dirpath = work_dir + 'intermediate/prot-families/all-by-all'

    name_table_filepath = '{}/name-table.txt'.format(work_dir + 'output/prot-families/representative')
    name_table          = pd.read_csv(name_table_filepath)
    prot_n              = len(name_table)

    hhr_table_filpath   =  '{}/table-hhr.txt'.format(work_dir + 'output/prot-families/all-by-all')
    ftable              = open(hhr_table_filpath, 'w')
    ftable.write('qname,qstart,qend,qlength,sname,sstart,send,slength,pident,bitscore,eval,prob,pval\n') # write header

    for fhhr in sorted(glob.glob(output_hhblits_dirpath + '/*.hhr')):
        qname    = fhhr.split('/')[-1].split('.')[0]
        parser   = HHOutputParser()
        hit_list = parser.parse_file(fhhr)
        for hit in hit_list:
            record = ','.join([ str(i) for i in [qname, hit.qstart, hit.qend,
                               hit.qlength, hit.id, hit.start, hit.end, hit.length,
                               int(hit.identity), hit.score, hit.evalue, (hit.probability * 100),
                               hit.pvalue]])
            ftable.write(record + '\n')
    ftable.close()


def setup_paths():
    '''load existing configuration from file or create new if no file'''

    data_dir  = ''
    bin_paths = {}
    temp_dir  = ''

    # check if file config.py exists
    # if True: load paths from there to variables
    # if False: create file and ask user to state paths
    # return loaded/created data

    return data_dir, bin_paths, temp_dir
