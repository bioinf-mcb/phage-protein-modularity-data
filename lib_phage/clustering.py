### clustering-related functions ###

from subprocess import call

def cluster_proteins(input_fasta_filepath, output_dirpath,
                     mmseqs_tempdir,
                     cluster_params_min_seqid = 0.3,
                     cluster_params_sensitivity = 7,
                     cluster_params_coverage = 0.95,
                     verbose=False):

    clustering_filepath = output_dirpath + '/clustering.tsv'

    mmseqs_createdb_cmd = 'mmseqs createdb {} {}/SEQDB'.format(input_fasta_filepath, mmseqs_tempdir)


    mmseqs_cluster_cmd = 'mmseqs cluster {0}/SEQDB {0}/CLUDB {0} --min-seq-id {1:2f} -s {2:d} -c {3:2f}'.format(
                                  mmseqs_tempdir,
                                  cluster_params_min_seqid,
                                  cluster_params_sensitivity,
                                  cluster_params_coverage)


    mmseqs_createtsv_cmd = 'mmseqs createtsv {0}/SEQDB {0}/SEQDB {0}/CLUDB {1}'.format(
                                    mmseqs_tempdir, clustering_filepath)


    # print(mmseqs_createdb_cmd)
    # print(mmseqs_cluster_cmd)
    # print(mmseqs_createtsv_cmd)

    if verbose:
        print('Creating db...', end=' ')
    call(mmseqs_createdb_cmd, shell=True)
    if verbose:
        print('Done!')
        print('Clustering...', end=' ')
    call(mmseqs_cluster_cmd, shell=True)
    if verbose:
        print('Done!')
        print('Generating a clustering table...', end=' ')
    call(mmseqs_createtsv_cmd, shell=True)
    if verbose:
        print('Done!')

    ###!!! TODO: check file creation success (maybe outside this function)
    '''
    output.successful <- file.exists(output.filename) & file.size(output.filename)>100
    if(!output.successful) cat("Output failed to create!\n")
    '''

    return clustering_filepath
