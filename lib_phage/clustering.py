### clustering-related functions ###

from subprocess import call

def cluster_proteins(input_fasta_filepath, output_dirpath,
                     mmseqs_tempdir, mmseqs_binpath,
                     cluster_params_min_seqid = 0.3,
                     cluster_params_sensitivity = 7,
                     cluster_params_coverage = 0.95,
                     verbose=False):

    """Cluster proteins with mmseqs2."""

    clustering_filepath = output_dirpath + '/clustering.tsv'

    mmseqs_createdb_cmd = '{} createdb {} {}/SEQDB'.format(mmseqs_binpath, input_fasta_filepath, mmseqs_tempdir)


    mmseqs_cluster_cmd = '{0} cluster {1}/SEQDB {1}/CLUDB {1} --min-seq-id {2:2f} -s {3:d} -c {4:2f}'.format(
                                  mmseqs_binpath,
                                  mmseqs_tempdir,
                                  cluster_params_min_seqid,
                                  cluster_params_sensitivity,
                                  cluster_params_coverage)


    mmseqs_createtsv_cmd = '{0} createtsv {0}/SEQDB {0}/SEQDB {0}/CLUDB {1}'.format(
                                    mmseqs_binpath, mmseqs_tempdir, clustering_filepath)

    ### DEVEL
    mmseqs_createfasta_cmd = '{0} result2flat {0}/SEQDB {0}/SEQDB {0}/CLUDB {1}'.format(
                                    mmseqs_binpath, mmseqs_tempdir, clustering_filepath.replace('.tsv', '.fasta'))


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
