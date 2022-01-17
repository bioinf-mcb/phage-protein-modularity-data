### clustering-related functions ###

from subprocess import call

def cluster_proteins(input_fasta_filepath, output_dirpath,
                     mmseqs_tempdir, mmseqs_binpath,
                     cluster_params_min_seqid = 0.3,
                     cluster_params_sensitivity = 7,
                     cluster_params_coverage = 0.95,
                     verbose=False):

    """Cluster proteins with mmseqs2."""

    clustering_filepath     = output_dirpath + '/clustering.tsv'
    clustering_msa_filepath = output_dirpath + '/clustering.msa'


    mmseqs_createdb_cmd = '{} createdb {} {}/SEQDB'.format(mmseqs_binpath, input_fasta_filepath, mmseqs_tempdir)

    mmseqs_cluster_cmd = '{0} cluster {1}/SEQDB {1}/CLUDB {1} --min-seq-id {2:2f} -s {3:d} -c {4:2f}'.format(
                                  mmseqs_binpath,
                                  mmseqs_tempdir,
                                  cluster_params_min_seqid,
                                  cluster_params_sensitivity,
                                  cluster_params_coverage)

    mmseqs_createtsv_cmd = '{0} createtsv {1}/SEQDB {1}/SEQDB {1}/CLUDB {2}'.format(
                                    mmseqs_binpath, mmseqs_tempdir, clustering_filepath)

    mmseqs_result2msa = '{0} result2msa {1}/SEQDB {1}/SEQDB {1}/CLUDB {1}/DB_clu_msa --msa-format-mode 3'.format(
                                  mmseqs_binpath,
                                  mmseqs_tempdir)

    mmseqs_createseqfiledb_cmd = '{0} createseqfiledb {1}/SEQDB {1}/CLUDB {1}/CLUSEQDB '.format(
                                  mmseqs_binpath,
                                  mmseqs_tempdir)

    mmseqs_apply_cmd = '{0} apply {1}/CLUSEQDB {1}/CLUSEQDB_msa -- clustalo -i - --threads=1 '.format(
                                  mmseqs_binpath,
                                  mmseqs_tempdir)

    mmseqs_copy_msa = 'cp {0}/CLUSEQDB_msa {1}'.format(mmseqs_tempdir, clustering_msa_filepath)


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
    if verbose:
        print('Creating MSA...', end=' ')
    call(mmseqs_createseqfiledb_cmd, shell=True)
    call(mmseqs_apply_cmd, shell=True)
    call(mmseqs_copy_msa, shell=True)
    if verbose:
        print('Done!')

    return clustering_filepath, clustering_msa_filepath
