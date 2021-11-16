### functions for step validations and logs writing ###

import os
import glob
import shutil

def check_input_repr_prot_selection():

    # check input files integrity & if this step was already executed:
    # if it was warn about data overwrite

    """
    input.file.ok <- file.exists(input.fasta) & file.size(input.fasta)>100
      if(!input.file.ok){
        cat("Something wrong with the input file! Exiting...\n")
      }

      else{
        if(!dir.exists(output.dir)) dir.create(output.dir, recursive = T, showWarnings = F)
        mmseqs.tmp <- tempdir()
        if(dir.exists(mmseqs.tmp)) unlink(mmseqs.tmp, recursive = T); dir.create(mmseqs.tmp)
        # cat(mmseqs.tmp,"\n")
        output.filename <- sprintf("%s/%s", output.dir, output.filename)
        output.file.exists <- file.exists(output.filename) & file.size(output.filename)>100
        if(output.file.exists){
          if(overwrite){
            unlink(output.filename)
          } else{
            stop("Output file already exists! Exiting...\n")
          }
        }
    """

    return True

def validate_output_repr_prot_selection(work_dir, output_dirpath,
                                    cluster_params_min_seqid,
                                    cluster_params_sensitivity,
                                    cluster_params_coverage):

    """Validate if clustering was successful."""

    # validate clustering and selection process - check files integrity, write to log:
    flog = open(work_dir + 'log/repr-seq.log', 'w')
    # check files present and not empty
    for f in ['/clustering.tsv', '/name-table.txt', '/repr-seqs-lengths.txt', '/repr-seqs.fa']:
        if os.path.isfile(output_dirpath + f):
            if os.path.getsize(output_dirpath + f) != 0:
                pass
            else:
                print(f + ' file empty')
                return False
        else:
            print(f + ' does not exist')
            return False

    # write validation result
    flog.write('STATUS:OK\n')

    # write params used
    flog.write('PARAMS USED:\n')
    flog.write('cluster_params_min_seqid=' + str(cluster_params_min_seqid) + '\n')
    flog.write('cluster_params_sensitivity=' + str(cluster_params_sensitivity) + '\n')
    flog.write('cluster_params_coverage=' + str(cluster_params_coverage) + '\n')

    # paths to inputs ?

    flog.close()

    print('Validation success, log file stored.')

def check_input_all_vs_all_HMM(work_dir, force=False):

    """Check if clustering successful before starting hhblits."""

    if os.path.isfile(work_dir + 'log/repr-seq.log'):
        flog = open(work_dir + 'log/repr-seq.log', 'r')
        # read status line
        status = flog.readline().split(':')[1]
        if status.strip() == 'OK':
            pass
        else:
            print('Log file reads: previous step failed. Aborting.')
    else:
        print('Failed to fetch log file, aborting.')
        return False

    # check if data already exist
    output_hhblits_dirpath = work_dir + 'intermediate/prot-families/all-by-all'

    if force == True:
        confirm = input('This will overwrite all data from this step. Proceed? [y/n]')
        if confirm == 'y':
            #!TODO clear dir
            # clear dir
            print('Clearing ' + output_hhblits_dirpath + '...')
            return True
        else:
            return False
    else:
        if len(os.listdir(output_hhblits_dirpath) ) != 0:
            print('This step was already executed. Run validation fucntion with force=True to overwrite.')
            return False

    return True

def save_params_hhblits(work_dir, n, mact, p, qid, cov):

    # validate hhblits - save params to log:
    flog = open(work_dir + 'log/hhblits.log', 'w')

    # write validation result
    flog.write('STATUS:RUNNING\n')

    # write params used
    flog.write('PARAMS USED:\n')
    flog.write('number of iterations[n]=' + str(n) + '\n')
    flog.write('posterior prob threshold for MAC realignment[mact]=' + str(mact) + '\n')
    flog.write('minimum probability in summary and alignment list[p]=' + str(p) + '\n')
    flog.write('minimum sequence identity with master sequence[qid]=' + str(qid) + '\n')
    flog.write('minimum coverage with master sequence (%)[cov]=' + str(cov) + '\n')

    # paths to inputs ?

    flog.close()

    print('Parameters saved, log file stored.')

def validate_output_hhblits(work_dir):

    # check if hhblits step complete
    flog = open(work_dir + 'log/hhblits.log', 'r')
    status = flog.readline().split(':')[-1].strip()
    params = flog.read()
    flog.close()
    if status == 'RUNNING':
        # check number of .hhr files compared to input files
        ind_seqs_dirpath       = work_dir + 'tmp/all-by-all/individual-seqs'
        output_hhblits_dirpath = work_dir + 'intermediate/prot-families/profiles'

        input_seqs_n = len(glob.glob(ind_seqs_dirpath + '/*.fa'))
        output_a3m_n = len(glob.glob(output_hhblits_dirpath + '/*.a3m'))

        if input_seqs_n == output_a3m_n:
            flog = open(work_dir + 'log/hhblits.log', 'w')
            flog.write('STATUS:OK\n')
            flog.write(params)
            flog.close()
            print('HHblits step complete. Updated status.')
            return True
        else:
            print('HHblits step not yet finished. Please wait and try again later.')
            return False

    elif status == 'OK':
        return True

def validate_create_db(work_dir):

    flog   = open(work_dir + 'log/hh-db.log', 'r')
    status = flog.readline().split(':')[-1].strip()
    flog.close()

    # write validation result
    if status == 'OK':
        return True
    else:
        return False

def validate_search_all_vs_all(work_dir):

    # check if hhblits step complete
    flog = open(work_dir + 'log/hhblits-all-vs-all.log', 'r')
    status = flog.readline().split(':')[-1].strip()
    params = flog.read()
    flog.close()
    if status == 'RUNNING':
        # check number of .hhr files compared to input files
        ind_profile_dir        = work_dir + 'intermediate/prot-families/profiles'
        output_hhblits_dirpath = work_dir + 'intermediate/prot-families/all-by-all'

        input_prof_n = len(glob.glob(ind_profile_dir + '/*.a3m'))
        output_hhr_n = len(glob.glob(output_hhblits_dirpath + '/*.hhr'))

        # first check filesize of all hhr files
        #!TODO

        if input_prof_n == output_hhr_n:
            flog = open(work_dir + 'log/hhblits-all-vs-all.log', 'w')
            flog.write('STATUS:OK\n')
            flog.write(params)
            flog.close()
            print('HHblits all-vs-all step complete. Updated status.')
            return True
        else:
            print('HHblits all-vs-all step not yet finished. Please wait and try again later.')
            return False

    elif status == 'OK':
        return True

def validate_input_ECF(work_dir):

    """Check input files integrity & if this step was already executed:
    if it was warn about data overwrite, otherwise create dir for results
    """

    # check if hhr table and annotation table files exist & are not empty
    output_dirpath  = work_dir + 'output/prot-families/all-by-all'
    ecf_out_dirpath = work_dir + 'output/ecf-search/'

    for f in ['/table-hhr.txt', '/repr-annot.txt']:
        if os.path.isfile(output_dirpath + f):
            if os.path.getsize(output_dirpath + f) != 0:
                pass
            else:
                print(f + ' file empty')
                return False
        else:
            print(f + ' does not exist')
            return False

    # check if ECF search step already executed (look for log file)
    if os.path.isfile(work_dir + '/log/ecf-search.log'):

        # if last time run in benchmark mode: clean all and run without prompt
        # if last time scan mode: ask before proceed
        flog   = open(work_dir + '/log/ecf-search.log', 'r')
        status = flog.readline().split(':')[-1].strip()

        if status != 'BENCHMARK':
            # warn about data overwrite
            print('Log file for ECF search already exists. Do you want to repeat analysis and overwrite existing data?')
            response = input('Proceed? [y/n]')
        else:
            response = 'y'
        if response == 'y':
            # perform data dir cleanup
            os.remove(work_dir + '/log/ecf-search.log')
            try:
                shutil.rmtree(ecf_out_dirpath)
            except:
                print('Error occured when performing output dir clean-up')
            # setup dirs for analysis output
            try:
                os.mkdir(ecf_out_dirpath)
                os.mkdir(ecf_out_dirpath + 'plots/')
            except FileExistsError:
                print('/output/ecf-search/ directory already set up')

    else:
        # setup dirs for analysis output
        try:
            os.mkdir(ecf_out_dirpath)
            os.mkdir(ecf_out_dirpath + 'plots/')
        except FileExistsError:
            print('/output/ecf-search/ directory already set up')

    return True

def validate_output_ECF(work_dir, mode, filters_params, filters_used, clust_eps, clust_min_sample,
                        min_coverage, domain_min_len, merge_size_cutoff,
                        merge_shared_cutoff, steps):

    ecf_out_dirpath = work_dir + 'output/ecf-search/'
    flog            = open(work_dir + '/log/ecf-search.log', 'w')

    # what exactly should be logged and saved in which mode?
    if (mode == 'scan' or mode == 'scan-min'):
        print('Validating ECF scan results...')

        if os.path.isfile(ecf_out_dirpath + 'ecfs_results'):
            if os.path.getsize(ecf_out_dirpath + 'ecfs_results') != 0:
                flog.write('STATUS:OK\n')

                # write params used
                flog.write('PARAMS USED:\n')
                flog.write('=Data filtering params=\n')
                flog.write('Probability threshold: ' + str(filters_params['prob_threshold']) + '\n')
                flog.write('E-value threshold: ' + str(filters_params['eval_threshold']) + '\n')
                flog.write('Coverage cutoff: ' + str(filters_params['coverage_cutoff']) + '\n')
                flog.write('Filters used: ' + ' '.join(filters_used) + '\n')
                flog.write('=Algorithm params=\n')
                flog.write('Clustering epsilon parameter: ' + str(clust_eps) + '\n')
                flog.write('Clustering min_sample parameter: ' + str(clust_min_sample) + '\n')
                flog.write('Minimal hit coverage at positions to consider it a domain in all steps of algorithm: ' + str(min_coverage) + '\n')
                flog.write('Minimal length [fraction of full query length] of domain that can be annotated: ' + str(domain_min_len) + '\n')
                flog.write('Ratio of domain1/domain2 and domain2/domain1 must be no smaller than this fraction: ' + str(merge_size_cutoff) + '\n')
                flog.write('Ratio of (set of common positions) to (set of union of positions) must be no smaller than this fraction: ' + str(merge_shared_cutoff) + '\n')
                flog.write('=Steps used in the algorith=\n')
                flog.write('Similiar regions will be merged according to params above: ' + str(steps['merge_similiar']) + '\n')
                flog.write('Regions with high coverage of hits and without domains will be annotated as domains at the end of algorithm: ' + str(steps['resolve_no_domain']) + '\n')

                print('Done. Validation success.')
            else:
                print('Error: ECF results file is empty.')
                flog.write('STATUS:ERROR\n')
                flog.close()
                return False
        else:
            print('Error: Results file does not exists.')
            flog.write('STATUS:ERROR\n')
            flog.close()
            return False

    elif mode == 'benchmark':

        # benchmark should be easy to quickly rerun with new params
        flog.write('STATUS:BENCHMARK\n')

    else:
        print('Unknown mode for running ECF detection. Alavilable options: benchmark, scan')
        flog.close()
        return False

    flog.close()
