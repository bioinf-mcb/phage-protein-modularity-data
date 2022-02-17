### Functions for all vs all protein comparison ###

from Bio import SeqIO
from subprocess import call
import os
import glob

def save_individual_seqs(work_dir):
    '''Save each of the renamed sequences into separate file.'''

    repr_prot_seqs_filepath_final = '{}/repr-seqs.fa'.format(work_dir + 'output/prot-families/representative')
    tmp_dir                       = work_dir + 'tmp/all-by-all/individual-seqs/'

    fasta = SeqIO.parse(repr_prot_seqs_filepath_final, "fasta")
    for record in fasta:
        fasta_filename = record.id + '.fa'
        SeqIO.write(record, tmp_dir + fasta_filename, "fasta")

def run_hhblits(work_dir, hhsuite_bins, hhsuite_scripts, cpu, uniref_db_path, n, mact, p, qid, cov):

    """Run hhblits in order to build profile for each of the representative proteins."""

    # in general two options: run on queue or in background (no reason to keep notebook open)
    # write bash runfile
    bash_script_filepath   = work_dir + 'tmp/all-by-all/helper-build-profiles.sh'
    output_hhblits_dirpath = work_dir + 'intermediate/prot-families/profiles/hhblits'
    ind_seqs_dirpath       = work_dir + 'tmp/all-by-all/individual-seqs'
    out_hhr                = output_hhblits_dirpath + '/${FILE}.hhr'
    out_a3m                = output_hhblits_dirpath + '/${FILE}.a3m'
    output                 = output_hhblits_dirpath + '/${FILE}.o'

    cmd0 = '#!/bin/bash\n\n'
    cmd0 = cmd0 + 'export PATH="{}:{}:$PATH"\n\n'.format(hhsuite_bins, hhsuite_scripts)
    cmd1 = 'FILE=$(basename "${1}")\nFILE=${FILE%.*}\n'
    cmd2 = 'hhblits -cpu 1 -i $1 -d {} -o {} -oa3m {} -n {} -mact {} -p {} -z 0 -v 0 -b 0 -qid {} -cov {} &> {}\n'.format(
           uniref_db_path, out_hhr, out_a3m, n, mact, p, qid, cov, output)
    cmd3 = 'rm -rf {} {}'.format(out_hhr, out_a3m) # to activate add to list below
    fb = open(bash_script_filepath, 'w')
    for cmd in [cmd0, cmd1, cmd2]:
        fb.write(cmd)
    fb.close()
    # set execute mode
    call('chmod a+x ' + bash_script_filepath, shell=True)

    # run script
    run_cmd = 'nohup find {} -name "reprseq*fa" | xargs -P {} -n 1 {} &'.format(ind_seqs_dirpath, cpu, bash_script_filepath)
    print('DEVEL: hhblits run step omitted for quick check. Uncomment to set on.')
    # call(run_cmd, shell=True)

def run_hhblits_dbs(work_dir, hhsuite_bins, hhsuite_scripts, cpu, db_path, db_name, run_mode, n, mact, p, qid, cov):

    """Run hhblits in order to search databases with profiles constructed for reprseqs."""

    # in general two options: run on queue or in background (no reason to keep notebook open)
    # write bash runfile
    bash_script_filepath   = work_dir + 'tmp/all-by-all/helper-search-{}.sh'.format(db_name)
    output_hhblits_dirpath = work_dir + 'intermediate/prot-families/annot/{}'.format(db_name)
    try:
        os.mkdir(output_hhblits_dirpath) # create db dir
    except FileExistsError:
        print('output directory already set up')

    ind_seqs_dirpath       = work_dir + 'intermediate/prot-families/profiles/{}'.format(run_mode)
    out_hhr                = output_hhblits_dirpath + '/${FILE}.hhr'
    out_a3m                = output_hhblits_dirpath + '/${FILE}.a3m'
    output                 = output_hhblits_dirpath + '/${FILE}.o'

    cmd0 = '#!/bin/bash\n\n'
    cmd0 = cmd0 + 'export PATH="{}:{}:$PATH"\n\n'.format(hhsuite_bins, hhsuite_scripts)
    cmd1 = 'FILE=$(basename "${1}")\nFILE=${FILE%.*}\n'
    cmd2 = 'hhblits -cpu 1 -i $1 -d {} -o {} -oa3m {} -n {} -mact {} -p {} -z 0 -v 0 -b 0 -qid {} -cov {} &> {}\n'.format(
           db_path, out_hhr, out_a3m, n, mact, p, qid, cov, output)
    cmd3 = 'rm -rf {} {}'.format(out_hhr, out_a3m) # to activate add to list below
    fb = open(bash_script_filepath, 'w')
    for cmd in [cmd0, cmd1, cmd2]:
        fb.write(cmd)
    fb.close()
    # set execute mode
    call('chmod a+x ' + bash_script_filepath, shell=True)

    # run script
    run_cmd = 'nohup find {} -name "reprseq*a3m" | xargs -P {} -n 1 {} &'.format(ind_seqs_dirpath, cpu, bash_script_filepath)
    call(run_cmd, shell=True)

def build_hh_db(work_dir, hhsuite_bins, hhsuite_scripts, run_mode, verbose=False):

    """Build HH-suite database from profiles of representative sequences."""

    def validate_db_creation(work_dir, db_dirpath):

        """Validate if database was created correctly."""

        flog = open(work_dir + 'log/hh-db.log', 'w')
        # check files present and not empty
        for f in ['a3m.ffdata', 'a3m.ffindex', 'cs219.ffdata', 'cs219.ffindex', 'hmm.ffdata', 'hmm.ffindex']:
            if os.path.isfile(db_dirpath + '/all_proteins_' + f):
                if os.path.getsize(db_dirpath + '/all_proteins_' + f) != 0:
                    pass
                else:
                    print('ERROR: ' + f + ' file empty')
                    return False
            else:
                print('ERROR: ' + f + ' does not exist')
                return False

        # write validation result
        flog.write('STATUS:OK\n')
        flog.close()
        return True

    ind_profile_dir = work_dir + 'intermediate/prot-families/profiles/' + run_mode
    db_dirpath      = work_dir + 'intermediate/prot-families/db/' + run_mode

    # if db exists: ask to overwrite
    db_status = len(glob.glob(db_dirpath + '/*'))
    if db_status > 0:
        confirm = input('Database already exists. Overwrite? [y/n]')
        if confirm == 'y':
            cmd_clean_db = 'rm {}/*'.format(db_dirpath)
            call(cmd_clean_db, shell=True)
            print('Database cleaned.')

    db_status = len(glob.glob(db_dirpath + '/*'))
    if db_status == 0:
        cmd0              = 'export PATH="{}:{}:$PATH"'.format(hhsuite_bins, hhsuite_scripts)
        cmd_clear_a3m_dir = 'rm {0}/*.hhr; rm {0}/*.o'.format(ind_profile_dir)
        cmd_make_a3m      = cmd0 + ';ffindex_build -s {1}/all_proteins_a3m.ff{{data,index}} {2}'.format(hhsuite_bins, db_dirpath, ind_profile_dir)
        cmd_make_hmm      = cmd0 + ';ffindex_apply {1}/all_proteins_a3m.ff{{data,index}} -i {1}/all_proteins_hmm.ffindex -d {1}/all_proteins_hmm.ffdata -- hhmake -i stdin -o stdout -v 0'.format(hhsuite_bins, db_dirpath)
        cmd_make_cs       = cmd0 + ';cstranslate -x 0.3 -c 4 -i {1}/all_proteins_a3m -o {1}/all_proteins_cs219 -I a3m -b -f'.format(hhsuite_bins, db_dirpath)

        cmd_sort1         = 'sort -k3 -n {0}/all_proteins_cs219.ffindex | cut -f1 > {0}/sorting.dat'.format(db_dirpath)
        cmd_sort2         = cmd0 + ';ffindex_order {0}/sorting.dat {0}/all_proteins_hmm.ff{{data,index}} {0}/all_proteins_hmm_ordered.ff{{data,index}}'.format(db_dirpath)
        cmd_sort3         = 'mv {0}/all_proteins_hmm_ordered.ffindex {0}/all_proteins_hmm.ffindex'.format(db_dirpath)
        cmd_sort4         = 'mv {0}/all_proteins_hmm_ordered.ffdata {0}/all_proteins_hmm.ffdata'.format(db_dirpath)
        cmd_sort5         = cmd0 + ';ffindex_order {0}/sorting.dat {0}/all_proteins_a3m.ff{{data,index}} {0}/all_proteins_a3m_ordered.ff{{data,index}}'.format(db_dirpath)
        cmd_sort6         = 'mv {0}/all_proteins_a3m_ordered.ffindex {0}/all_proteins_a3m.ffindex'.format(db_dirpath)
        cmd_sort7         = 'mv {0}/all_proteins_a3m_ordered.ffdata {0}/all_proteins_a3m.ffdata'.format(db_dirpath)
        cmd_sort8         = 'rm {0}/sorting.dat'.format(db_dirpath)

        if len(glob.glob(ind_profile_dir + '/*.hhr')) > 0:
            call(cmd_clear_a3m_dir, shell=True)
            if verbose:
                print('Cleared .hhr and .o files from profiles dir.')
        call(cmd_make_a3m, shell=True)
        if verbose:
            print('Concatenated a3m alignments.')
        call(cmd_make_hmm, shell=True)
        if verbose:
            print('Created HMM profiles.')
        call(cmd_make_cs, shell=True)
        if verbose:
            print('Created column state (CS) sequence database.')

        for cmd in [cmd_sort1, cmd_sort2, cmd_sort3, cmd_sort4, cmd_sort5, cmd_sort6, cmd_sort7, cmd_sort8]:
            call(cmd, shell=True)
        if verbose:
            print('DB sorted.')

        if validate_db_creation(work_dir, db_dirpath):
            print('DB successfuly created.')

def run_all_vs_all(work_dir, hhsuite_bins, hhsuite_scripts, cpu, n, p, a3m_wildcard, run_mode):

    """Run HMM comparison of all representative proteins agains each other."""

    def save_log(work_dir, n, p):

        """Save log file: save parameters and set status to running (computations in background)."""

        # validate hhblits - save params to log:
        flog = open(work_dir + 'log/hhblits-all-vs-all.log', 'w')

        # set status to running
        flog.write('STATUS:RUNNING\n')

        # write params used
        flog.write('PARAMS USED:\n')
        flog.write('number of iterations[n]=' + str(n) + '\n')
        flog.write('minimum probability in summary and alignment list[p]=' + str(p) + '\n')
        flog.close()

        print('Parameters saved, log file stored.')


    bash_script_filepath = work_dir + 'tmp/all-by-all/helper-search-all.sh'
    ind_profile_dir      = work_dir + 'intermediate/prot-families/profiles/' + run_mode
    db_dirpath           = work_dir + 'intermediate/prot-families/db/' + run_mode
    output_dir           = work_dir + 'intermediate/prot-families/all-by-all/' + run_mode

    out_hhr = output_dir + '/${FILE}.hhr'
    output  = output_dir + '/${FILE}.o'

    cmd0 = '#!/bin/bash\n\n'
    cmd0 = cmd0 + 'export PATH="{}:{}:$PATH"\n\n'.format(hhsuite_bins, hhsuite_scripts)
    cmd1 = 'FILE=$(basename "${1}")\nFILE=${FILE%.*}\n'
    cmd2 = '{}/hhblits -cpu 1 -i $1 -d {}/all_proteins -o {} -n {} -p {} -z 0 -Z 32000 -v 0 -b 0 -B 32000 &> {}'.format(
    hhsuite_bins, db_dirpath, out_hhr, n, p, output)

    fb = open(bash_script_filepath, 'w')
    for cmd in [cmd0, cmd1, cmd2]:
        fb.write(cmd)
    fb.close()
    # set execute mode
    call('chmod a+x ' + bash_script_filepath, shell=True)

    run_cmd = 'nohup find {} -name "{}" | xargs -P {} -n 1 {} &'.format(ind_profile_dir, a3m_wildcard, cpu, bash_script_filepath)
    call(run_cmd, shell=True)

    save_log(work_dir, n, p)
