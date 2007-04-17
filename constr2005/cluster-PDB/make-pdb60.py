from modeller import *

log.verbose()
env = environ()
sdb = sequence_db(env, seq_database_file='pdball.pir',
                  seq_database_format='PIR',
                  chains_list='ALL', minmax_db_seq_len=(30, 3000),
                  clean_sequences=True)
sdb.filter(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, -50), 
           matrix_offset = -450, max_diff_res=30, seqid_cut=60, 
           output_grp_file='pdb_60.grp', output_cod_file='pdb_60.cod')

sdb = sequence_db(env, chains_list='pdb_60.cod', seq_database_file='pdball.pir',
                  seq_database_format='PIR')
sdb.write(chains_list='ALL', seq_database_file='pdb_60.pir',
          seq_database_format='PIR')
