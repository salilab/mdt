from modeller import *
import re

log.verbose()
env = environ()
sdb = sequence_db(env, seq_database_file='pdball.pir',
                  seq_database_format='PIR',
                  chains_list='ALL', minmax_db_seq_len=(30, 3000),
                  clean_sequences=True)
sdb.filter(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, -50), 
           matrix_offset = -450, max_diff_res=30, seqid_cut=60, 
           output_grp_file='pdb_60.grp', output_cod_file='pdb_60.cod')

# Make pdb_60.pir file by copying every sequence listed in pdb_60.cod
# from pdball.pir:
out = file("pdb_60.pir", "w")
codes = [line.rstrip('\r\n') for line in file("pdb_60.cod")]
codes = dict.fromkeys(codes)

pirhead = re.compile(">P1;(.*)$")
printline = False
for line in file("pdball.pir"):
    m = pirhead.match(line)
    if m:
        printline = m.group(1) in codes
    if printline:
        out.write(line)
