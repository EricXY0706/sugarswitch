from utils import MsaFileGenerator

msa = MsaFileGenerator()
msa.run_mmseqs2(x="ACGTACGTACGTACGTAGCTGATCGTAGCTGATCGTAGTCGTAGTGTGATGCTAGTCGAGCTGACT", prefix="/sddn/yyf_work/sugarswitch/test")