import pytest
import subprocess

modules = ["CONS", "DNA", "FIB", "FIBD", "GC", "GRPH", "HAMM",
           "IPRB", "LCSM", "LEXF", "MPRT", "MRNA", "ORF", "PERM",
           "PPER", "PROB", "PROT", "PRTM", "REVC", "REVP", "RNA", "SIGN",
           "SPLC", "SSEQ", "SUBS", "TRAN"]

need_blas = ["IEV", "LIA"]


def compile_idfn(module):
    return "compile_" + module


@pytest.mark.parametrize("module", modules, ids=compile_idfn)
def test_basic_compilation(module):
    lowercase = module.lower()
    sourcecode = lowercase+'.chpl'
    executable = lowercase

    ret_val = subprocess.check_call(['chpl', '-o', executable, sourcecode],
                                    cwd='./'+module)
    assert ret_val == 0
