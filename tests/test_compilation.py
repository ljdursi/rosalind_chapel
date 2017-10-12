import pytest
import subprocess

modules = ["CONS", "DNA", "FIB", "FIBD", "GC", "GRPH", "HAMM",
           "IPRB", "LCSM", "LEXF", "MPRT", "MRNA", "ORF", "PERM",
           "PPER", "PROB", "PROT", "PRTM", "REVC", "REVP", "RNA", "SIGN",
           "SPLC", "SSEQ", "SUBS", "TRAN"]

modules_need_blas = ["IEV", "LIA"]


def compile_idfn(module):
    return "compile_" + module


def compile_blas_idfn(module):
    return "compile_" + module + "_blas"


@pytest.mark.parametrize("module", modules, ids=compile_idfn)
def test_basic_compilation(module):
    lowercase = module.lower()
    sourcecode = lowercase+'.chpl'
    executable = lowercase

    ret_val = subprocess.check_call(['chpl', '-o', executable, sourcecode],
                                    cwd='./'+module)
    assert ret_val == 0


@pytest.mark.parametrize("module", modules_need_blas, ids=compile_blas_idfn)
def test_basic_compilation_blas(module):
    lowercase = module.lower()
    sourcecode = lowercase+'.chpl'
    executable = lowercase

    BLASLIB = '/usr/local/opt/openblas/lib'
    BLASINC = '/usr/local/opt/openblas/include'
    ret_val = subprocess.check_call(['chpl', '-o', executable, sourcecode,
                                     '-L'+BLASLIB, '-llapack', '-lblas',
                                     '-I'+BLASINC],
                                    cwd='./'+module)
    assert ret_val == 0