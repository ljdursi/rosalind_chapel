import subprocess
import os.path
import pytest


modules = ["CONS", "FIBD", "HAMM", "LIA", "PROT",
           "DNA", "GC", "IEV", "LEXF", "MPRT", "PRTM", "RNA", "SUBS",
           "FIB", "IPRB", "LGIS", "MRNA", "REVC", "SPLC", "PMCH", "TREE",
           "LONG", "SSET", "LEXV", "INOD", "PDST", "KMER", "KMP",
           "GRPH", "LCSM", "ORF", "PERM", 'PPER', 'PROB', 'REVP',
           'SIGN', 'SSEQ', 'TRAN', "REAR", "SORT", "RSTR", "LCSQ",
           "CAT", "MMCH", "CORR", "EDIT", "EVAL", "SPEC", "SCSP",
           "TRIE", "MOTZ", "ASPC", "NWCK", "DBRU"]


def run_idfn(module):
    return "run_" + module


def is_float(s):
    try:
        float(s)
        if s == "NAN":    # there has to be a better way
            return False
        return True
    except ValueError:
        return False


def is_equal_abs(x, y, abs_error=0.001):
    return abs(x-y) <= abs_error


def compare_lines(runoutput, goodoutput, anyorder=False):
    if anyorder:
        runoutput, goodoutput = sorted(runoutput), sorted(goodoutput)
    for goodline, outline in zip(goodoutput, runoutput):
        gooditems = goodline.split()
        outitems = outline.split()
        assert len(gooditems) == len(outitems)
        for gooditem, outitem in zip(gooditems, outitems):
            if is_float(gooditem):
                assert is_equal_abs(float(gooditem), float(outitem))
            else:
                assert gooditem == outitem


@pytest.mark.parametrize("module", modules, ids=run_idfn)
def test_python_chapel_run(module):
    lowercase = module.lower()
    testout = module + '/' + 'test.out'
    script = lowercase+'.py'
    executable = './' + lowercase

    # if executable exists, run the chapel
    assert os.path.isfile(module+'/'+executable)
    chploutput = subprocess.check_output([executable, '--infile', 'test.in'],
                                         cwd='./'+module)
    chploutput = chploutput.splitlines()

    # run the python
    pyoutput = subprocess.check_output(['python', script, 'test.in'],
                                       cwd='./'+module)
    pyoutput = pyoutput.splitlines()

    # if 'test.anyorder' exists, we can sort the lines
    # before comparing
    anyorder = False
    anyorderfile = module + '/' + 'test.anyorder'
    if os.path.isfile(anyorderfile):
        anyorder = True

    # if a given solution exists, compare against it
    # otherwise, compare against each other
    if os.path.isfile(testout):
        with open(testout, 'r') as outfile:
            goodoutput = outfile.readlines()
            compare_lines(pyoutput, goodoutput, anyorder)
            compare_lines(chploutput, goodoutput, anyorder)
    else:
        compare_lines(pyoutput, chploutput, anyorder)
