import pytest
import subprocess

modules = ["DNA", "RNA", "REVC", "FIB", "GC", "HAMM", "IPRB", "PROT", "SUBS",
           "CONS"]


def run_idfn(module):
    return "run_" + module


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_equal_abs(x, y, abs_error=0.001):
    return abs(x-y) <= abs_error


def compare_lines(runoutput, goodoutput):
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
def test_python_run(module):
    lowercase = module.lower()
    script = lowercase+'.py'

    output = subprocess.check_output(['python', script, 'test.in'],
                                     cwd='./'+module)
    with open(module+'/test.out', 'r') as f:
        goodoutput = f.readlines()

    compare_lines(output.splitlines(), goodoutput)


@pytest.mark.parametrize("module", modules, ids=run_idfn)
def test_chapel_run(module):
    lowercase = module.lower()
    executable = './'+lowercase

    output = subprocess.check_output([executable, '--infile', 'test.in'],
                                     cwd='./'+module)
    with open(module+'/test.out', 'r') as f:
        goodoutput = f.readlines()

    compare_lines(output.splitlines(), goodoutput)
