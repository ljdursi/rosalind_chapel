import pytest
import subprocess

modules = ["DNA", "RNA", "REVC", "FIB"]


def run_idfn(module):
    return "run_" + module


@pytest.mark.parametrize("module", modules, ids=run_idfn)
def test_python_run(module):
    lowercase = module.lower()
    script = lowercase+'.py'

    output = subprocess.check_output(['python', script, 'test.in'],
                                     cwd='./'+module)
    with open(module+'/test.out', 'r') as f:
        goodoutput = f.readlines()

    for goodline, outline in zip(goodoutput, output.splitlines()):
        assert outline.strip() == goodline.strip()


@pytest.mark.parametrize("module", modules, ids=run_idfn)
def test_chapel_run(module):
    lowercase = module.lower()
    executable = './'+lowercase

    output = subprocess.check_output([executable, '--infile', 'test.in'],
                                     cwd='./'+module)
    with open(module+'/test.out', 'r') as f:
        goodoutput = f.readlines()

    for goodline, outline in zip(goodoutput, output.splitlines()):
        assert outline.strip() == goodline.strip()
