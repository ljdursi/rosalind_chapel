import pytest
import subprocess
import os
import os.path
import time


def time_chapel(subdir):
    lowercase = subdir.lower()
    testin = 'timing.in'
    executable = './' + lowercase

    # if executable exists, run the chapel
    assert os.path.isfile(subdir+'/'+executable)
    start_time = time.time()
    chploutput = subprocess.check_call([executable, '--infile', testin],
                                       cwd='./'+subdir)
    assert chploutput == 0
    end_time = time.time()
    return end_time - start_time


def time_python(subdir):
    lowercase = subdir.lower()
    script = lowercase+'.py'
    testin = 'timing.in'

    # if executable exists, run the chapel
    assert os.path.isfile(subdir+'/'+script)
    start_time = time.time()
    pyoutput = subprocess.check_call(['python', script, testin],
                                     cwd='./'+subdir)
    assert pyoutput == 0
    end_time = time.time()
    return end_time - start_time


def test_timings():
    subdirs = [d for d in os.listdir('.')
               if os.path.isdir(os.path.join('.', d))]

    with open('timings.txt', 'w') as f:
        for subdir in subdirs:
            print(subdir)
            testdir = os.path.join('.', subdir)
            testinput = os.path.join(testdir, 'timing.in')
            print(testinput)
            if os.path.isfile(testinput):
                chpltime = time_chapel(subdir)
                pytime = time_python(subdir)
                f.write("\t".join([subdir, str(chpltime), str(pytime)])+"\n")
