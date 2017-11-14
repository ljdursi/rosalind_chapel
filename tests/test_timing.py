import pytest
import subprocess
import os
import os.path


def time_chapel(subdir):
    lowercase = subdir.lower()
    testin = 'timing.in'
    executable = './' + lowercase

    # if executable exists, run the chapel
    assert os.path.isfile(subdir+'/'+executable)
    chploutput = subprocess.Popen(['/usr/bin/time', '-f', '%x %U %M',
                                   executable, '--infile', testin],
                                  cwd='./'+subdir,
                                  stderr=subprocess.PIPE)

    result = chploutput.communicate()[1]
    items = result.split()
    assert int(items[0]) == 0
    time = float(items[1])
    maxmem = int(items[2])
    return time, maxmem


def time_python(subdir):
    lowercase = subdir.lower()
    script = lowercase+'.py'
    testin = 'timing.in'

    # if executable exists, run the chapel
    assert os.path.isfile(subdir+'/'+script)
    pyoutput = subprocess.Popen(['/usr/bin/time', '-f', '%x %U %M',
                                 'python', script, testin],
                                cwd='./'+subdir,
                                stderr=subprocess.PIPE)
    result = pyoutput.communicate()[1]
    items = result.split()
    assert int(items[0]) == 0
    time = float(items[1])
    maxmem = int(items[2])
    return time, maxmem


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
                chpltime, chplmem = time_chapel(subdir)
                pytime, pymem = time_python(subdir)
                f.write("\t".join([subdir, str(chpltime), str(chplmem),
                                   str(pytime), str(pymem)])+"\n")
