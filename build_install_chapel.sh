#!/bin/bash
cd || exit
VERSION="1.17.1"
if [[ ! -f ~/chapel-${VERSION}/util/setchplenv.bash ]]
then
    wget https://github.com/chapel-lang/chapel/releases/download/${VERSION}/chapel-${VERSION}.tar.gz
    tar xzf chapel-${VERSION}.tar.gz
    rm chapel-${VERSION}.tar.gz
    cd chapel-${VERSION} || exit
    source util/setchplenv.bash
    make
    export CHPL_AUX_FILESYS=curl
    make
fi
