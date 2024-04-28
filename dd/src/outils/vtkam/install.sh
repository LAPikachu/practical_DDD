pathadd() {
    if [ -d "$1" ] && [[ ":$PYTHONPATH:" != *":$1:"* ]]; then
        PYTHONPATH="${PYTHONPATH:+"$PYTHONPATH:"}$1"
        export PYTHONPATH
    fi
}


pathadd $PWD/install/lib/python

python setup.py install --home=./install

export PYTHONPATH

