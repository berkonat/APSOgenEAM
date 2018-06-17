#!/bin/bash
source /RS/progs/intel/impi/4.0.3.008/bin64/mpivars.sh
source /RS/progs/intel/composerxe/bin/compilervars.sh intel64
source /RS/progs/intel/composerxe/mkl/bin/mklvars.sh intel64
make -f Makefile.parallel clean
make -f Makefile.parallel
