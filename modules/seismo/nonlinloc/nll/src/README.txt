======================================================
Complete NonLinLoc distribution software package
======================================================
Unpack source files: unpack: NLL[VER]_src.tgz
To build:
Solaris:
	make all
Linux:
	make -R all
See http://alomax.net/nlloc and http://alomax.net/nlloc -> tutorials for further information
======================================================


===========================================================================
NLLoc_func_test program demonstrating running NLLoc through a function call
===========================================================================
Unpack source files: unpack: NLL[VER]_src.tgz
To build:
Solaris (not tested):
	make -f Make_NLLoc_func_test
Linux:
	gmake -Rf Make_NLLoc_func_test
Unpack demo files: unpack: NLL[VER]_func.tgz
To run:
	cd nll_func
	./run_func.sh
	# clean
	rm -rf out/*
	cd ..
See nll_func/run_func.sh for more detail.
===========================================================================
