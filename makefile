CC=g++

tfim.out: TFIM_proj.cpp basis.h head_proj.h matrix.h measure.h simparam.h
	$(CC) TFIM_proj.cpp -o tfim.out
