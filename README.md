# CCCgen
Scripts for lattice generation

Default is to always generate diagonal.dat file, i.e. the diagonal corrections.
The switches in 'paras' are as follows:

if_gen_offdiag, if_gen_offsite (not supported at the moment), number of dopants, correction type (dop, 1NN, 2NN)
pos of 1st dopant
pos of 2nd dopant (if any) and so on

To use, need to generate a proper lattice.dat file (or use the provide 40x40x40 lattice)
diagcorr and offdiagcorr are user-generated correction values. Right now the optimal values are provided up to second nearest-neighbor

You would only need numpy installed on your machine. 
