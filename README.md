# Lizard
### The c code of the Public Key Encryption scheme Lizard

Lizard is a public key encryption scheme based on the Learning with Errors(LWE) problem and the Learning with Rounding(LWR) problem.

The code consists of IND-CPA Lizard, IND-CCA Lizard and Ring Lizard, which are included in a single c file.

To build the code, the following commands are sufficient :

`
$ make clean 
$ make all
`

Or equivalently,

`
$ make new
`

Note that the AVX2 optimization is included in our compile options. For the implementation on the machines without AVX2 instruction, you can erase an AVX2 option -mavx2 before running the code.

