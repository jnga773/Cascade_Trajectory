# Hamiltonian in matrix form

Instead of defining the evolution part in each k-vetor loop I've defined the Hamiltonian at the very beginning of the code.

The easiest way I found was to define each operator as a separate matrix and then use MATMUL to find the final H matrix. The non-zero elements from H will be called in each k-vector loop.

I tried defining H without the other matrices but I made a mistake somewhere and couldn't figure out where each element should lie so, even though defining each matrix is slower and takes more memory, seeing as it only happens once at the beginning of the code and they aren't very big matrices I think it should be okay.
