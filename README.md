qsieve_omp
==========

Prima implementazione del crivello quadratico (quadratic sieve) parallela con openmp

Per compilare l' applicazione sono necessarie le librerie di omp, e gmp (GNU Multiple Precision Arithmetic Library), queste ultime dovranno essere installate con le interfacce con il linguaggio c++.

comando con cui compilo io:
g++ -03 -fopenmp -o qsieve_openmp qsieve_openmp.cpp -lgmpxx -lgmp

Prima di lanciare l applicazione è consigliato impostare i seguanti parametri:

- $ export GOMP_CPU_AFFINITY=0-"NUM TOTALI PROCESSORI", serve a bloccare un certo thread su di un processore.

- $ export OMP_NUM_THREADS="NUM PROCESSORI DA USARE", serve a decidere quanti processori usare, ad esempio impostando ad 1, l algoritmo è usato come in versione seriale.

Lancio Applicazione:
./qsieve_openmp "NUM RSA"

è consigliato reindirizzare lo standard output su di un file per memorizzare i tempi di risoluzione.

