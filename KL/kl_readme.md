----compiling and executing OpenMP codes

gcc -O3 -fopenmp filename.c

export OMP_NUM_THREADS=<number of threads>

time ./a.out dataset_filename.net

----compiling and executing CUDA codes

nvcc filename.c

time ./a.out dataset_filename.net

-----------profiling cuda code

nvprof ./a.out dataset_filename.net

*****keep the dataset file in the same directory where code files are present
