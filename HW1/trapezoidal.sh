gcc serial_trapezoidal_cos.c -lm
echo "------------SERIAL MODE------------"
time ./a.out
time ./a.out
time ./a.out
time ./a.out
time ./a.out
gcc -fopenmp omp_trapezoidal_cos.c -lm
export OMP_NUM_THREADS=1
echo "------------THREAD = 1------------"
time ./a.out
time ./a.out
time ./a.out
time ./a.out
time ./a.out
export OMP_NUM_THREADS=2
echo "------------THREAD = 2------------"
time ./a.out
time ./a.out
time ./a.out
time ./a.out
time ./a.out
export OMP_NUM_THREADS=4
echo "------------THREAD = 4------------"
time ./a.out
time ./a.out
time ./a.out
time ./a.out
time ./a.out
export OMP_NUM_THREADS=6
echo "------------THREAD = 6------------"
time ./a.out
time ./a.out
time ./a.out
time ./a.out
time ./a.out
export OMP_NUM_THREADS=8
echo "------------THREAD = 8------------"
time ./a.out
time ./a.out
time ./a.out
time ./a.out
time ./a.out
