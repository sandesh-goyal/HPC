gcc -O3 -fopenmp omp_combine.c
export OMP_NUM_THREADS=1
echo "------------THREAD = 1 (SERIAL)------------"
echo "------------ITERATION 1"
time ./a.out
echo "------------ITERATION 2"
time ./a.out
echo "------------ITERATION 3"
time ./a.out
echo "------------ITERATION 4"
time ./a.out
echo "------------ITERATION 5"
time ./a.out
export OMP_NUM_THREADS=2
echo "------------THREAD = 2------------"
echo "------------ITERATION 1"
time ./a.out
echo "------------ITERATION 2"
time ./a.out
echo "------------ITERATION 3"
time ./a.out
echo "------------ITERATION 4"
time ./a.out
echo "------------ITERATION 5"
time ./a.out
export OMP_NUM_THREADS=4
echo "------------THREAD = 4------------"
echo "------------ITERATION 1"
time ./a.out
echo "------------ITERATION 2"
time ./a.out
echo "------------ITERATION 3"
time ./a.out
echo "------------ITERATION 4"
time ./a.out
echo "------------ITERATION 5"
time ./a.out
export OMP_NUM_THREADS=6
echo "------------THREAD = 6------------"
echo "------------ITERATION 1"
time ./a.out
echo "------------ITERATION 2"
time ./a.out
echo "------------ITERATION 3"
time ./a.out
echo "------------ITERATION 4"
time ./a.out
echo "------------ITERATION 5"
time ./a.out
export OMP_NUM_THREADS=8
echo "------------THREAD = 8------------"
echo "------------ITERATION 1"
time ./a.out
echo "------------ITERATION 2"
time ./a.out
echo "------------ITERATION 3"
time ./a.out
echo "------------ITERATION 4"
time ./a.out
echo "------------ITERATION 5"
time ./a.out
