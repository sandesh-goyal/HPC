mpicc -O3 mpi_combine.c
echo "------------PROCESS = 2 (SERIAL: ROOT+1)------------"
echo "------------ITERATION 1"
time mpirun -np 2 ./a.out
echo "------------ITERATION 2"
time mpirun -np 2 ./a.out
echo "------------ITERATION 3"
time mpirun -np 2 ./a.out
echo "------------ITERATION 4"
time mpirun -np 2 ./a.out
echo "------------ITERATION 5"
time mpirun -np 2 ./a.out
echo "------------PROCESS = 3 (ROOT+2)------------"
echo "------------ITERATION 1"
time mpirun -np 3 ./a.out
echo "------------ITERATION 2"
time mpirun -np 3 ./a.out
echo "------------ITERATION 3"
time mpirun -np 3 ./a.out
echo "------------ITERATION 4"
time mpirun -np 3 ./a.out
echo "------------ITERATION 5"
time mpirun -np 3 ./a.out
echo "------------PROCESS = 5 (ROOT+4)------------"
echo "------------ITERATION 1"
time mpirun -np 5 ./a.out
echo "------------ITERATION 2"
time mpirun -np 5 ./a.out
echo "------------ITERATION 3"
time mpirun -np 5 ./a.out
echo "------------ITERATION 4"
time mpirun -np 5 ./a.out
echo "------------ITERATION 5"
time mpirun -np 5 ./a.out
echo "------------PROCESS = 7 (ROOT+6)------------"
echo "------------ITERATION 1"
time mpirun -np 7 ./a.out
echo "------------ITERATION 2"
time mpirun -np 7 ./a.out
echo "------------ITERATION 3"
time mpirun -np 7 ./a.out
echo "------------ITERATION 4"
time mpirun -np 7 ./a.out
echo "------------ITERATION 5"
time mpirun -np 7 ./a.out
echo "------------PROCESS = 9 (ROOT+8)------------"
echo "------------ITERATION 1"
time mpirun -np 9 ./a.out
echo "------------ITERATION 2"
time mpirun -np 9 ./a.out
echo "------------ITERATION 3"
time mpirun -np 9 ./a.out
echo "------------ITERATION 4"
time mpirun -np 9 ./a.out
echo "------------ITERATION 5"
time mpirun -np 9 ./a.out
