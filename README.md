# Parallel Sparse Boolean Matrix Multiplication

Boolean Matrix Multiplication, using a custom blocked data structure similar to the CSR/CSC, is implemented for three different parallel configurations: 1) OpenMP, 2) OpenMPI and 3) OpenMP/OpenMPI. More about the project can be found at this [report](https://github.com/thodkatz/sparse-bmm/blob/master/report/report.pdf).

## Serial

```shell
make serial
./bin/serial matrices/A.mtx matrices/B.mtx matrices/F.mtx
```

## OpenMP

```shell
make openmp
export OMP_NUM_THREADS=<threads>
./bin/openmp matrices/A.mtx matrices/B.mtx matrices/F.mtx
```

## OpenMPI

```shell
make openmpi
mpirun -n <processes>./bin/openmpi matrices/A.mtx matrices/B.mtx matrices/F.mtx
```

## Hybrid OpenMP/OpenMPI

```shell
make hybrid
export OMP_NUM_THREADS=<threads>
mpirun -n <processes>./bin/openmpi matrices/A.mtx matrices/B.mtx matrices/F.mtx
```

## Input
To create random matrices:
```shell
mkdir matrices
cd test/

# edit the python file to choose the dimensions of the arrays
python mtxCreate.mtx 
```

## Validation

After successfully executing one of the available versions, for validation run:

```shell
cd test/
python spgemm.py
```