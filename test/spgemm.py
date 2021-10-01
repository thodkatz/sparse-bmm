from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, find
from scipy.io import mmread
from pytictoc import TicToc
from termcolor import colored
import numpy as np
import sys
import time


if(len(sys.argv) != 4):
    print("Usage: python spgemm.py <mtx filename A> <mtx filename B> <mtx filename C>")
    exit()

# Python SpGEMM

#cooA = mmread(sys.argv[1])

# Read MM format
# Skip banner
def skipBanner(filename):
    readFile = open(filename, 'r')
    lines = readFile.readlines()
    numComments = 0;
    for line in lines:
        if line[0] == "%":
            numComments = numComments+1
        else:
            metadata = [int(i) for i in line.split(" ")]
            numRow, numCol, nnz = metadata[0], metadata[1], metadata[2]
            break

    return numRow, numCol, nnz, numComments

numRowA, numColA, nnzA, numCommentsA = skipBanner(sys.argv[1])
numRowB, numColB, nnzB, numCommentsB = skipBanner(sys.argv[2])
numRowF, numColF, nnzF, numCommentsF = skipBanner(sys.argv[3])

mmA = np.genfromtxt(sys.argv[1], delimiter=' ', skip_header=numCommentsA+1)
mmB = np.genfromtxt(sys.argv[2], delimiter=' ', skip_header=numCommentsB+1)
mmF = np.genfromtxt(sys.argv[3], delimiter=' ', skip_header=numCommentsF+1)

rowA,colA = mmA[:,0], mmA[:,1]
rowB,colB = mmB[:,0], mmB[:,1]
rowF,colF = mmF[:,0], mmF[:,1]

# adjust index
rowA = rowA - 1
colA = colA - 1
rowB = rowB - 1
colB = colB - 1
rowF = rowF - 1
colF = colF - 1

cooA = coo_matrix((np.ones(nnzA), (rowA,colA)), shape=(numRowA, numColA))
cooB = coo_matrix((np.ones(nnzB), (rowB,colB)), shape=(numRowB, numColB))
cooF = coo_matrix((np.ones(nnzF), (rowB,colF)), shape=(numRowF, numColF))

csrA = cooA.tocsr()
csrB = cooB.tocsr()
csrF = cooF.tocsr()

t = TicToc()
t.tic()
csrMulRef = csrF.multiply(csrA*csrB)
csrMulRef = (csrMulRef>0)*1
t.toc()

print("Reference C.nnz: " + str(csrMulRef.getnnz()))

# C++ SpGEMM

matFile="csrMul.txt"
csrFile = open(matFile)
lines = csrFile.readlines()
csrFile.close()

# check special case, zero C matrix
if(len(lines) == 1):
    pointer = np.fromstring(lines[0],dtype=int, sep=',')
    if(not np.array_equiv(pointer,csrMulRef.indptr)):
        print(colored("Test failed :(","red", attrs=["bold"]))
        print("Pointer arrays are different")
        exit()
    else:
        print(colored("Hooray! Test passed :)", "green", attrs=["bold"]))
        exit()

if(len(lines) != 2):
    print("Couldn't parse C-CSR file")
    exit()

pointer = np.fromstring(lines[0],dtype=int, sep=',')
indices = np.fromstring(lines[1],dtype=int, sep=',')
nnz     = len(indices)
n       = len(pointer)
data    = np.ones(nnz, dtype=int)
csrMul = csr_matrix((data, indices, pointer), shape=(n-1,n-1))

# Validation

if(csrMul.getnnz() != csrMulRef.getnnz()):
    print(colored("Test failed :(","red", attrs=["bold"]))
    print(colored("Incompatible size"))
    print("Reference: " + str(csrMulRef.getnnz()))
    print("Got: " + str(csrMul.getnnz()))
    exit()

if(not np.array_equiv(csrMul.indptr,csrMulRef.indptr)):
    print(colored("Test failed :(","red", attrs=["bold"]))
    print("Pointer arrays are different")
    exit()

if(not np.array_equiv(csrMul.indices,csrMulRef.indices)):
    print(colored("Test failed :(","red", attrs=["bold"]))
    print("Pointer arrays are different")
    exit()

print(colored("Hooray! Test passed :)", "green", attrs=["bold"]))
