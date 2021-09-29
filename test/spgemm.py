from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, find
from scipy.io import mmread
import numpy as np
import sys
from termcolor import colored
import time
from pytictoc import TicToc


if(len(sys.argv) != 2):
    print("Usage: python spgemm.py <mtx filename>")
    exit()

# Python SpGEMM

coo = mmread(sys.argv[1])
csr = coo.tocsr()

t = TicToc()
t.tic()
csrMulRef = (csr*csr>0)*1
t.toc()

# C++ SpGEMM

matFile="csrMul.txt"
csrFile = open(matFile)
lines = csrFile.readlines()
csrFile.close()
if(len(lines) != 2):
    print("Invalid file\n")
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

csrDiffs = find(csrMul-csrMulRef)
if len(csrDiffs[0])==0 and len(csrDiffs[1])==0 and len(csrDiffs[2])==0:
    print(colored("Hooray! Test passed","green", attrs=["bold"]))
else:
    print(find(csrDiff))
    print(colored("Test failed :(","red", attrs=["bold"]))
