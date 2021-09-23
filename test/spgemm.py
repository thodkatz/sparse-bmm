from scipy.sparse import coo_matrix, csr_matrix, csc_matrix
from scipy.io import mmread
import numpy as np
import sys

if(len(sys.argv) != 2):
    print("Usage: python spgemm.py <mtx filename>")
    exit()

coo = mmread(sys.argv[1])
coo = np.array(coo.toarray(), dtype=np.uint32)
# csr = coo.tocsr()
refMul = np.matmul(coo,coo)
boolRefMul = (refMul > 0)*1;

# print("\nTest array A:\n")
# print(coo)
# print("\nMultiplication A*A:\n")
# print(refMul)
# print("\nBoolean Multiplication A*A:\n")
# print(boolRefMul)

csrFile = open("csrMul.txt")
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
csrMul = csrMul.toarray()
#print("\nMultiplication A*A reading C++ output")
#print(csrMul)

if(np.array_equal(csrMul,boolRefMul)):
    print("Hooray! Test passed")
else:
    print("Test failed :(")
