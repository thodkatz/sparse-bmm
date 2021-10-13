from scipy.sparse import random
from numpy.random import default_rng
import numpy as np

print("Creating random sparse coo matrix A,B,F .mtx format")
print("Dir: $(PROJECT_ROOT)/matrices")

nRowA,nColA = int(4),int(4)
nRowB,nColB = nColA,nRowA
nRowF,nColF = nRowA,nColB

d = 1;
den = d/nRowA;

def printRowCol(name,nRow,nCol):
    print(name + ":" + str(nRow) + "," + str(nCol))

printRowCol('A',nRowA,nColA)
printRowCol('B',nRowB,nColB)
printRowCol('F',nRowF,nColF)

rng = default_rng()
#rng = 1
cooAsci = random(nRowA, nColA, format='coo', density=den, random_state=rng)
cooBsci = random(nRowB, nColB, format='coo', density=den, random_state=rng)
cooFsci = random(nRowF, nColF, format='coo', density=den, random_state=rng)

# print(((cooAsci.todense() > 0))*1)
# print("\n")
# print(((cooBsci.todense() > 0))*1)
# print("\n")
# print(((cooFsci.todense() > 0))*1)

nnzA, nnzB, nnzF = len(cooAsci.row), len(cooBsci.row), len(cooFsci.row)

cooA = np.empty((nnzA,2))
indA = np.lexsort((cooAsci.col,cooAsci.row))
cooA[:,0] = cooAsci.row[indA]
cooA[:,1] = cooAsci.col[indA]

cooB = np.empty((nnzB,2))
indB = np.lexsort((cooBsci.col,cooBsci.row))
cooB[:,0] = cooBsci.row[indB]
cooB[:,1] = cooBsci.col[indB]

cooF = np.empty((nnzF,2))
indF = np.lexsort((cooFsci.col,cooFsci.row))
cooF[:,0] = cooFsci.row[indF]
cooF[:,1] = cooFsci.col[indF]

cooA = cooA + 1
cooB = cooB + 1
cooF = cooF + 1

def headerInfo(nRow,nCol,nnz):
    meta = str(nRow) + " " + str(nCol) + " " + str(nnz)
    return meta

np.savetxt('../matrices/A.mtx', cooA, fmt='%u', header=headerInfo(nRowA,nColA,nnzA), comments='', delimiter=' ')
np.savetxt('../matrices/B.mtx', cooB, fmt='%u', header=headerInfo(nRowB,nColB,nnzB), comments='', delimiter=' ')
np.savetxt('../matrices/F.mtx', cooF, fmt='%u', header=headerInfo(nRowF,nColF,nnzF), comments='', delimiter=' ')

print("\nFinished")
