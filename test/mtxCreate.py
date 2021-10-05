from scipy.sparse import random
from numpy.random import default_rng
import numpy as np

print("Creating random sparse coo matrix A,B,F .mtx format")
print("Dir: $(PROJECT_ROOT)/matrices")

nRowA,nColA = 10,20
nRowB,nColB = nColA,nRowA
nRowF,nColF = nRowA,nColB
def printRowCol(name,nRow,nCol):
    print(name + ":" + str(nRow) + "," + str(nCol))

printRowCol('A',nRowA,nColA)
printRowCol('B',nRowB,nColB)
printRowCol('F',nRowF,nColF)

rng = default_rng()
#rng = 1
cooAsci = random(nRowA, nColA, format='coo', density=0.25, random_state=rng)
cooBsci = random(nRowB, nColB, format='coo', density=0.25, random_state=rng)
cooFsci = random(nRowF, nColF, format='coo', density=0.25, random_state=rng)

# print(((cooAsci.todense() > 0))*1)
# print("\n")
# print(((cooBsci.todense() > 0))*1)
# print("\n")
# print(((cooFsci.todense() > 0))*1)

nnzA, nnzB, nnzF = len(cooAsci.row), len(cooBsci.row), len(cooFsci.row)

# cooAtuple = list(zip(cooAsci.row,cooAsci.col))
# cooAtuple.sort(key=lambda x:x[1])
# cooBtuple = list(zip(cooBsci.row,cooBsci.col))
# cooBtuple.sort(key=lambda x:x[1])
# cooFtuple = list(zip(cooFsci.row,cooFsci.col))
# cooFtuple.sort(key=lambda x:x[1])

cooA = np.empty((nnzA,2))
indA = np.lexsort((cooAsci.row,cooAsci.col))
cooA[:,0] = cooAsci.row[indA]
cooA[:,1] = cooAsci.col[indA]

cooB = np.empty((nnzB,2))
indB = np.lexsort((cooBsci.row,cooBsci.col))
cooB[:,0] = cooBsci.row[indB]
cooB[:,1] = cooBsci.col[indB]

cooF = np.empty((nnzF,2))
indF = np.lexsort((cooFsci.row,cooFsci.col))
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
