import scipy.sparse
nodeCount = 3000
print(nodeCount)
graph = scipy.sparse.rand(m=nodeCount,n=nodeCount,density=0.01).tocsr()
for i in range(nodeCount):
    for j in range(nodeCount):
        if graph[i,j] != 0:
            print i+1,j+1,int(graph[i,j]*200)

