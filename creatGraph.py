import random
nodeNum = 100
edge = 0
print nodeNum
for i in range(1,nodeNum+1):
    for j in range(1,nodeNum+1):
        if i != j:
            edge = random.randint(1,20)
            print i,j,edge
