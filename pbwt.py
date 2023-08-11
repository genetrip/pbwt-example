import numpy as np

X = [[0, 1, 0, 1, 0, 1],
    [1, 1, 0, 0, 0, 1],
    [1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 0, 1],
    [0, 1, 0, 1, 1, 0]]

# M haplotypes
M = len(X)

# N variable sites
N = len(X[0])
ppa_matrix_np = np.zeros((M,N))
div_matrix_np = np.zeros((M,N))


# initialise positional prefix array
ppa = list(range(M))
div = [0] * M

# iterate over variants
for k in range(N):

    # setup intermediates
    a = list()
    b = list()
    d = list()
    e = list()
    p = q = k + 1

    # iterate over haplotypes in reverse prefix sorted order
    for index, match_start in zip(ppa, div):
        # current haplotype
        haplotype = X[index]
        
        # allele for current haplotype
        allele = haplotype[k]
        # update intermediates
        if match_start > p:
            p = match_start
        
        if match_start > q:
            q = match_start
        if allele == 0:
            a.append(index)
            d.append(p)
            p = 0
        else:
            b.append(index)
            e.append(q)
            q = 0

        
    ppa = a + b
    ppa_matrix_np[:,k] = np.array(ppa)
    div = d + e
    div_matrix_np[:,k] = np.array(div)
print(ppa_matrix_np)
print(div_matrix_np)
