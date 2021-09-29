import numpy as np

# general note, everything has been converted to 0-indexed 

def readpara():

    with open('paras') as f:
        ifoffdiag, ifoffsite, numdop, ctype = f.readline().split()
        para = {
            'ifoffdiag': int(ifoffdiag),
            'ifoffsite': int(ifoffsite),
            'numdop': int(numdop),
            'ctype': ctype
        }

        dops = []
        for _ in range(para['numdop']):
            dops += [[float(val) for val in f.readline().split()]]

        para['dops'] = dops
    
    para['CONST'] = 13.605692 * 2.0 * 0.5291772
    para['latcon'] = 5.431
    return para

def readlattice(para):
    
    def testoutput():
        dop = searchdict[(0.0, 0.0, 0.0)]
        print('dopant site: {}'.format(dop))
        print('dopant NN: {}'.format(NN[dop]))

    with open('lattice.dat') as f:
    # hardcoded, the 2nd entry on the first line denotes the total # of sites
    # Then the first 13 lines are discarded

        size = int(f.readline().split()[1])

        # discard lines
        for _ in range(12):
            f.readline()

        # read site
        site = []
        for _ in range(size):
            site +=  [[ float(val) for val in f.readline().split()]]

        site = np.array(site)[:, :3]
        searchdict = { tuple(s) : i for i, s in enumerate(site)}

        
        # read NN
        NN = [[0] for _ in range(size) ]

        for s in range( size):
            # 0-indexed
            NN[s] = [int(val) - 1 for val in f.readline().split()]
            f.readline()
            f.readline()
            f.readline()

    #np.savetxt('NN', NN, fmt='%i')
    testoutput()

    para['size'] = size
    return site, NN, searchdict


def onsite_gen(sites, NN, searchdict, para):
    size = para['size']
    dops = para['dops']
    const = para['CONST']
    latcon = para['latcon']
    ctype = para['ctype']

    with open('diagcorr') as f:
        epsilon = float(f.readline())

        # 17 = 1 + 4 + 3 * 4, up to second nearest neighbors
        cccdiag = []

        for _ in range(17):

            cccdiag  += [[float(val) for val in f.readline().split()]]


    with open('offdiagcorr') as f:
        cccoff = []

        for _ in range(17):

            cccoff  += [[float(val) for val in f.readline().split()]]

    diag = np.zeros((size, 10))
    offdiag = [ 0 for _ in range(size)]
    corrchart = [ 0 for _ in range(size)]

    for dop in dops:

        # sets up the NN sequence to add correction terms
        dopid = searchdict[tuple(dop)]

        if ctype == 'dop':
            dopNN = [dopid]

        elif ctype == 'NN':
            dopNN = [dopid] + NN[dopid]

        elif ctype == '2NN':

            dopNN = [dopid] + NN[dopid]

            for NNid in NN[dopid]:
                for rNNid in NN[NNid]:
                    if rNNid != dopid:
                        dopNN += [rNNid]

        print(dopNN)

        for i, site in enumerate(sites):
            
            if i not in dopNN:
                
                # Coulomb corrections
                r = np.linalg.norm(site - dop) * latcon 
                diag[i] += [ const / (epsilon * r)] * 10

            else:
                diag[i] += cccdiag[dopNN.index(i)]

                if offdiag[i] == 0:

                    offdiag[i] = cccoff[dopNN.index(i)]

                else:

                    offdiag[i] += cccoff[dopNN.index(i)]

                corrchart[i] = 1

    
    np.savetxt('diagonal.dat', diag)

    with open('offdiag.dat', 'w') as f:

        for site in range(size):
            f.write(str(offdiag[site]) + '\n')

    np.savetxt('corrchart.dat', corrchart, fmt='%i')


# currently no offsite
def offsite_gen(para):
    size = para['size']
    offsite = np.zeros((size, 1))

    np.savetxt('offsite', offsite)

if __name__ == '__main__':
    para = readpara()
    sites, NN, searchdict = readlattice(para)
    onsite_gen(sites, NN, searchdict, para)
    offsite_gen(para)
