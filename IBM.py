'''IBM.py
----------------------------------
Individual-based simulation of Lotka-Volterra dynamics with
intraspecific trait variation


Requirements:
-------------
python 2.7.6 with the packages sys, numpy, and matplotlib


Usage:
------
From the command line, invoke

python IBM.py [inputfile]

The command-line parameter [inputfile] is a text file containing the
parameters with which the simulation is run. It should always look
like this:

S: 2
w: 0.1
theta: 0.5
randomseed: 54321
max_steps: 100000
N0: 1000 1000
mu: -0.3 0.3
sigma: 0.1 0.15

Explanation:
- S: number of species (positive integer)
- w: competition width (positive float)
- theta: half-width of effective range of niche axis (positive float)
- randomseed: set random generator for full replicability (positive integer)
- max_steps: number of iterations to run simulation for (positive integer)
- N0: initial abundances (S positive integers, separated by spaces)
- mu: species' trait means (S floats, separated by spaces)
- sigma: species' trait std devs (S positive floats, separated by spaces)


Output:
-------

1) A table with individuals in the rows, and their species identity
   and trait value in the columns. One can redirect the output to a
   text file for future analysis, like this:
   
   python IBM.py [inputfile] > [outputfile]

2) A histogram plot of the trait distributions

'''


# import packages
import sys
import numpy
import matplotlib.pyplot as plt


# define ingredient functions

def alpha(trait1, trait2, w): # competition kernel (Gaussian function
    return numpy.exp(-((trait1 - trait2)**2) / (w**2)) # of trait difference)


def b(trait, theta): # intrinsic growth rate of individual with given trait,
    b_ret = 1000 # equal to 1000 if trait is in the interval [-theta, theta],
    if (trait > theta) or (trait < -theta):
        b_ret = 0 # and to 0 otherwise
    return b_ret


def main(argv):
    
    # Get input data
    inputfile = sys.argv[1] # input file
    f = open(inputfile, 'r') # open input file for reading
    inputstr = f.read().split() # read contents and split into list by words
    f.close() # close input file
    
    S = int(inputstr[1]) # number of species
    w = float(inputstr[3]) # competition width
    theta = float(inputstr[5]) # half-width of effective range of trait axis
    randomseed = int(inputstr[7]) # for setting random seed
    max_steps = int(inputstr[9]) # number of simulation steps
    N0 = [int(i) for i in inputstr[11:(11+S)]] # initial abundances
    mu = [float(i) for i in inputstr[(12+S):(12+2*S)]] # trait means / spp
    sigma = [float(i) for i in inputstr[(13+2*S):(13+3*S)]] # trait std devs / spp
    
    # Initialize community
    N = sum(N0) # total number of individuals initially
    sppid = [] # initialize vector of individuals' species identities
    trait = numpy.zeros(N) # initialize vector of individuals' trait values
    birth = numpy.zeros(N) # initialize vector of individuals' birth rates
    death = numpy.zeros(N) # initialize vector of individuals' death rates
    for sp in range(S): # create initial vector of species identities
        sppid = sppid + [sp for i in range(N0[sp])]
    sppid = numpy.array(sppid) # convert into numpy array
    for i in range(N):
        trait[i] = mu[sppid[i]] # initially, no spread in trait values
        birth[i] = b(trait[i], theta) # get intrinsic growth
    for i in range(N):
        death[i] = numpy.sum(alpha(trait[i], trait, w)) # death
    numpy.random.seed(randomseed) # seed random number generator

    # Main simulation
    for t in range(max_steps):
        L = numpy.sum(birth + death) # sum of rates of all events
        r1 = numpy.random.random() # generate random number to pick the individual
                                   # to whom the birth or death event will happen
        # first entry less than r1 in list of cumulative sums of rates?
        try: # random individual chosen, weighted by rates
            rind = numpy.where(numpy.cumsum(birth+death)/L >= r1)[0][0]
        except: # the previous line returns an error if we try to pick the very last
            rind = len(sppid) - 1 # individual; here we treat that case manually
        r2 = numpy.random.random() # does selected individual die or give birth?
        if r2 < birth[rind] / (birth[rind] + death[rind]): # birth event:
            sppid = numpy.append(sppid, sppid[rind]) # append new ind's spp identity
            trait = numpy.append(trait, numpy.random.normal(mu[sppid[rind]],
                                 sigma[sppid[rind]])) # generate its trait
            birth = numpy.append(birth, b(trait[-1], theta)) # get its birth rate
            death = numpy.append(death, 0) # and its death rate; calculated below
            death += alpha(trait, trait[-1], w) # modify death rates
            death[-1] = numpy.sum(alpha(trait, trait[-1], w)) # rate of new ind
        else: # death event:
            death -= alpha(trait, trait[rind], w) # modify death rates
            sppid = numpy.delete(sppid, rind) # erase individual from all 4 lists
            trait = numpy.delete(trait, rind) # erase individual from all 4 lists
            birth = numpy.delete(birth, rind) # erase individual from all 4 lists
            death = numpy.delete(death, rind) # erase individual from all 4 lists
    
    # Generate output
    print 'sppid' + '\t' + 'trait' # header for command-line output
    for i in range(len(sppid)): # print species identity and trait value of all
        print str(sppid[i]) + '\t' + str(trait[i]) # individuals in separate rows
    
    # Plotting; comment this whole section out if plot is not needed
    cols = ['b', 'y', 'g', 'r', 'c', 'm', 'k'] # color palette
    trait_by_spp = [[] for sp in range(S)] # for gathering traits of individuals
                                           # of only species sp in sp-th row
    n = [[] for sp in range(S)] # for storing total abundance of each species
    for sp in range(S):
        for i in range(len(sppid)):
            if sppid[i]==sp: # if individual i belongs to species sp:
                trait_by_spp[sp].append(trait[i]) # add trait to sp-th row
        if len(trait_by_spp[sp]) > 0: # execute if species sp is not extinct
            n[sp], bins, patches = plt.hist(trait_by_spp[sp], 50, # histogram
                                            range=[-1, 1], color=cols[sp%7],
                                            histtype='stepfilled', alpha=0.3)
    bmax = numpy.array(n).max() # where to draw top of intrinsic growth function
    plt.plot([-1, -0.4999, -0.5, 0.4999, 0.5, 1], [0, 0, bmax, bmax, 0, 0],
             color='r', linestyle='dashed') # draw intrinsic growth function
    plt.xlabel('Trait value')
    plt.ylabel('Number of individuals')
    plt.show() # show plot

    
if (__name__ == '__main__'):
    status = main(sys.argv)
