from .libsmop import *


@function
def randsample(n=None, k=None, replace=None, w=None):
    # RANDSAMPLE Random sample, with or without replacement.
    #   Y = RANDSAMPLE(N,K) returns Y as a 1-by-K vector of values sampled
    #   uniformly at random, without replacement, from the integers 1:N.
    #
    #   Y = RANDSAMPLE(POPULATION,K) returns K values sampled uniformly at
    #   random, without replacement, from the values in the vector POPULATION.
    #
    #   Y = RANDSAMPLE(...,REPLACE) returns a sample taken with replacement if
    #   REPLACE is true, or without replacement if REPLACE is false (the default).
    #
    #   Y = RANDSAMPLE(...,true,W) returns a weighted sample, using positive
    #   weights W, taken with replacement.  W is often a vector of probabilities.
    #   This function does not support weighted sampling without replacement.
    #
    #   Example:  Generate a random sequence of the characters ACGT, with
    #   replacement, according to specified probabilities.
    #
    #      R = randsample('ACGT',48,true,[0.15 0.35 0.35 0.15])
    #
    #   See also RAND, RANDPERM.
    #
    #   Copyright 1993-2004 The MathWorks, Inc.
    #   $Revision: 1.1.4.1 $  $Date: 2003/11/01 04:28:51 $

    nargin = randsample.nargin

    if nargin < 2:
        error('stats:randsample:TooFewInputs', 'Requires two input arguments.')
    else:
        if numel(n) == 1:
            population = copy([])
        else:
            population = n
            n = numel(population)
            if length(population) != n:
                error('stats:randsample:BadPopulation', 'POPULATION must be a vector.')

    if nargin < 3:
        replace = false

    if nargin < 4:
        w = copy([])
    else:
        if logical_not(isempty(w)):
            if length(w) != n:
                if isempty(population):
                    error('stats:randsample:InputSizeMismatch', 'W must have length equal to N.')
                else:
                    error('stats:randsample:InputSizeMismatch', 'W must have the same length as the population.')
            else:
                p = ravel(w).T / sum(w)

    if cellarray([true, 'true', 1]) == replace:
        if isempty(w):
            y = ceil(multiply(n, rand(k, 1)))
        else:
            dum, y = histc(rand(k, 1), concat([0, cumsum(p)]), nargout=2)
        # Sample without replacement
    else:
        if cellarray([false, 'false', 0]) == replace:
            if k > n:
                if isempty(population):
                    error('stats:randsample:SampleTooLarge',
                          'K must be less than or equal to N for sampling without replacement.')
                else:
                    error('stats:randsample:SampleTooLarge', 'K must be less than or equal to the population size.')
            if isempty(w):
                # If the sample is a sizeable fraction of the population,
                # just randomize the whole population (which involves a full
                # sort of n random values), and take the first k.
                if dot(4, k) > n:
                    rp = randperm(n)
                    y = rp(arange(1, k))
                    # is wasteful.  Repeatedly sample with replacement until there are
                # k unique values.
                else:
                    x = zeros(1, n)
                    sumx = 0
                    while sumx < k:
                        x[ceil(dot(n, rand(1, k - sumx)))] = 1
                        sumx = sum(x)

                    y = find(x > 0)
                    y = y(randperm(k))
            else:
                error('stats:randsample:NoWeighting', 'Weighted sampling without replacement is not supported.')
        else:
            error('stats:randsample:BadReplaceValue', 'REPLACE must be either true or false.')

    if logical_not(isempty(population)):
        y = population(y)
    else:
        y = ravel(y)
    return y
