#!/usr/local/bin/sage
# vim: set fileencoding=utf-8 :

## suppose we estimate a branching process with two types of people
## so there are four transmission probabilities p_12, etc.
## Or rather expected numbers of offspring m_12 etc.
## Given a current infective population (v_1, v_2)
## The expected next one is (m_11 v_1 + m_12 v_2, m_21 v_1 + m_22 v_2)
## or M v.
## Over the full course of an outbreak the final numbers of 1 and 2
## should be given by (1-M)^{-1} v.
## Given an index case v = either (1, 0) or (0, 1), leading to
## final numbers (c_1, c_2), we have some information about (1-M)^-1,
## but not enough to reconstruct M.
## Given an ensemble of pairs v, c, we have too much information to
## determine M but can get a best-fit M by basically linear regression.

## Transpose for comparison to standard linear regression:
##  C = V β
## Where each row of C is an outcome, each row of V is an index case,
## and β is (1-M)^-T.
## The least squares estimator then is β = (V^T V)^-1 V^T C.
## So let's try that.  V^T V won't have an inverse,
## but let's see what we can do.

## some test data
## 2 index cases
V = matrix( [
    [ 1, 0 ],
    [ 0, 1 ]
] )

## 2 final size vectors
C = matrix( [
    [ 5, 3 ],
    [ 4, 5 ]
] )

## to estimate M matrix
def estimate_M( V, C ):
    ## estimate β
    beta = (V.transpose() * V).inverse() * V.transpose() * C
    print(beta)
    print('')
    ## now if β = (1 - M)^-1, what is M
    ## M = 1 - β^-1
    M = identity_matrix(2) - beta.inverse()
    return(M)

print(estimate_M(V,C))

## try with more clusters
V = matrix( [
    [ 1, 0 ],
    [ 1, 0 ],
    [ 0, 1 ]
] )

C = matrix( [
    [ 5, 3 ],
    [ 6, 3 ],
    [ 4, 8 ]
] )

print('\n====')
print(estimate_M(V,C))

## what can we do if we have multiple index cases in one outbreak,
## so that we don't know all the rows of C but only their sum?
## We could average beta over all the combinatoric possibilities of C,
## I guess. Seems like we ought to be able to do better.
## Well maybe in that case we have a sum of V rows, and 1 dimension
## of info about M.
## What does that do?

V = matrix( [ [ 30, 25 ] ] )

C = matrix( [ [ 80, 51 ] ] )

print('\n====')
print(estimate_M(V,C))

