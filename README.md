Abstract
========

We present the `diagonals` R package, which implements functions for
dealing with **fat diagonals**. Fat diagonals are block matrix-diagonals
that occur when two or more dimensions are mapped along a single edge of
a matrix. For an asymetric network graph (e.g. a dyadic social network)
to be mapped to a matrix, we would need each node along each edge of
matrix, however we would also need to map the direction of the tie,
which is an additional dimension. Typically these would be represented
as higher-order arrays (i.e. \>= 3). In order to effectively visualise
such arrays, it can be helpful to do so in a matrix. This could for
instance be represented as in the following matrix (where the `i` and
`o` suffices prepresent incoming and outgoing repectively).

    ##   Ai Ao Bi Bo Ci Co Di Do
    ## A  1  1  0  0  1  0  1  0
    ## B  0  1  1  1  0  1  0  1
    ## C  1  0  0  0  1  1  0  0
    ## D  1  0  1  0  1  1  1  1

Sometimes the ties of a node to itself are not particularly meaningful
(e.g. feeling of amiability towards oneself) and can be removed. For a
symetric network this can simply be done using the function `diag()` in
R's `base` package, e.g.

    sm <- matrix(1, nrow=4, ncol=4)
    diag(sm) <- NA
    sm

    ##      [,1] [,2] [,3] [,4]
    ## [1,]   NA    1    1    1
    ## [2,]    1   NA    1    1
    ## [3,]    1    1   NA    1
    ## [4,]    1    1    1   NA

However, for higher-order matrices this does not work well.

    diag(m) <- NA
    m

    ##   Ai Ao Bi Bo Ci Co Di Do
    ## A NA  1  0  0  1  0  1  0
    ## B  0 NA  1  1  0  1  0  1
    ## C  1  0 NA  0  1  1  0  0
    ## D  1  0  1 NA  1  1  1  1

In comes the `diagonals` package and its workhorse `fatdiag()` function.
The function is designed to mimmick the behaviour of the `diag()` as
closely as possible, but with then for **fat diagonals**.

    library(diagonals)

    ## 
    ## D I
    ## A G
    ##     O N
    ##     A L
    ##         S

    fatdiag(m, steps=4) <- NA
    m

    ##   Ai Ao Bi Bo Ci Co Di Do
    ## A NA NA  0  0  1  0  1  0
    ## B  0  1 NA NA  0  1  0  1
    ## C  1  0  0  0 NA NA  0  0
    ## D  1  0  1  0  1  1 NA NA
