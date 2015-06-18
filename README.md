# diagonals: Fat Diagonals in R
Bastiaan Quast  
June 18, 2015  

# Abstract
We present the `diagonals` R package, which implements functions for dealing with **fat diagonals**.
Fat diagonals are block matrix-diagonals that occur when two or more dimensions are mapped along a single edge of a matrix.
An asymetric network graph (e.g. a dyadic social network) in a matrix, we would need each node along each edge of matrix,
however we would also need to map the direction of the tie, which is an additional dimension.
This could be presented as such (where the `i` and `o` suffices prepresent incoming and outgoing repectively). 


```
##   Ai Ao Bi Bo Ci Co Di Do
## A  1  1  0  0  1  0  1  0
## B  0  1  1  1  0  1  0  1
## C  1  0  0  0  1  1  0  0
## D  1  0  1  0  1  1  1  1
```

