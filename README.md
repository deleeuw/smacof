# smacof


This is a code repository for the smacof project -- the many variants of smacof for metric, constrained, non-metric, individual differences are developed in C, using a uniform programming strategy for each instance. 

The programs are all called smacofXXXX(), where XXXX is a four letter
code. 

* The first X is either S for single matrix or M for multiple matrices.
* The second X is either S for symmetric or A for asymmetric.
* The third X is either U for unweighted or W for weighted.
* The fourth X is either R for ratio, I for interval, P for polynomial, 
  O for ordinal, S for splinical, Q for ordinal polynomial, T for ordinal
  splinical.
  
That makes 56 programs. They share some code, which will eventually all
be in the common directory, and they share a common include file 
smacof.h. For now all smacofXXX.c file have a main() as a driver
for the basoc subroutines, applied to an example.

The code accompanies the book "Least Squares Euclidean Multidimensional Scaling".

The book is still very incomplete and sketchy -- but you can look at updates on my web page https://jansweb.netlify.app. 
