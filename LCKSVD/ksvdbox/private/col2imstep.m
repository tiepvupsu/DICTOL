%COL2IMSTEP Rearrange matrix columns into blocks.
%  A = COL2IMSTEP(B,[MM NN],[N1 N2]) rearranges the columns of B into
%  sliding N1-by-N2 blocks producing the matrix A of size MM-by-NN. B is
%  usually the result of calling IM2COLSTEP(...) with a stepsize of 1, or
%  using Matlab's IM2COL(..,'sliding'). Overlapping blocks are summed in A.
%
%  A = COL2IMSTEP(B,[MM NN],[N1 N2],[S1 S2]) arranges the blocks in A with
%  a step size of (S1,S2) between them. The first block is at A(1:N1,1:N2),
%  and the rest are at A((1:N1)+i*S1,(1:N2)+j*S2). Overlapping blocks are
%  summed in A. Note that B is usually the result of calling
%  IM2COLSTEP(...) with a stepsize of [S1 S2].
%
%  A = IM2COLSTEP(B,[MM NN KK],[N1 N2 N3],[S1 S2 S3]) generates a 3-D
%  output matrix A. The step size [S1 S2 S3] may be ommitted, and defaults
%  to [1 1 1].
%
%  See also IM2COLSTEP, IM2COL, COUNTCOVER.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
