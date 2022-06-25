function A=fillA(A,EVS)
% make A feasible 
%
% A=fillA(A,EVS)
%
% Input:
%   A:      (infeasible) (k,k)-matrix
%   EVS:    (N,k)-matrix with eigenvectors of A for eigenvalues close to one
%
% Output:
%   A:      feasible (k,k)-matrix A
%
% Written by Susanna Roeblitz and Marcus Weber, Zuse Institute Berlin, 2012

[N,k]=size(EVS);

% compute 1st column of A by row sum condition 
A(2:k,1)=-sum(A(2:k,2:k),2);

% compute 1st row of A by maximum condition
for j=1:k
    A(1,j)=- EVS(1,2:k)*A(2:k,j);
    for l=2:N
        dummy = - EVS(l,2:k)*A(2:k,j);
        if (dummy > A(1,j))
            A(1,j) = dummy;
        end
    end
end

% reskale A to be in the feasible set
A=A/sum(A(1,:));

  