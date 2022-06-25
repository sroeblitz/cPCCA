function optval=objective(alpha,EVS)
% compute objective function value
%
% optval=objective(alpha,EVS)
%
% Input:
%   alpha:  (k-1)^2-vector with current iterate
%           (it columnwise contains A(2:k,2:k))
%   EVS:    (N,k)-matrix with eigenvectors column by column
%
% Output:
%    optval:    current value of the objective function k-trace(S)
%
% written by Susanna Roeblitz and Marcus, Zuse Institute Berlin, Takustrasse 7, 14195 Berlin

k=size(EVS,2);
A=zeros(k,k);

% rebuild transformation matrix A
for i=1:k-1
    for j=1:k-1
      A(i+1,j+1)=alpha(j + (i-1)*(k-1));
    end
end


%make A feasible
A=fillA(A, EVS);

%compute value of objective function
optval=k - trace(diag(1./A(1,:))*(A'*A));  

% Note: other choices are possible here, e.g.:
% optval=-trace(log((diag(1./A(1,:))*(A'*A))));  %[White/Shalloway, 2009]
    

