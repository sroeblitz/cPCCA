function varargout = problem_pcca_nlscon(x,flag,par)
%
%     Problem function.
%
%     VAR = PROBLEM(X,FLAG,par) computes the problem function f(X),
%       Jacobian df/dx(X),startvector for the Gauss-Newton
%       iteration and vector of to be fitted values.
%
% FX = PROBLEM(X,'',par)          
%    returns the right hand side f(X) for a input column vector X.
% JAC = PROBLEM(X,'jacobian',t) 
%    returns the Jacobian df/dx(X) for a input column vector X.


 switch flag
 case ''                                 % Return y = f(x).
   [varargout{1:2}] = f_rhs(x,par) ;
 case 'jacobian'                         % Return Jacobian matrix df/dx.
   [varargout{1:2}] = jacobian(x,par);
 otherwise
   error(['Unknown flag ''' flag '''.']) ;
 end

% --------------------------------------------------------------------------
 

function [f,ifail] = f_rhs(x,par)

ifail=1;

EVS=par.evs;
k=size(EVS,2);

A=zeros(k,k);

for i=1:k-1
    A(i+1,2:k)=x(((i-1)*(k-1)+1):i*(k-1));
end

A=fillA(A, EVS);

F=trace(diag(1./A(1,:))*(A'*A));

f=k-F;

ifail = 0;


% --------------------------------------------------------------------------

function [JF,ifail] = jacobian(x2,par)

EVS=par.evs;
k=size(EVS,2);

  
ifail=1;

A=zeros(k,k);


for i=1:k-1
    A(i+1,2:k)=x2(((i-1)*(k-1)+1):i*(k-1));
end

A=fillA(A,EVS);

JF=zeros(k,k);
for i=1:k
    for j=1:k
        JF(i,j)=A(i,j)/A(1,j);
    end
end

JF(1,:)=[];
JF(:,1)=[];

JF=JF(:);
JF=JF';

ifail = 0;

