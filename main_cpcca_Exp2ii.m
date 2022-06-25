%% Main cPCCA+ run file for Example 2, case (ii)
close all; clear all;

tic

%Load data: Choose between Case(i): set d=1 or Case(ii): set d=2
load('Exp2_x09.mat')
n=size(P,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Coarse-graining and clustering via PCCA+
%density for normalization of eigenvectors
ovec=ones(n,1)/n; % (!)

% parameters for PCCA+
target=1+eps;           % target eigenvalue       
kmin=3;                 % minimum number of clusters
kmax=3;                 % maximum number of clusters
solver2='gauss-newton'; % optimization method for optimal cluster number
opts.disp=0;            % no display of eigenvalue solver

% the following two parameters are not used if kmin==kmax
flag=1;                 % flag=1: full optimization, flag=0: unscaled initial guess
solver1='nelder-mead';  % optmization method for finding optimal cluster number

% Calculate the eigenvalues:
[E,L]=eigs(P,3,target);
la=diag(L)
%separate real and imaginary parts
evsmod=preprocessEVS(E,la);

%PCCA+ routine
[chi,A,optval,EVS]=pcca(evsmod,ovec,kmin,kmax,flag,solver1,solver2);

%% Sort chi vectors 
[~,m]=max(chi);
[i,idx]=sort(m);
chi=chi(:,idx);

fprintf('\nList of states and their cluster with highest membership:\n')
for k=1:n
    [~,idx]=max(chi(k,:));
    fprintf('%d: %d\n',k,idx)
end

disp (' ')
disp ('Coarse-grainded transition matrix:') 
Pc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*P*chi  

%plot membership functions
figure(1)
plot(chi,'-s','LineWidth',4,'Markersize',20)
set(gca,'FontSize',20)
%title('Membership functions chi','FontSize',20)
xlabel('States','Fontsize',20)
ylabel('Membership to cluster in range [0  1]','FontSize',20)
axis([1 9 0 1])
xticks([1 2 3 4 5 6 7 8 9])
xticklabels({'1','2','3','4','5','6','7','8','9'})

%End elapse time
toc
et=toc;




