function EVS=orthogon(EVS,pi)
% Eigenvectors computed by an eigenvalue solver are normed to 1 w.r.t.
% the Euclidean metric. Theory requires pi-orthonormal eigenvectors,
% i.e. EVS'*diag(pi)*EVS=0.
% Moreover, the 1st eigenvectors is not constant in general. This is 
% readjusted here.
%
% EVS=orthogon(EVS,pi,lambda)
%
% Input:
%   EVS:    (N,k)-matrix with eigenvectors columnwise
%   pi:     N-vector with stationary density    
%
% Output:
%   EVS:    (N,k)-matrix with orthogonal eigenvectors columnwise
%
% written by Marcus Weber, Zuse Institute Berlin, Takustrasse 7, 14195 Berlin


[N,k]=size(EVS);
EVS_old=EVS;
%perron=1;

% scale eigenvectors to correct length
pi=pi/sum(pi);
for i=1:k
    EVS(:,i)=EVS(:,i)/(sqrt((EVS(:,i).*pi)'*EVS(:,i))*sign(EVS(1,i)+eps));
end

%subspace(EVS_old,EVS)
% % count possible degenerations
% for i=1:k
%     if lambda(i) > 0.9999
%       perron=perron+1;
%     else
%       break;
%     end
% end

%perron=min(perron,size(EVS,2));
perron=size(EVS,2);

if (perron > 1) % search for constant eigenvector
    maxscal=0.0;
    for i=1:perron
        scal=sum(pi.*EVS(:,i));
        if (abs(scal) > maxscal)
            maxscal = abs(scal);
            maxi=i;
        end
    end
    % shift non-constant EV to the right 
    EVS(:,maxi)=EVS(:,1);
    EVS(:,1)=ones(N,1);
    %subspace(EVS_old,EVS)
    %EVS(find(pi<=0),:)=0;
    % pi-orthogonalize all other EVs
    for i=2:k
        for j=1:i-1
            scal = (EVS(:,j).*pi)'*EVS(:,i);
            EVS(:,i) = EVS(:,i)- scal * EVS(:,j);
        end
        sumval =sqrt((EVS(:,i).*pi)'*EVS(:,i));
        EVS(:,i)=EVS(:,i)/sumval;
    end
    %subspace(EVS_old,EVS)
end

