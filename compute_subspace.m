function [EVS,lambda,kmax]=compute_subspace(M,kmax,target,opts)
%compute invariant subspace of dimension kmax for matrix M corresponding to
%eigenvalues closest to target

% Compute eigenvectors and sort them according to distance from 1
[EVS,la]=eigs(M,kmax,target,opts);
la=diag(la)
%evsmod=preprocessEVS(E,la);

% la=abs(diag(la)); %Same as in circular FPCCA+
[~,index]=sort(abs(target-la),'ascend');
EVS=EVS(:,index);
lambda=la(index);

% extract real and imaginary parts of eigenvectors
i=1;
while (i<size(EVS,2))
    if (imag(lambda(i))~=0)
        EVS(:,i)=real(EVS(:,i));
        EVS(:,i+1)=imag(EVS(:,i+1));
        lambda(i)=real(lambda(i));
        lambda(i+1)=imag(lambda(i+1));
        i=i+2;
    else
        i=i+1;
    end
end
if imag(EVS(:,end))~=0
    %disp('reduce kmax by 1; otherwise split of complex eigenpair');
    kmax=kmax-1;
    lambda=lambda(1:end-1)
    EVS(:,end)=[];
end