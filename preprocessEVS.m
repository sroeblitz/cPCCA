function EVS=preprocessEVS(EVS,lambda)
% separates real and imaginary parts of complex eigenvectors

i=1;
while (i<size(EVS,2))
    if (imag(lambda(i))~=0)
        EVS(:,i)=real(EVS(:,i));
        EVS(:,i+1)=imag(EVS(:,i+1));
        i=i+2;
    else
        i=i+1;
    end
end
if imag(EVS(:,end))~=0
    disp('reduce kmax by 1; otherwise split of complex eigenpair');
    kmax=kmax-1;
    lambda=lambda(1:end-1);
    EVS(:,end)=[];
end


end
