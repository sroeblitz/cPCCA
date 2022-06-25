function index=indexsearch(evs)
% find a simplex structure in the data 
%
% index=indexsearch(evs)
%
% Input:
%   evs:    (N,k)-matrix with eigenvectors columnwise
%
% Output:
%   index:  k-vector with indices of objects that build the simplex vertices
%
% written by Marcus Weber, Zuse Institute Berlin, Takustrasse 7, 14195 Berlin

maxdist=0.0;

[N,k]=size(evs);

OrthoSys=evs;
index=zeros(1,k);

% first vertex: row with largest norm
for l=1:N
    dist = norm(OrthoSys(l,:));
    if (dist > maxdist)
        maxdist=dist;
	    index(1)=l;
    end
end


OrthoSys=OrthoSys-ones(N,1)*evs(index(1),:);

% all further vertices as rows with maximum distance to existing subspace
for j=2:k
    maxdist=0.0;
    temp=OrthoSys(index(j-1),:);
    for l=1:N
        sclprod=OrthoSys(l,:)*temp';
        OrthoSys(l,:)=OrthoSys(l,:)-sclprod*temp;
        distt=norm(OrthoSys(l,:));
	    if distt > maxdist  %&& ~ismember(l,index(1:j-1))
            maxdist=distt;
            index(j)=l;
        end
    end
    OrthoSys = OrthoSys/maxdist;
end


   