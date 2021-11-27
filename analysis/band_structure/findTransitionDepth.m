function freq=findTransitionDepth(depth,n1,n2,k)
%FINDLATTICEDEPTH Summary of this function goes here
%   Detailed explanation goes here
% depths=linspace(1,20,100);
% paramsBANDS;

switch nargin
    case 1
        n1=1;
        n2=3;
        k=0;
end   

nfo = struct;
nfo.depth = depth;
nfo.k = k;
nfo.numStates = 31;

H1=makeH(nfo);
[~,b1]=eig(full(H1));
lambda1=b1*ones(nfo.numStates,1); 
eng=abs(lambda1(n1)-lambda1(n2));

freq=eng;
end

function Hmatrix=makeH(nfo)
% Calculates the hamiltonian.
depth=      nfo.depth;                    % static lattice depth
k=          nfo.k;                    % quasimomentum
numStates=  nfo.numStates;    % number of basis states

%%%%%%%%%%%%%% crystal momentum operator %%%%%%%%%%%%%%%%
C=0;
for ii=2:2:numStates
   C=[C; ii; 0]; 
end
C(end)=[];
D=zeros(numStates,1);
M1=gallery('tridiag',C,D,-C);

%%%%%%%%%%%%% kinetic energy operator %%%%%%%%%%%%%%
M2=0;
for ii=2:2:(numStates-1)
    M2=[M2 ii^2 ii^2];
end
M2=sparse(-diag(M2));

%%%%%%%%%%%%%%%%%% lattice operator %%%%%%%%%%%%%%%%%%%%
rowsL=[3:1:numStates];
colsL=[1:1:numStates-2];
valsL=[sqrt(2) ones(1,numStates-3)];
M3=sparse(rowsL,colsL,valsL,numStates,numStates)+...
   sparse(colsL,rowsL,valsL,numStates,numStates);

%%%%%%%%%%%%%% Put it all together %%%%%%%%%%%%%%%%%
Hmatrix=(k^2*speye(numStates)-2*1i*k*M1-M2)...
    -depth/4*M3;
end
