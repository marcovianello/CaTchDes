function [cpts,cw,geff,momerr] = NORD(deg,X,gefftol,maxit); 

% computes a Near Optimal Regression Design with a given G-efficiency on   
% a discrete set X in 2d or 3d, by a basic multiplicative algorithm and 
% Caratheodory-Tchakaloff discrete measure concentration 

% by Len Bos (Univ. of Verona) and Marco Vianello (Univ. of Padova)
% v1.0, July 2019

% INPUT
% deg: polynomial regression degree
% X: discrete set (2- or 3-column array)
% gefftol: G-efficiency tolerance 
% maxit: maximum number of iterations of the multiplicative algorithm 

% OUTPUT
% cpts,cw: compressed near optimal design points and weights 
% Y,w: near optimal design points and weights before compression
% geff: G-efficiency of the output design 
% momerr: moment reconstruction error by the concentrated measure

% FUNCTION BODY 

% initializing the probability weights
M=length(X(:,1));
w=ones(M,1)/M;

% Chebyshev-Vandermonde matrix at X for degree deg
VX=chebvand2d3d(deg,X,[]);

% rank computation and polynomial basis determination
N=rank(VX);
if N<length(VX(1,:))
[Q0,R0,pm]=qr(VX,'vector');
p=pm(1:N);
VX=VX(:,p);
end;

% initializing nodes and weights tolerance 
Y=X; VY=VX; wtol=1e-12;
nit=0; go=1;

% multiplicative iteration up to the given G-efficiency 
while go==1
nit=nit+1;
% eliminating tiny weights and corresponding nodes 
ind=find(w>wtol/length(w));
w=w(ind)/sum(w(ind));
VY=VY(ind,:);
Y=Y(ind,:);
% scaling the matrix rows by the sqrt of the weights 
for k=1:length(VY(1,:))
A(:,k)=VY(:,k).*sqrt(w);
end;
[Q,R]=qr(A,0);
% scaling the Christoffel polynomial at Y by the weights 
K=sum((Q.*Q)');
% computing G-efficiency of the weights on X 
U=VX/R;
geff=N/max(sum((U.*U)'));
% updating the cicle exit condition
go=(geff<gefftol & nit<maxit);
if go==1
% updating the weights
w=K'/N;
end;
clear A;

end;

% Caratheodory-Tchakaloff design concentration 
[cpts,cw,momerr]=CTDC(2*deg,Y,w);

end


