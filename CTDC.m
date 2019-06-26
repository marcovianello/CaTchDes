function [pts,w,momerr] = CTDC(deg,X,omega)

% Caratheodory-Tchakaloff discrete measure concentration in 2d and 3d, 
% for example probability measures (designs) or quadrature formulas;
% the moments are invariant (close to machine precision) up to degree deg;  
% adapts to the (numerical) dim of the polynomial space on the support; 
% works satisfactorily for low/moderate degrees

% by Len Bos (Univ. of Verona) and Marco Vianello (Univ. of Padova)
% v1.0, July 2019

% INPUT:   
% deg: polynomial exactness degree
% X: 2 or 3 column array of mass points coordinates
% omega: 1-column array of weights 

% OUTPUT:
% pts: 2 or 3 column array of extracted mass points
% w: 1-column array of corresponding new weights 
% momerr: moment reconstruction error
 

% FUNCTION BODY 

% total-degree Chebyshev-Vandermonde matrix at X
U=chebvand2d3d(deg,X,[]);

% dimension of the polynomial space on the support 
N=rank(U);

% scaling the matrix rows by the sqrt of the weights 
for k=1:length(U(1,:))
B(:,k)=U(:,k).*sqrt(omega);
end;
% polynomial basis orthogonalization
if N<length(U(1,:))
[Q0,R0,pm]=qr(B,'vector');
p=pm(1:N);R=R0(1:N,1:N);
Q=U(:,p)/R; 
else 
[Q0,R]=qr(B,0);
Q=U/R;
end;

% moments of the orthogonal basis
orthmom=Q'*omega;

% Caratheodory-Tchakaloff points and weights via NNLS 
weights=lsqnonneg(Q',orthmom);
% indexes of nonvanishing weights and compression  
ind=find(abs(weights)>0);
pts=X(ind,:);
w=weights(ind);

% moment reconstruction error
momerr=norm(U(ind,:)'*w-U'*omega);

end


