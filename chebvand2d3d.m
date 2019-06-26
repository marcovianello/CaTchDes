function V = chebvand2d3d(deg,gmesh,rect);

% computes by recurrence the Chebyshev-Vandermonde matrix on a 2d or 3d 
% arbitrarily located mesh, in the total-degree product Chebyshev basis 
% of a given rectangle with the graded lexicographical order  

% by Len Bos (Univ. of Verona) and Marco Vianello (Univ. of Padova)
% v1.0, July 2019


% INPUT:
% deg = polynomial degree
% gmesh = 2- or 3-column array of mesh point coordinates
% rect = 4- or 6-component vector such that the rectangle 
% [rect(1),rect(2)] x [rect(3),rect(4)] in 2d
% or [rect(1),rect(2)] x [rect(3),rect(4)] x [rect(5),rect(6)] in 3d
% contains the mesh


% OUTPUT:
% V = Chebyshev-Vandermonde matrix 

% FUNCTION BODY

if length(gmesh(1,:)) == 2 

% rectangle containing the mesh 
if isempty(rect) 
rect=[min(gmesh(:,1)) max(gmesh(:,1)) min(gmesh(:,2)) max(gmesh(:,2))];
end;
    
% couples with length less or equal to deg
% graded lexicographical order
j=(0:1:deg);
[j1,j2]=meshgrid(j);
dim=(deg+1)*(deg+2)/2;
couples=zeros(dim,2);
for s=0:deg
good=find(j1(:)+j2(:)==s);
couples(1+s*(s+1)/2:(s+1)*(s+2)/2,:)=[j1(good) j2(good)];
end;

% mapping the mesh in the square [-1,1]^2
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*gmesh(:,1)-b-a)/(b-a) (2*gmesh(:,2)-d-c)/(d-c)];

% Chebyshev-Vandermonde matrix on the mesh
T1=chebpolys(deg,map(:,1));
T2=chebpolys(deg,map(:,2));
V=T1(:,couples(:,1)+1).*T2(:,couples(:,2)+1);

end;


if length(gmesh(1,:)) == 3 
        
% parallelepiped containing the mesh  
if isempty(rect)
rect=[min(gmesh(:,1)) max(gmesh(:,1)) min(gmesh(:,2)) max(gmesh(:,2)) ... 
 min(gmesh(:,3)) max(gmesh(:,3))];
end;
    
% triples with length less or equal to deg
% graded lexicographical order
j=(0:1:deg);
[j1,j2,j3]=meshgrid(j,j,j);
dim=(deg+1)*(deg+2)*(deg+3)/6;
triples=zeros(dim,3);
for s=0:deg
good=find(j1(:)+j2(:)+j3(:)==s);
triples(1+s*(s+1)*(s+2)/6:(s+1)*(s+2)*(s+3)/6,:)=...
    [j1(good) j2(good) j3(good)];
end;

%mapping the mesh in the cube [-1,1]^3
a=rect(1);b=rect(2);
c=rect(3);d=rect(4);
e=rect(5);f=rect(6);
map=[(2*gmesh(:,1)-b-a)/(b-a) (2*gmesh(:,2)-d-c)/(d-c) ... 
(2*gmesh(:,3)-f-e)/(f-e)];

%Chebyshev-Vandermonde matrix on the mesh
T1=chebpolys(deg,map(:,1));
T2=chebpolys(deg,map(:,2));
T3=chebpolys(deg,map(:,3));
V=T1(:,triples(:,1)+1).*T2(:,triples(:,2)+1).*T3(:,triples(:,3)+1);

end; 

end


function T=chebpolys(deg,x)

% computes the Chebyshev-Vandermonde matrix on the real line by recurrence

% INPUT:
% deg = maximum polynomial degree
% x = 1-column array of abscissas

% OUTPUT
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg

% inizialization
T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

% 3-term recurrence 
for j=2:deg
t2=2*x.*t1-t0;
T(:,j+1)=t2;
t0=t1;
t1=t2;
end;

end

