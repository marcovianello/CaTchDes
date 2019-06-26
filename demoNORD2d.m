clear;

% NORD: 2D DEMO

% POLYGON: array of vertices coordinates   
% notice: for any polygon the vertices must be in counterclockwise or 
% clockwise order and the first vertex must be repeated at the end

% EXAMPLE: France-shaped polygon
Q=[0 0;8.2 -2;8 0;10 0.8;13 -1;15.5 1;14 4;15 4.8;14.5 7.2;13 6.5; ...
    15 10;16 10.3;17 14.5;12.5 15.7;7.2 19.6;6 19.2;6 17.5;3 16.1; ...
    3.2 15.7; 1.1 15.5; 1 16.1; 0 16.1;0.3 13.7;-4 13.8;-4.5 12; ...
    -1 10.5; 1 8;0 0];

% initializing regression degree, G-efficiency and max iterations
n=8; 
tol=0.95;
maxit=10000;

% minimal rectangle containing Q  
a=min(Q(:,1));b=max(Q(:,1));
c=min(Q(:,2));d=max(Q(:,2));
q1=(a+b)/2;  
q2=(c+d)/2;
r=max(b-a,c-d)/2;

% intersecting an mxm grid with Q 
m=100;
u=linspace(q1-r,q1+r,m);
v=linspace(q2-r,q2+r,m);
[xx,yy]=meshgrid(u,v);
p=[xx(:) yy(:)];
ind=find(inpolygon(p(:,1),p(:,2),Q(:,1),Q(:,2)));
pts=p(ind,:);
 
% computing a near optimal design with cputime  
tic
[cpts,cw,geff,momerr]=NORD(n,pts,tol,maxit);
toc
  
% displaying results
fprintf('\n initial design cardinality = %4.0f \n',length(pts(:,1)));
fprintf('\n concentrated support cardinality = %4.0f \n',length(cw));
fprintf('\n compression ratio = %4.0f \n',length(pts(:,1))/length(cw));
fprintf('\n g-efficiency = %1.4f \n',geff);
fprintf('\n moment reconstruction error = %4.2e \n \n',momerr);

% plotting the region, the grid and the concentrated design support 
plot(Q(:,1),Q(:,2),'b','MarkerSize',5);
axis square;
axis off;
hold on;
plot(p(:,1),p(:,2),'k.','MarkerSize',0.1);
plot(pts(:,1),pts(:,2),'b.','MarkerSize',0.1);
plot(cpts(:,1),cpts(:,2),'ro');
hold off;

