clear;

% NORD: 3D DEMO

% CUBE [-1,1]^3 
% mxmxm grid
m=30;
u=linspace(-1,1,m);
[a,b,c]=meshgrid(u);
pts=[a(:) b(:) c(:)];

% initializing regression degree, G-efficiency and max iterations
n=4; 
tol=0.95;
maxit=10000;

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

% plotting the cube edges, the grid and the concentrated design support
plot3(cpts(:,1),cpts(:,2),cpts(:,3),'ro','MarkerSize',3);
axis square;
hold on;
plot3(pts(:,1),pts(:,2),pts(:,3),'g.','MarkerSize',0.01);
vert=[-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 -1];
plot3(vert(:,1),vert(:,2),vert(:,3),'k-','MarkerSize',2);
vert=[-1 -1 1;1 -1 1;1 1 1;-1 1 1;-1 -1 1];
plot3(vert(:,1),vert(:,2),vert(:,3),'k-','MarkerSize',2);
vv1=[-1 -1 -1;-1 -1 1];vv2=[1 -1 -1;1 -1 1];
vv3=[1 1 -1;1 1 1];vv4=[-1 1 -1;-1 1 1];
plot3(vv1(:,1),vv1(:,2),vv1(:,3),'k-','MarkerSize',2);
plot3(vv2(:,1),vv2(:,2),vv2(:,3),'k-','MarkerSize',2);
plot3(vv3(:,1),vv3(:,2),vv3(:,3),'k-','MarkerSize',2);
plot3(vv4(:,1),vv4(:,2),vv4(:,3),'k-','MarkerSize',2);
hold off;

