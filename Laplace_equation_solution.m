%%%%% Define the rectangle and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx=1; Ly=1; %rectangle dimensions
Nx=100; Ny=100; %# of intervals
nx=Nx+1; ny=Ny+1; %# of gridpoints in x,y directions including boundaries
dx=Lx/Nx; dy=Ly/Ny; %grid size in x,y directions
x=(0:Nx)*dx; y=(0:Ny)*dy; %x,y values on the grid
%%%%% Define the indices associated with the boundaries %%%%%%%%%%%%%%%%%%%
% boundary_index = [bottom, left, top, right]
boundary_index=[          1:nx,   1:nx:1+(ny-1)*nx, ...
             1+(ny-1)*nx:nx*ny,   nx:nx:nx*ny           ];  
%%%%% Set up matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagonals = [4*ones(nx*ny,1), -ones(nx*ny,4)];
A=spdiags(diagonals,[0 -1 1 -nx nx], nx*ny, nx*ny); %use sparse matrices
I=speye(nx*ny);
A(boundary_index,:)=I(boundary_index,:);
%%%%% SET-UP RIGHT HAND SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(nx,ny);
b(:,1)=0;          %bottom
b(1,:)=0;          %left
b(ny,:) = x.*(2 - x);   %top
b(:,nx) = y.*(2 - y);  %right
b=reshape(b,nx*ny,1); %make column vector
%%%%% Solve the Laplace equation using Gaussian elimination %%%%%%%%%%%%%%%
Phi=A\b; %solution step (all the computational time is here)
Phi=reshape(Phi,nx,ny); %make matrix
%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(x,y);
v = 0:0.1:1; 
contour(X,Y,Phi',v,'ShowText','on');%requires transpose (read the notes)
axis equal;
set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1]);
set(gca, 'XTick', [0 0.2 0.4 0.6 0.8 1]);
xlabel('$x$','Interpreter','latex','FontSize',14 );
ylabel('$y$','Interpreter','latex','FontSize',14);
title('Solution of the Laplace equation','Interpreter','latex','FontSize',16);
