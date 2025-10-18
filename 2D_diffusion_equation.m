%%%%% Define the square and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=1;  %square is 2L x 2L 
N=100; %# of intervals in x and y directions
n=N+1; %# of gridpoints in x,y directions including boundaries
h=2*L/N;  %grid size in x,y directions
x=-L + (0:N)*h; %x values on the grid
y=-L + (0:N)*h; %y values on the grid
[X,Y]=meshgrid(x,y);
%%%%% Define the indices associated with the boundaries %%%%%%%%%%%%%%%%%%%
% boundary_index = [bottom, left, top, right]
boundary_index=[          1:n,       1:n:1+(n-1)*n, ...
                1+(n-1)*n:n*n,   n:n:n*n           ]; 
%%%%% Diffusion constant and time-step parameters
D=1;
dt=h^2/(2*D); %borderline stability of FTCS scheme
alpha=dt*D/h^2; %equation parameter
nsteps=1000; %number of time steps
%%%%% CONSTRUCT THE MATRIX AND COMPUTE LU DECOMPOSITION %%%%%%%%%%%%%%%%%%%%
% Build 2D Laplacian with natural ordering (x-index varies fastest)
NN = n*n;
e  = ones(n,1);
T  = spdiags([e -2*e e],[-1 0 1],n,n);  % 1D second-difference
I  = speye(n);
L2 = kron(I,T) + kron(T,I);             % 2D Laplacian (5-pt stencil)

A  = speye(NN) - (alpha/2)*L2;          % CN left matrix
B  = speye(NN) + (alpha/2)*L2;          % CN right matrix

% Enforce Dirichlet u=0 on the boundary: rows -> identity in A, zeros in B
bnd = unique(boundary_index(:));
A(bnd,:) = 0; 
A(sub2ind([NN,NN],bnd,bnd)) = 1;
B(bnd,:) = 0;

% LU factorization (P*A = L*U). Solve via U\(L\(P*b))
[LF,UF,Pperm] = lu(A);
%%%%% Define initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(n,n,nsteps);
sigma=L/4;
u(:,:,1)=1/(2*pi*sigma^2)*exp(-0.5*(X.^2+Y.^2)/sigma^2); 
u(1,:,1)=0; u(n,:,1)=0; u(:,1,1)=0; u(:,n,1)=0; %b.c.
%%%%% ADVANCE SOLUTION u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=2:nsteps
    uprev = u(:,:,m-1);
    rhs   = B * uprev(:);               % right-hand side vector
    y         = LF \ (Pperm * rhs);     % forward solve
    unew_vec  = UF \ y;                 % back solve
    unew = reshape(unew_vec,[n,n]);     % back to grid
    % re-enforce boundary for robustness
    unew(1,:) = 0; unew(n,:) = 0; unew(:,1) = 0; unew(:,n) = 0;
    u(:,:,m) = unew;
end
% %%%%% Plot with animation:  UNCOMMENT TO RUN ON MATLAB ONLINE OR DESKTOP %%%
figure('units','normalized','outerposition',[0 0 1 1])
s=surf(X,Y,u(:,:,1)); zlim([0, 2.6]);
xlabel('$x$','Interpreter','latex','FontSize',14); 
ylabel('$y$','Interpreter','latex','FontSize',14); 
zlabel('$u(x,y,t)$','Interpreter','latex','FontSize',14); 
title('Solution of the 2D diffusion equation','Interpreter','latex','FontSize',16);
pause(1)
for j=2:nsteps
     s.ZData=u(:,:,j); pause(0.01);
end
