%%%%% Define the rectangle and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx=1; Ly=1; %rectangle dimensions
Nx=100; Ny=100; %# of intervals
nx=Nx+1; ny=Ny+1; %# of gridpoints in x,y directions including boundaries
dx=Lx/Nx; dy=Ly/Ny; %grid size in x,y directions
x=(0:Nx)*dx; y=(0:Ny)*dy; %x,y values on the grid
%%%%% Define the iteration parameters and initial condition %%%%%%%%%%%%%%%
eps=1.e-6; %convergence criteria for each value of Phi
index_x=2:nx-1;  index_y=2:ny-1; %internal grid points 
Phi=zeros(nx,ny);%matrix with solution and boundary conditions
%%%%% Define the boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the boundary conditions
Phi(:,1)=0;         %bottom
Phi(1,:)=0;         %left
Phi(:,ny)=x.*(2-x);        %top
Phi(nx,:)=y.*(2-y);        %right
%%%%% Jacobi iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_old=Phi;
error=2*eps; ncount=0;
while (error > eps)
   ncount=ncount+1;
   Phi(index_x,index_y)=0.25*(Phi(index_x+1,index_y) ...
    +Phi(index_x-1,index_y)+Phi(index_x,index_y+1)+Phi(index_x,index_y-1));
   error=max(abs(Phi(:)-Phi_old(:)));
   if any(isnan(Phi(:))) || any(isinf(Phi(:)))
        fprintf('iterations diverge\n');
        return;
    end
    Phi_old=Phi;
    %fprintf('%g %e\n',ncount, error);
end
fprintf('%g\n',ncount);   
%%%%% graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(x,y);
v=0:0.1:1; %SET THE CONTOUR LEVELS
contour(X,Y,Phi',v,'ShowText','on');%requires transpose (read the notes)
axis equal;
set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1]);
set(gca, 'XTick', [0 0.2 0.4 0.6 0.8 1]);
xlabel('$x$','Interpreter','latex','FontSize',14 );
ylabel('$y$','Interpreter','latex','FontSize',14);
title('Solution of the Laplace equation','Interpreter','latex','FontSize',16);
