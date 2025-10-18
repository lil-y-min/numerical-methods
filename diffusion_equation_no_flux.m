D=1; %diffusion coefficient
%%%%% Define the x-domain and x-grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx=1; %domain: -Lx < x < Lx
Nx=500; %# of intervals
nx=Nx+1;%# of gridpoints in x
dx=2*Lx/Nx; %grid length in x
x=-Lx + (0:Nx)*dx; %x values on the grid
%%%%% Time step parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsteps=10000; %number of time steps
nout=500; %plot every nout time steps
dt=(dx)^2/(2*D); %borderline stability of FTCS scheme
alpha=dt*D/dx^2; %equation parameter
%%%%% Construct the matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagonals = [2*(1+alpha)*ones(nx,1), -alpha*ones(nx,2)];
A=spdiags(diagonals,[0 -1 1], nx, nx);
A(1,:) = [1+alpha, -alpha, zeros(1,nx-2)];     %NO-FLUX BOUNDARY CONDITION
A(nx,:) = [zeros(1,nx-2), -alpha, 1+alpha];    %NO-FLUX BOUNDARY CONDITION
%%%%% Define initial conditions and plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma=Lx/16;
u=1/(sigma*sqrt(2*pi))*exp(-0.5*(x/sigma).^2); u=u';
plot(x,u,'r'); hold on;
xlabel('$x$','Interpreter','latex','FontSize',14); 
ylabel('$u(x, t)$','Interpreter','latex','FontSize',14); 
title('Solution of the diffusion equation','Interpreter','latex','FontSize',16);
%%%%% Advance solution and plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:nsteps
    b=[0; [alpha*u(1:nx-2) + 2*(1-alpha)*u(2:nx-1) + alpha*u(3:nx)]; 0];
    u=A\b;
    if mod(m,nout)==0, plot(x,u,'b'), end
end
