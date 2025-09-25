r=28; sigma=10; beta=8/3; 
x1=0; y1=0; z1=0;
x2=sqrt(beta*(r-1)); y2=sqrt(beta*(r-1)); z2=r-1;
x3=-sqrt(beta*(r-1)); y3=-sqrt(beta*(r-1)); z3=r-1;
nx=500; nz=500;
xmin=-40; xmax=40; zmin=-40; zmax=40;
x_grid=linspace(xmin,xmax,nx); z_grid=linspace(zmin,zmax,nz);
[X,Z]=meshgrid(x_grid,z_grid);

y0 = 3*sqrt(2);

% Newton settings
tolF = 1e-10;        % residual tolerance
tolStep = 1e-10;     % step tolerance
maxIter = 40;        % a bit above 33 for safety

X_out = zeros(size(X));  

for ii = 1:nz          
    for jj = 1:nx
        % initial guess at this grid point
        x = X(ii,jj);
        y = y0;
        z = Z(ii,jj);

        it = 0;
        while it < maxIter
            % F and J for Lorenz fixed points
            f1 = sigma*(y - x);
            f2 = x*(r - z) - y;
            f3 = x*y - beta*z;
            F  = [f1; f2; f3];

            if max(abs(F)) < tolF, break; end

            J = [-sigma,   sigma,   0;
                  (r - z),   -1,   -x;
                   y,         x,   -beta];

            % Newton step and update
            delta = -J\F;
            x = x + delta(1);
            y = y + delta(2);
            z = z + delta(3);
            it = it + 1;

            if max(abs(delta)) < tolStep, break; end
        end

        % snap to the nearest of the three analytic roots (robust coloring)
        d1 = (x-x1)^2 + (y-y1)^2 + (z-z1)^2;
        d2 = (x-x2)^2 + (y-y2)^2 + (z-z2)^2;
        d3 = (x-x3)^2 + (y-y3)^2 + (z-z3)^2;
        if d1 <= d2 && d1 <= d3
            X_out(ii,jj) = x1;
        elseif d2 <= d3
            X_out(ii,jj) = x2;
        else
            X_out(ii,jj) = x3;
        end
    end
end

% overwrite X with the converged x-values for the provided coloring code
X = X_out;


eps=1.e-03;
X1 = abs(X-x1) < eps; X2 = abs(X-x2) < eps; X3 = abs(X-x3) < eps;
X4 = ~(X1+X2+X3);
figure; 
map = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; colormap(map); %[red;green;blue;black]
X=(X1+2*X2+3*X3+4*X4); 
image([xmin xmax], [zmin zmax], X); set(gca,'YDir','normal');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize',14);
title('Fractal from the Lorenz Equations', 'Interpreter', 'latex','FontSize', 16)  
