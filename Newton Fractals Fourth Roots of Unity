f = @(Z) Z.^4 - 1;
fp = @(Z) 4 * Z.^3;
root1 = 1; 
root2 = 1i;
root3 = -1;
root4 = -1i;

nx=1000; ny=1000;
xmin=-2; xmax=2; ymin=-2; ymax=2;

x=linspace(xmin,xmax,nx); y=linspace(ymin,ymax,ny);
[X,Y]=meshgrid(x,y);
Z=X+1i*Y;

nit=50;
for n=1:nit
    Z = Z - f(Z) ./ fp(Z);
end

eps=0.001;
Z1 = abs(Z-root1) < eps; Z2 = abs(Z-root2) < eps;
Z3 = abs(Z-root3) < eps; Z4 = abs(Z-root4) < eps;
Z5 = ~(Z1+Z2+Z3+Z4);

figure;
map = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 0 0]; colormap(map); %[red;blue;green;yellow;black]
Z=(Z1+2*Z2+3*Z3+4*Z4+5*Z5); 
image([xmin xmax], [ymin ymax], Z); set(gca,'YDir','normal');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize',14);
title('Fractal from $f(z)=z^4-1$', 'Interpreter', 'latex','FontSize', 16)
