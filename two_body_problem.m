e=0.7; m1=1; m2=4;
T=2*pi./(1-e).^1.5; tspan=linspace(0,T,1000);
options=odeset('RelTol',1.e-6);
%%%%% Solve differential equations for x and y 

%%%%% initial conditions in COM coords (k=1, r_min=1)
x0 = -1; y0 = 0; vx0 = 0; vy0 = sqrt(1+e);
z0 = [x0; y0; vx0; vy0];

[t,z] = ode45(@two_body, tspan, z0, options);   % z = [x y vx vy]
x = z(:,1);  y = z(:,2);

%%%%% Determine x1, y1 and x2, y2
M = m1 + m2;
x1 = (m2/M)*x;   y1 = (m2/M)*y;      % position of m1
x2 = -(m1/M)*x;  y2 = -(m1/M)*y;     % position of m2

%%%%% Graphics %%%%%%%%%%%%%%
k=0.1;
R1=k*(m1)^(1/3); R2=k*(m2)^(1/3); %radius of masses
theta = linspace(0,2*pi); 
figure; axis equal; hold on; set(gcf,'color','w');
axis off; 
xlim([-2,5]); ylim([-2.5,2.5]);
planet=fill(R1*cos(theta)+x1(1), R1*sin(theta)+y1(1),'b'); 
sun=fill(R2*cos(theta)+x2(1), R2*sin(theta)+y2(1),'r'); 
pause(1);
nperiods=5; %number of periods to plot
for j=1:nperiods
    for i=1:length(t)
        planet.XData=R1*cos(theta)+x1(i); planet.YData=R1*sin(theta)+y1(i); 
        sun.XData=R2*cos(theta)+x2(i); sun.YData=R2*sin(theta)+y2(i); 
        drawnow;
    end
end

%%%%% Local function for differential equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dzdt = two_body(~, z)
    x = z(1); y = z(2); vx = z(3); vy = z(4);
    r2 = x*x + y*y;
    r3 = r2^(3/2);
    ax = -x / r3;                 % k = 1
    ay = -y / r3;
    dzdt = [vx; vy; ax; ay];
end
