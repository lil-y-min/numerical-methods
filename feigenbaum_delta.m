% Compute the Feigenbaum delta

num_doublings=11; delta=zeros(1,num_doublings); delta(1)=5;

m = zeros(1, num_doublings+1);   
m(1) = 2;                        % m_0
m(2) = 1 + sqrt(5);              % m_1

% Iteration controls
tol_g  = 1e-14;                  
tol_mu = 1e-14;                  
max_newton = 50;                 

for n = 2:num_doublings          
    N = 2^n;                     

    % initial guess for m_n using last two m's and current delta
    mu = m(n) + (m(n) - m(n-1)) / delta(n-1);   
    
    % keep guess in the desired interval (m_{n-1}, 4]
    if ~(mu > m(n)), mu = 0.5*(m(n) + 4); end
    if (mu > 4),     mu = 0.5*(m(n) + 4); end

    % Newton iterations for this n
    for it = 1:max_newton
        
        % iterate logistic map and its mu-derivative N times from x0=1/2
        x  = 0.5;                % x_0
        xp = 0.0;                % x'_0 
        
        for k = 1:N
            x_old = x;
            x  = mu * x_old * (1 - x_old);
            xp = x_old*(1 - x_old) + mu * xp * (1 - 2*x_old);  
        end

        g  = x - 0.5;            % g(mu) = x_N - 1/2
        gp = xp;                 % g'(mu) = x'_N

        if abs(g) < tol_g, break; end
        if abs(gp) < eps, gp = sign(gp) * eps; end  % avoid zero divisor

        dmu    = g / gp;
        mu_new = mu - dmu;

        if abs(mu_new - mu) < tol_mu, mu = mu_new; break; end
        mu = mu_new;
    end

    m(n+1)    = mu;
    delta(n)  = (m(n) - m(n-1)) / (m(n+1) - m(n));
end

% Output results
fprintf('n        delta(n)\n');
for n=1:num_doublings
    fprintf('%2g %18.15f\n',n,delta(n));
end
