% CuckerSmaleEverything.m

% This is a modification of the CS algorithm to add collision avoidance,
% interparticle bonding forces and a 7 nearest neighbour interaction.

clear all
close all

% Set up variables
N=20;           % Number of birds
M=1;            % Number of tests
wx=10;          % Length of square to uniformly generate positions on
wv=10;          % Length of square to uniformly generate velocities on
numIt=1000;     % Number of iterations
K=1;            % Coupling constant
K2=1;
K3=1;
R=1;          % Separation Radius
beta=1;         % Rate of Decay
h=0.01;         % Time step

% Uniformly generate M positions and M velocites in R^2 on the square of
% length w. Convert to CoM coordinates.
x = zeros(N,2*M);
v = zeros(N,2*M);
for n = 1:M
    x(:,2*n-1) = wx*(-1/2 + rand(N,1));
    x(:,2*n) = wx*(-1/2 + rand(N,1));
    v(:,2*n-1) = wv*(-1/2 + rand(N,1));
    v(:,2*n) = wv*(-1/2 + rand(N,1)); 
end

% For interest, quiver plot the first test

figure(3)
hold on
for i=1:N
    quiver(x(i,1),x(i,2),v(i,1),v(i,2));
end
hold off
% Convert to Centre of Mass Coordinates.

xcom = (1/N)*sum(x);
vcom = (1/N)*sum(v);

x = x - xcom;
v = v - vcom;

% Create the matrix of standard deviations

xnorm = zeros(numIt+1,M);
vnorm = zeros(numIt+1,M);

for m = 1:M
    xnorm(1,m) = sqrt(sum(x(1:N,2*m-1).^2 + x(1:N,2*m).^2));
    vnorm(1,m) = sqrt(sum(v(1:N,2*m-1).^2 + v(1:N,2*m).^2));
end

% Now perform the Cucker Smale test on each pair of columns.

for n = 1:numIt
    for m = 1:M
        a = zeros(N,N);
        for i=1:N
            % Calculate 7 nearest neighbours to agent i (8 is used since
            % this also counts i itself as a nearest neighbour, in which
            % case the entry a(i,i) is just 1)
            u = knnsearch(x(:,(2*m-1):2*m),x(i,(2*m-1):2*m),'K',8,'Distance','euclidean');
            for j =1:N
                % Use ismember as an indicator function
                a(i,j) = ismember(x(j,(2*m-1):2*m),x(u,(2*m-1):2*m),'rows')*K*(1 + norm(x(i,(2*m-1):2*m)-x(j,(2*m-1):2*m))^2)^(-beta/2);
            end
        end
        a;
        L = diag(sum(a,2)) - a;
        
        % Add collision avoidance and inter-particle bonding forces
        
        va = zeros(N,2);
        vb = zeros(N,2);
        vc = zeros(N,2);
        vd = zeros(N,2);
        
        for i=1:N
            for j =1:N
                if j == i
                else
                    va(j,:) = (K2/(2*(norm(x(i,(2*m-1):2*m)-x(j,(2*m-1):2*m)))^2))*dot(v(i,(2*m-1):2*m)-v(j,(2*m-1):2*m),x(i,(2*m-1):2*m)-x(j,(2*m-1):2*m))*(x(j,(2*m-1):2*m)-x(i,(2*m-1):2*m));
                    vb(j,:) = (K3/(2*norm(x(i,(2*m-1):2*m)-x(j,(2*m-1):2*m))))*((norm(x(i,(2*m-1):2*m)-x(j,(2*m-1):2*m)))-2*R)*(x(j,(2*m-1):2*m)-x(i,(2*m-1):2*m));
                end
            end
            vc(i,:) = sum(va);
            vd(i,:) = sum(vb);
        end
        % Do CS
        x(1:N,2*m-1) = x(1:N,2*m-1) + h*v(1:N,2*m-1);
        x(1:N,2*m) = x(1:N,2*m) + h*v(1:N,2*m);
        v(1:N,2*m-1) = (eye(N) - h*L)*v(1:N,2*m-1)+h*(vc(:,1)+vd(:,1));
        v(1:N,2*m) = (eye(N) - h*L)*v(1:N,2*m)+h*(vc(:,2)+vd(:,2));  
        
        % Calculate norms
        xnorm(n+1,m) = sqrt(sum(x(1:N,2*m-1).^2 + x(1:N,2*m).^2));
        vnorm(n+1,m) = sqrt(sum(v(1:N,2*m-1).^2 + v(1:N,2*m).^2));
    end
end

% Take average over M tests per iteration

ivector = zeros(numIt+1,1);
for i=1:numIt
    ivector(i+1)=i;
end
xnormav = (1/M)*sum(xnorm,2);
vnormav = (1/M)*sum(vnorm,2);

% Scatter plot the averages per iteration. Recommend commenting the x_m,
% x_M plots as they are very generous!

figure(1)
hold on
scatter(ivector, xnormav,2, 'filled');
%yline(xm_min, '--','Minimum x_m');
%yline(xm_max, '--','Maximum x_m');
%yline(xm_av, '--','Average x_m');
%yline(xM_min, '--','Minimum x_M');
%yline(xM_max, '--','Maximum x_M');
%yline(xM_av, '--','Average x_M');
title('Position Standard Deviation per Iteration');
xlabel('Number of Iterations') ;
ylabel('Standard Deviation');
hold off
figure(2)
hold on
scatter(ivector, vnormav,2, 'filled');
%scatter(ivector, vM_min,1/2);
%scatter(ivector, vM_max,1/2);
%scatter(ivector, vM_av,1/2, 'filled');
title('Velocity Standard Deviation per Iteration')
xlabel('Number of Iterations') 
ylabel('Standard Deviation')
hold off

% Finally, we look at a quiver plot of the final velocities
x = x + xcom + vcom*numIt*h;
v = v + vcom;
figure(4)
hold on
for i=1:N
    quiver(x(i,1),x(i,2),v(i,1),v(i,2));
end
hold off

% Note that the standard deviation is scaled by a factor of sqrt(N).
