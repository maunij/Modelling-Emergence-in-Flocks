% Cucker Smale Algorithm in R^2
%
% Script for running a variable number of tests on the Cucker Smale model
% and calculating the decay of average position/velocity standard
% deviations per iteration. Variables are currently set to produce the
% nicest scatter plot but can be adjusted as necessary. Keep K*h small
% though as this is only a first order approximation (Euler method).

clear all
close all

% Set up variables

N=20;           % Number of birds
M=1;            % Number of tests
wx=10;          % Length of square to uniformly generate positions on
wv=10;          % Length of square to uniformly generate velocities on
numIt=1000;     % Number of iterations
K=1;            % Coupling constant
beta=1;         % Decay Constant, Ha-Liu Tests will be wrong if not 1
h=0.005;        % Time step

% Uniformly generate M positions and M velocites in R^2 on the square of
% length w. Positions have mean (0,0), and velocities have mean (wv/4,0).
x = zeros(N,2*M);
v = zeros(N,2*M);
for n = 1:M
    x(:,2*n-1) = wx*(-1/2 + rand(N,1));
    x(:,2*n) = wx*(-1/2 + rand(N,1));
    v(:,2*n-1) = wv*(-1/4 + rand(N,1));
    v(:,2*n) = wv*(-1/2 + rand(N,1)); 
end

% Quiver plot the initial conditions on the first test.

figure(3)
hold on
for i=1:N
    quiver(x(i,1),x(i,2),v(i,1),v(i,2),'LineWidth',2, 'MaxHeadSize', 1);
end
title('Initial positions and velocities');
xlabel('x') ;
ylabel('y');
hold off

% Convert to Centre of Mass Coordinates.

xcom = (1/N)*sum(x);
vcom = (1/N)*sum(v);

x = x - xcom;
v = v - vcom;

% Create the matrix of standard deviations.

xnorm = zeros(numIt+1,M);
vnorm = zeros(numIt+1,M);

for m = 1:M
    xnorm(1,m) = sqrt(sum(x(1:N,2*m-1).^2 + x(1:N,2*m).^2));
    vnorm(1,m) = sqrt(sum(v(1:N,2*m-1).^2 + v(1:N,2*m).^2));
end

% Get Ha-Liu Bounds. THIS IS ONLY TRUE FOR beta = 1! Keep wx, wv
% small if you want close bounds, otherwise these will look quite extreme.

x0min = min(xnorm(1,:));
x0max = max(xnorm(1,:));
x0av = (1/M)*sum(xnorm(1,1:M));
v0min = min(vnorm(1,:));
v0max = max(vnorm(1,:));
v0av = (1/M)*sum(vnorm(1,1:M));

am_min = (sqrt(2)*x0min + sqrt(2*(x0min)^2 + 1))*exp(-sqrt(2)*v0min*1/(N*K));
bm_min = max(1,am_min);
am_max = (sqrt(2)*x0max + sqrt(2*(x0max)^2 + 1))*exp(-sqrt(2)*v0max*1/(N*K));
bm_max = max(1,am_max);
am_av = (sqrt(2)*x0av + sqrt(2*(x0av)^2 + 1))*exp(-sqrt(2)*v0av*1/(N*K));
bm_av = max(1,am_av);

aM_min = (sqrt(2)*x0min + sqrt(2*(x0min)^2 + 1))*exp(sqrt(2)*v0min*1/(N*K));
aM_max = (sqrt(2)*x0max + sqrt(2*(x0max)^2 + 1))*exp(sqrt(2)*v0max*1/(N*K));
aM_av = (sqrt(2)*x0av + sqrt(2*(x0av)^2 + 1))*exp(sqrt(2)*v0av*1/(N*K));

xm_min = ((bm_min)^2 - 1)/(2*sqrt(2)*bm_min);
xm_max = ((bm_max)^2 - 1)/(2*sqrt(2)*bm_max);
xm_av = ((bm_av)^2 - 1)/(2*sqrt(2)*bm_av);
xM_min = ((aM_min)^2 - 1)/(2*sqrt(2)*aM_min);
xM_max = ((aM_max)^2 - 1)/(2*sqrt(2)*aM_max);
xM_av = ((aM_av)^2 - 1)/(2*sqrt(2)*aM_av);

vM_min = zeros(numIt+1,1);
vM_max = zeros(numIt+1,1);
vM_av = zeros(numIt+1,1);

for n = 0:numIt
    vM_min(n+1,1) = v0min*exp(-(K/(1+2*(xM_min)^2)^(1/2))*n*h);
    vM_max(n+1,1) = v0max*exp(-(K/(1+2*(xM_max)^2)^(1/2))*n*h);
    vM_av(n+1,1) = v0av*exp(-(K/(1+2*(xM_av)^2)^(1/2))*n*h);
end

% Now perform the Cucker Smale test on each pair of columns. To speed up
% calculations I have used the Laplacian matrix method introduced in Cucker
% & Smale's original paper (2007).

for n = 1:numIt
    for m = 1:M
        a = zeros(N,N);
        for i=1:N
            for j =i+1:N
                a(i,j) = K*(1 + norm(x(i,(2*m-1):2*m)-x(j,(2*m-1):2*m))^2)^(-beta/2);
            end
        end
        A = a + a';
        L = diag(sum(A,2)) - A;
        
        x(1:N,2*m-1) = x(1:N,2*m-1) + h*v(1:N,2*m-1);
        x(1:N,2*m) = x(1:N,2*m) + h*v(1:N,2*m);
        v(1:N,2*m-1) = (eye(N) - h*L)*v(1:N,2*m-1);
        v(1:N,2*m) = (eye(N) - h*L)*v(1:N,2*m);  
        
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
%yline(xm_min, '--','Minimum m');
%yline(xm_max, '--','Maximum m');
yline(xm_av, '-','Minimum Bound m','color',[0.8500 0.3250 0.0980]); %change back to avg
%yline(xM_min, '--','Minimum M');
%yline(xM_max, '--','Maximum M');
yline(xM_av, '-','Maximum Bound M','color',[0.8500 0.3250 0.0980]);
title('Position Standard Deviation per Iteration');
xlabel('Number of Iterations') ;
ylabel('Standard Deviation');
hold off
figure(2)
hold on
scatter(ivector, vnormav,2, 'filled');
%scatter(ivector, vM_min,1/2);
%scatter(ivector, vM_max,1/2);
scatter(ivector, vM_av,1/2, 'filled');
title('Velocity Standard Deviation per Iteration')
xlabel('Number of Iterations') 
ylabel('Standard Deviation')
hold off

% Finally, we look at a quiver plot of the final velocities
x = x + xcom + vcom*h*numIt;
v = v + vcom;
figure(4)
hold on
for i=1:N
    quiver(x(i,1),x(i,2),v(i,1),v(i,2), 'LineWidth',2, 'MaxHeadSize', 1);
end
title('Final positions and velocities');
xlabel('x') ;
ylabel('y');
hold off

% Note that the standard deviation is scaled by a factor of sqrt(N).
