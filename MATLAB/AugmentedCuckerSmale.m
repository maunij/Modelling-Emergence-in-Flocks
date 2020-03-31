% Augmented Cucker-Smale Model

clear all
close all

% Set up variables
N=5;            % Number of birds
wx=5;           % Length of square to uniformly generate positions on
wv=5;           % Length of square to uniformly generate velocities on
numIt=2000;    % Number of iterations
K=1;            % Coupling constant
K2=1;
K3=1;
R=0.5;          % Separation Radius
beta=1;         % Rate of Decay
h=0.005;         % Time step

x = zeros(N,2);
v = zeros(N,2);

x(:,1) = wx*(-1/2 + rand(N,1));
x(:,2) = wx*(-1/2 + rand(N,1));
v(:,1) = wv*(-1/4 + rand(N,1));
v(:,2) = wv*(-1/2 + rand(N,1)); 



figure(3)
hold on
for i=1:N
    quiver(x(i,1),x(i,2),v(i,1),v(i,2));
end
title('Initial positions and velocities');
xlabel('x') ;
ylabel('y');
hold off

xcom = (1/N)*sum(x);
vcom = (1/N)*sum(v);

x = x - xcom;
v = v - vcom;

xnorm = zeros(numIt+1,1);
vnorm = zeros(numIt+1,1);

xnorm(1,1) = sqrt(sum(x(1:N,1).^2 + x(1:N,2).^2));
vnorm(1,1) = sqrt(sum(v(1:N,1).^2 + v(1:N,2).^2));

v_CS = zeros(N,2);
v_CSsum = zeros(N,2);
v_Relv = zeros(N,2);
v_Relx = zeros(N,2);
v_Relvsum = zeros(N,2);
v_Relxsum = zeros(N,2);

for n = 1:numIt
    for i = 1:N
        for j = 1:N
            if j == i
            else
                v_CS(j,:) = K*((1 + norm(x(i,:)-x(j,:))^2)^(-beta/2))*(v(j,:)-v(i,:));
                v_Relv(j,:) = (K2/(2*N))*((norm(x(i,:)-x(j,:))^(-2)))*dot(v(i,:)-v(j,:),x(i,:)-x(j,:))*(x(j,:)-x(i,:));
                v_Relx(j,:) = ((K3/N)*(2*norm(x(i,:)-x(j,:)))^(-1))*((norm(x(i,:)-x(j,:)))-2*R)*(x(j,:)-x(i,:));
            end
        end
        v_CSsum(i,:) = sum(v_CS);
        v_Relvsum(i,:) = sum(v_Relv);
        v_Relxsum(i,:) = sum(v_Relx);
    end
    x(1:N,:) = x(1:N,:) + h*v(1:N,:);
    v(1:N,:) = v(1:N,:) + h*(v_CSsum + v_Relvsum + v_Relxsum)
    
    xnorm(n+1,1) = sqrt(sum(x(1:N,1).^2 + x(1:N,2).^2));
    vnorm(n+1,1) = sqrt(sum(v(1:N,1).^2 + v(1:N,2).^2));
    
    n

end

ivector = zeros(numIt+1,1);
for i=1:numIt
    ivector(i+1)=i;
end

xnormav = sum(xnorm,2);
vnormav = sum(vnorm,2);

figure(1)
hold on
scatter(ivector, xnormav,2, 'filled');
title('Position Standard Deviation per Iteration');
xlabel('Number of Iterations') ;
ylabel('Standard Deviation');
hold off
figure(2)
hold on
scatter(ivector, vnormav,2, 'filled');
title('Velocity Standard Deviation per Iteration')
xlabel('Number of Iterations') 
ylabel('Standard Deviation')
hold off

x = x + xcom + vcom*numIt*h;
v = v + vcom;
figure(4)
hold on
for i=1:N
    quiver(x(i,1),x(i,2),v(i,1),v(i,2));
end
title('Final positions and velocities');
xlabel('x') ;
ylabel('y');
hold off