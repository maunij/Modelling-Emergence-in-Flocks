% Boids.m
% 
% A topological Boids model, in which the alignment vector is calculated
% using a k nearest neighbour search. This outputs a .mat file that can be
% used in AnimationBoids.m to create a GIF of the animation.

clear all
close all

N = 30;           % Number of boids
w = 80;          % Length of square in which to generate boids
vlim = 10;       % Maximum velocity
xMax = 75;       
xMin = -75;      
yMax = 75;       
yMin = -75;      % Boundaries
numIt = 1000;     % Number of Iterations.

% Initial Conditions, uniformly generated on the square.

x = w*(-1/2 + rand(N,2));
v = (-1/2 + rand(N,2));

% Vector of average velocites (for calculating alignment parameter)
v_unit = zeros(N,2);
phi = zeros(numIt, 1);

% figure
% for i = 1:N
%     quiver(x(i,1),x(i,2),v(i,1),v(i,2),'LineWidth',2, 'Color', 'k')
%     hold on
% end
% hold off

v1 = zeros(N,2); % Alignment Vector
v2 = zeros(N,2); % Cohesion Vector
v3 = zeros(N,2); % Separation Vector
v4 = zeros(N,2); % Boundary Avoidance Vector

for n = 1:numIt
    % Calculate Centre of Mass of Flock (sum columns and divide by N)
    C = sum(x)/N;
    % Move boids to new positions according to 3 rules
    
    % Alignment, align with close neighbours. Iterate through each boid and
    % find 7 neighbours. If j is a nearest neighbour we add its velocity to v1.
    
    for i = 1:N
        neighbours = 0;
        a = zeros(1,2);
        u = knnsearch(x(:,:),x(i,:),'K',8,'Distance','euclidean');
        for j = 1:N
            if ismember(x(j,:),x(u,:),'rows') == 1
                a = a + v(j,:);
            end
        end
        v1(i,:) = a/norm(a);
        
     % This is an alternate, metric based alignment model that takes into
     % account all boids within a set radius.
     
%         for j = setdiff(1:N,i)
%             if norm(x(j,:) - x(i,:)) < 10;
%                 neighbours = neighbours + 1;
%                 a = a + v(j,:);
%             end
%         end
%         if neighbours == 0;
%             v1(i,:) = 0;
%         else
%             v1(i,:) = (a/neighbours)/norm(a/neighbours);
%         end
    end   
    
% Global Cohesion, move unit distance towards the centre of the flock
%     for i = 1:N
%         v2(i,:) = (C - x(i,:))/norm(C - x(i,:));
%     end

% Topological Cohesion, move unit distance towards the centre of mass of
% nearest 7 neighbours.
    
    for i = 1:N
        neighbours = 0;
        b = zeros(N,2);
        u2 = knnsearch(x(:,:),x(i,:),'K',8,'Distance','euclidean');
        for j = 1:N
            if ismember(x(j,:),x(u2,:),'rows') == 1
                b(j,:) = x(j,:);
            end
        end
        B = sum(b)/N;
        v2(i,:) = (B - x(i,:))/norm(B - x(i,:));
    end
 
    % Separation, avoid if within a set distance. Similar to alignment but
    % distance gets taken away instead.
    for i = 1:N
        b = zeros(1,2);
        for j = setdiff(1:N,i)
            if norm(x(j,:) - x(i,:)) < 6
                b = b - (x(j,:)-x(i,:));
            end
        end
        if b == 0
            v3(i,:) = b;
        else
        v3(i,:) = b;
        end
    end   

% Finally, we'll add some sort of boundary condition so we can view the
% animation without boids flying off the screen. If boids fly outside of
% the square they will turn back.

% Boundary Conditions create a milling phenomenon.
% 
%     for i = 1:N
%         c = zeros(1,2);
%         if x(i,1) > xMax
%             c = vlim*[-sqrt(2),sqrt(2)];
%         else if x(i,1) < xMin
%             c = vlim*[sqrt(2),-sqrt(2)];
%             else if x(i,2) > yMax
%                 c = vlim*[-sqrt(2),-sqrt(2)];
%                     else if x(i,2) < yMin
%                         c = vlim*[sqrt(2),sqrt(2)];
%                     end
%                 end
%             end
%         end
%         if any(c) == 1
%             v4(i,:) = c/norm(c);
%         end
%     end
% 
% Standard Boundary conditions
%
%     for i = 1:N
%         c = zeros(1,2);
%         if x(i,1) > xMax
%             c(1,1) = c(1,1) - vlim;
%             else if x(i,1) < xMin
%               c(1,1) = c(1,1) + vlim;
%               else if x(i,2) > yMax
%                   c(1,2) = c(1,2) - vlim;
%                       else if x(i,2) < yMin
%                           c(1,2) = c(1,2) + vlim;
%                     end
%                 end
%             end
%         end
%         if any(c) == 1
%             v4(i,:) = c/norm(c);
%         end
%     end
%     
    
% Now update accordingly
     for i = 1:N
          v(i,:) = v(i,:) + v1(i,:) + v2(i,:) + v3(i,:);
%Limit velocities
          if norm(v(i,:)) > vlim
              v(i,:) = vlim*v(i,:)/norm(v(i,:));
          end
          x(i,:) = x(i,:) + v(i,:);
          
          v_unit(i,:) = v(i,:)/norm(v(i,:));
          
          % Save data in a 3 dimensional array
          Res(i,:,n) = [x(i,:),v(i,:)];
          
     end
     hold off
     phi(n,1) = norm(mean(v_unit));
     phi(n,1);
%     hold off
%     pause(0.001);

end

% Plot the results

X = zeros(numIt,1);
for m = 1:numIt
    X(m,1) = m;
end

plot(X,phi)
mean(phi)
std(phi)
save('Boids Simulation.mat', 'Res' );
