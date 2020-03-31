% Animation.m
%
% This is designed to work with the .mat file created using Boids.m in
% order to create a gif animation of results. This was also used to create
% the animations for the Cucker-Smale model.

load( 'Boids Simulation.mat' );

numIt = 100;
N = 30;
h = figure;
axis off
set(gca,'XTick',[], 'YTick', [])
set(gca, 'Color', 'White')
set(gcf, 'Color', 'White')
filename =('BoidsAnimated.gif');
for n = 1:numIt
    for i = 1:N
        axis off
        set(gca,'XTick',[], 'YTick', [])
        q = quiver(Res(i,1,n),Res(i,2,n),Res(i,3,n),Res(i,4,n), 'LineWidth',2, 'MaxHeadSize', 1, 'Marker','o', 'MarkerSize', 20);
        hold on
    end
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
        imwrite(imind,cm,filename,'gif', 'DelayTime',0.001, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.001); 
       end 
    hold off
end