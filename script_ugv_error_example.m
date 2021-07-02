% Scratch Paper Examples for Polytope Error Bubble Generation

% Revision history:
% 2021_03_26 - Nick Carder
% -- First write of the function, 
% 2021_06_28 - S. Brennan
% -- Reworked to be compatible with MapGen library

close all

%% Set up workspace
clear flag_was_run_before  % Force init to always run?

if ~exist('flag_was_run_before','var')
    
    clc
    close all
    
    % add necessary directories
    addpath([pwd '\Functions'])
    %     addpath([pwd '\GeomClassLibrary\Functions'])
    %     addpath([pwd '\MapGenClassLibrary\Functions'])
    %     addpath([pwd '\Plotting'])
    %     addpath([pwd '\Map_Generation\polytope_generation'])
    %     addpath([pwd '\Map_Generation\polytope_editing'])
    %     addpath([pwd '\Map_Generation\polytope_calculation'])
    
    flag_was_run_before = 1;
end



%% Generating Error Bubbles and Plotting

%determine the size of the map with stretch
stretch = [200, 200]; % stretch in the x and y directions
%stretch = [1, 1]; % stretch in the x and y directions

% %create polytopes
% polytopes = fcn_polytope_generation_halton_voronoi_tiling(1,1000,stretch);

% Generate a set of polytopes from the Halton set
fig_num = 12;
Halton_range = [1 1000]; % range of Halton points to use to generate the tiling
polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,stretch,fig_num);
title('Halton set');


%plot polytopes to verify map size, polytopes will still be voronoi tiles
fcn_plot_polytopes(polytopes,1,'b',2,[0 200 -100 100],'square');

%trim polytopes to fit desired map size
trim_polytopes = fcn_polytope_editing_remove_edge_polytopes(polytopes,0,200,0,200);

%verify
fcn_plot_polytopes(trim_polytopes,2,'b',2);

%shrink polytopes to create space
scale = ones(1,size(trim_polytopes,2));
des_rad = 1; sigma_radius = 0.5; min_rad = 0.25;
shrunk_polytopes2=fcn_polytope_editing_shrink_to_average_max_radius_with_variance(trim_polytopes,des_rad,sigma_radius,min_rad);

%verify
fcn_plot_polytopes(shrunk_polytopes2,3,'b',2);







%generate error bubbles
[err] = err_ugv_bubble_v3(shrunk_polytopes2, 0, 5);

%verify
figure('name','UGV Positioning Bubbles')
ax.ugv1=gca;
hold on
for ii=1:length(shrunk_polytopes2)
    plot(shrunk_polytopes2(ii).vertices(:,1),shrunk_polytopes2(ii).vertices(:,2),'k')
    plot(ax.ugv1,err(ii).circ_x(err(ii).bubble),err(ii).circ_y(err(ii).bubble),'r')
end
hold off
xlabel('X Distance [m]')
ylabel('Y Distance [m]')
title('UGV Positioning Bubbles')
legend('Real Object','Perceived Object')
grid on
axis equal

%% Testing work with contour plots

%generating R and beta values for a grid of points
x=linspace(0,200,200);
y=linspace(100,-100,200);
[X,Y] = meshgrid(x,y);

R=sqrt(X.^2+Y.^2);  % perceived distance, nearly constant for flat object
beta=rad2deg(atan(Y./X));
kappa=zeros(size(R));

[err.x, err.y] = err_ugv_v3 ({R, beta, kappa}, ...
    {0.08, 0.08, 0.08}, {0.03, 0.03, 0.4}, {0.02, -0.05});

err.R = sqrt(err.x.^2 + err.y.^2);

figure('name','Error in X')
contour(x,y,err.x,'ShowText','on')
h=colorbar;
title('Error in X Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Error in Y')
contour(x,y,err.y,'ShowText','on')
h=colorbar;
title('Error in Y Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Total Error (Bubble Radius) [m]')
contour(x,y,err.R,'ShowText','on')
h=colorbar;
title('Total Error (Bubble Radius) [m]')
xlabel('X [m]')
ylabel('Y [m]')