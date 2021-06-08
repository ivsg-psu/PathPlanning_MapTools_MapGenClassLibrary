% script_demo_MapGenLibrary.m
% This is a script that shows the capabilities of the "MapGen" class of
% functions via demonstrations.

% Revision history:
% 2021_06_07:
% -- First write of the function, using the "Vis" library demo script as
% starter

% TO-DO:
% -- Reformat fcn_MapGen_polytopeShrinkToRadius

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

%% Show how inputs are checked
Twocolumn_of_numbers_test = [4 1; 3 0; 2 5];
fcn_MapGen_checkInputsToFunctions(Twocolumn_of_numbers_test, '2column_of_numbers');

%% Generate a set of polytopes from the Halton set
fig_num = 1;
Halton_range = [1 200]; % range of Halton points to use to generate the tiling
tiled_polytopes = fcn_MapGen_haltonVoronoiTiling(Halton_range,[1 1],fig_num);


%% Plot the polytopes
fig_num = 2;
line_width = 2;
axis_limits = [0 1 0 1];
fcn_MapGen_plotPolytopes(tiled_polytopes,fig_num,'r',line_width,axis_limits);

%% remove the edge polytopes that extend past the high and low points
fig_num = 3;
xlow = 0; xhigh = 1; ylow = 0; yhigh = 1;
bounding_box = [xlow ylow; xhigh yhigh];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(tiled_polytopes,bounding_box,fig_num);

%% shink the polytopes so that they are no longer tiled
des_radius = 0.03; % desired average maximum radius
sigma_radius = 0.002; % desired standard deviation in maximum radii
min_rad = 0.0001; % minimum possible maximum radius for any obstacle
shrink_seed = 1111; % seed used for randomizing the shrinking process
fig_num = 1;

[map_polytopes,all_pts] = ...
    fcn_MapGen_polytopeMapGen(...
    Halton_range,bounding_box,...
    des_radius,sigma_radius,min_rad,shrink_seed,fig_num);


%% Demonstrate visibility via ray casting

%% Demonstrate visibility via linear permutation

%% Demonstrate visibility via visibility matrix


%% Demonstrate radial occlusion function
fig_num = 3;
sensor_location = [0 0];
[is_visible_start,is_visible_end] = ...
    fcn_Vis_visViaRadialOcclusion(...
    wall_start,...
    wall_end,...
    sensor_location,...
    fig_num);

%% Demonstrate visibility on polytope maps

% Map generation

% Visibility matrix

% Radial occlusion

% Elimination of opposite-facing edges

% Compression of polytopes to line equivalents



%% 
if 1==0
    % generate map
    % generate Voronoi tiling from Halton points
    low_pt = 1; high_pt = 200; % range of Halton points to use to generate the tiling
    % remove the edge polytope that extend past the high and low points
    xlow = 0; xhigh = 1; ylow = 0; yhigh = 1;
    % shink the polytopes so that they are no longer tiled
    des_radius = 0.03; % desired average maximum radius
    sigma_radius = 0.002; % desired standard deviation in maximum radii
    min_rad = 0.0001; % minimum possible maximum radius for any obstacle
    shrink_seed = 1111; % seed used for randomizing the shrinking process
    [map_polytopes,all_pts] = ...
        fcn_MapGen_polytopeMapGen(...
        low_pt,high_pt,xlow,xhigh,ylow,yhigh,...
        des_radius,sigma_radius,min_rad,shrink_seed);  
    
    
    
    % starting (A) and finish (B) coordinates
    A.x = 0; A.y = 0.5; B.x = 1; B.y = 0.5;
    start = [A.x A.y size(all_pts,1)+1 -1 0];
    finish = [B.x B.y size(all_pts,1)+2 0 0];     

end

% Initiate the plot
fig = 103; % figure to plot on
line_spec = 'b-'; % edge line plotting
line_width = 2; % linewidth of the edge
axes_limits = [0 1 0 1]; % x and y axes limits
axis_style = 'square'; % plot axes style
fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);


% %% Test Case 1: Visibility from start
% cur_pt = start;
% end_pts = [all_pts; finish];
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(map_polytopes,cur_pt,end_pts);
% time1 = toc;
%
% % plot the map
% fig = 101; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
%
% %% Test Case 2: Visibility from finish
% cur_pt = finish;
% end_pts = [all_pts; start];
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(map_polytopes,cur_pt,end_pts);
% time2 = toc;
%
% % plot the map
% fig = 102; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
%
% %% Test Case 3: Intermediate point 1
% int_pt = 10;
% cur_pt = all_pts(int_pt,:);
% end_pts = [all_pts; start; finish];
% end_pts(int_pt,:) = []; % remove the intermediate point
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(map_polytopes,cur_pt,end_pts);
% time3 = toc;
%
% % plot the map
% fig = 103; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)

%% Test Case 3: Iterate through all points

% Initiate the plot
fig = 103; % figure to plot on
line_spec = 'b-'; % edge line plotting
line_width = 2; % linewidth of the edge
axes_limits = [0 1 0 1]; % x and y axes limits
axis_style = 'square'; % plot axes style
fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);

all_points_plus_start_and_end = [all_pts; start; finish];
N_points = length(all_pts(:,1));
for int_pt = 1:N_points
    cur_pt = all_points_plus_start_and_end(int_pt,:);
    end_pts = all_points_plus_start_and_end;  % start; finish];
    end_pts(int_pt,:) = []; % remove the intermediate point
    
    [clear_pts,~] = fcn_visibility_clear_and_blocked_points(map_polytopes,cur_pt,end_pts);
    
    % Prep lines for plotting
    N_clear = length(clear_pts(:,1));
    visibility_lines_start = ones(N_clear,1)*cur_pt;
    visibility_lines_end   = clear_pts(:,1:2);    
    visibility_lines_x = [visibility_lines_start(:,1) visibility_lines_end(:,1) NaN*visibility_lines_start(:,1)];
    visibility_lines_y = [visibility_lines_start(:,2) visibility_lines_end(:,2) NaN*visibility_lines_start(:,2)];
    visibility_lines_x = reshape(visibility_lines_x',N_clear*3,1);
    visibility_lines_y = reshape(visibility_lines_y',N_clear*3,1);
    
    
    plot(visibility_lines_x,visibility_lines_y,'g-','Linewidth',0.2)
    %     plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
    %     plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
    %     plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
    if 0==mod(int_pt,50)
        fprintf(1,'Now at: %.f of %.f\n',int_pt,N_points);
    end
end

% %% Test Case 4: Intermediate point 2
% int_pt = 50;
% cur_pt = all_pts(int_pt,:);
% end_pts = [all_pts; start; finish];
% end_pts(int_pt,:) = []; % remove the intermediate point
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(map_polytopes,cur_pt,end_pts);
% time4 = toc;
%
% % plot the map
% fig = 104; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
%
% %% Test Case 5: Intermediate point 3
% int_pt = 100;
% cur_pt = all_pts(int_pt,:);
% end_pts = [all_pts; start; finish];
% end_pts(int_pt,:) = []; % remove the intermediate point
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(map_polytopes,cur_pt,end_pts);
% time5 = toc;
%
% % plot the map
% fig = 105; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(map_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
%
%
% %% Custom map
% x1 = [0.4 0.6 0.5]; y1 = [0.5 0.5 0.7];
% x2 = [0.6 0.4 0.5]; y2 = [0.5 0.5 0.3];
% Xs = [x1 x2]; Ys = [y1 y2]; indices = [1 3; 4 6];
% [cust_polytopes,all_pts]=fcn_generate_custom_polytopes(Xs,Ys,indices);
%
%
% %% Test Case 6: adjacent polytopes
% cur_pt = start;
% end_pts = [all_pts; finish];
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(cust_polytopes,cur_pt,end_pts);
% time6 = toc;
%
% % plot the map
% fig = 106; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(cust_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
%
% %% Test Case 7: Merged Polytopes
% Xs = [0.4 0.5 0.6 0.5]; Ys = [0.5 0.3 0.5 0.7]; indices = [1 4];
% [cust_polytopes,all_pts]=fcn_generate_custom_polytopes(Xs,Ys,indices);
%
% cur_pt = start;
% end_pts = [all_pts; finish];
% tic
% [clear_pts,blocked_pts] = fcn_visibility_clear_and_blocked_points(cust_polytopes,cur_pt,end_pts);
% time7 = toc;
%
% % plot the map
% fig = 107; % figure to plot on
% line_spec = 'b-'; % edge line plotting
% line_width = 2; % linewidth of the edge
% axes_limits = [0 1 0 1]; % x and y axes limits
% axis_style = 'square'; % plot axes style
% fcn_plot_polytopes(cust_polytopes,fig,line_spec,line_width,axes_limits,axis_style);
% plot(clear_pts(:,1),clear_pts(:,2),'gx','linewidth',2)
% plot(blocked_pts(:,1),blocked_pts(:,2),'rx','linewidth',2)
% plot(cur_pt(1),cur_pt(2),'kx','linewidth',2)
%
