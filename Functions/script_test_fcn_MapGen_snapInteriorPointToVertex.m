% script_test_fcn_visibility_clear_and_blocked_points_global
% Tests: fcn_visibility_clear_and_blocked_points_global

%
% REVISION HISTORY:
%
% 2022_10_28 by S. Harnett
% -- first write of script
%%%%%%%%%%%%%%§

clear
clc
close all

%% add necessary directories
addpath([pwd '\..\Example_Map_Generation_Code'])
addpath([pwd '\..\PathPlanning_MapTools_MapGenClassLibrary\Functions'])
addpath([pwd '\..\PathPlanning_GeomTools_GeomClassLibrary\Functions'])

flag_do_plot = 1;
%% convex polytope
convex_polytope(1).vertices = [0 0; 1 1; -1 2; -2 1; -1 0; 0 0];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(convex_polytope);
pts_to_test = [0 0.5; -1 -1];
output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test)

% plot the map
if flag_do_plot
    fig = 111; % figure to plot on
    line_spec = 'b-'; % edge line plotting
    line_width = 2; % linewidth of the edge
    axes_limits = [-3 3 -3 3]; % x and y axes limits
    axis_style = 'square'; % plot axes style
    fcn_plot_polytopes(polytopes,fig,line_spec,line_width,axes_limits,axis_style);
    hold on
    box on
    xlabel('x [km]')
    ylabel('y [km]')
    plot(pts_to_test(1,1), pts_to_test(1,2),'rd')
    plot(pts_to_test(2,1), pts_to_test(2,2),'bd')
    plot(output_pts(1,1), output_pts(1,2),'rx')
    plot(output_pts(2,1), output_pts(2,2),'bx')
    legend('','polytope','pt. 1 init.','pt. 2 init.','pt. 1 snapped','pt. 2 snapped')
end
