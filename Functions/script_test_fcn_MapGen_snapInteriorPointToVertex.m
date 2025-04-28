% script_test_fcn_MapGen_snapInteriorPointToVertex
% tests function: fcn_MapGen_snapInteriorPointToVertex

%
% REVISION HISTORY:
%
% 2024_04_19 by S. Harnett
% -- first write of script
% 2025_04_28 by S. Harnett
% -- fix legends
%%%%%%%%%%%%%%§


%% run snapInteriorPointToVertex function
flag_do_plot = 1;
% convex polytope
convex_polytope(1).vertices = [0 0; 1 1; -1 2; -2 1; -1 0; 0 0];
convex_polytope(2).vertices = [convex_polytope(1).vertices(:,1) + 4, convex_polytope(1).vertices(:,2) - 2];
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(convex_polytope);
pts_to_test = [0 0.5; -1 -1; 4 -1; 4.1 -1];
output_pts = fcn_MapGen_snapInteriorPointToVertex(polytopes, pts_to_test);

% plot the map
if flag_do_plot
    fig = 111; % figure to plot on
    line_spec = 'b-'; % edge line plotting
    line_width = 2; % linewidth of the edge
    axes_limits = [-3 5 -3 5]; % x and y axes limits
    axis_style = 'square'; % plot axes style
    figure
    fcn_MapGen_plotPolytopes(polytopes,fig,line_spec,line_width,axes_limits,axis_style);
    hold on
    box on
    title('function result')
    xlabel('x [km]')
    ylabel('y [km]')
    plot(pts_to_test(1,1), pts_to_test(1,2),'rd')
    plot(pts_to_test(2,1), pts_to_test(2,2),'bd')
    plot(pts_to_test(3,1), pts_to_test(3,2),'gd')
    plot(pts_to_test(4,1), pts_to_test(4,2),'kd')
    plot(output_pts(1,1), output_pts(1,2),'rx')
    plot(output_pts(2,1), output_pts(2,2),'bx')
    plot(output_pts(3,1), output_pts(3,2),'gx')
    plot(output_pts(4,1), output_pts(4,2),'kx')
    legend('polytope','pt. 1 init.','pt. 2 init.','pt. 3 init.','pt. 4 init.','pt. 1 snapped','pt. 2 snapped','pt. 3 snapped','pt. 4 snapped')
end

assert(true)
