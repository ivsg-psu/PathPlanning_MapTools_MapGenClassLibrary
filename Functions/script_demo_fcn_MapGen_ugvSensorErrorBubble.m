% script_test_fcn_MapGen_ugvSensorErrorBubble
% Tests: fcn_MapGen_ugvSensorErrorBubble

% 
% REVISION HISTORY:
% 
% 2021_03_26 - Nick Carder
% -- First write of the function, 
% 2021_06_28 - S. Brennan
% -- Reworked to be compatible with MapGen library
% 2021_07_07 by Sean Brennan
% -- additional edits on script
% 2025_07_11 by S. Brennan
% -- still have no idea what this script does!

close all


%% Generating Error Bubbles and Plotting

% create polytopes
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 1000];
AABBs = [0 0 1 1];
mapStretchs = [200 200];
[polytopes] = fcn_MapGen_voronoiTiling(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));

% Plot the polytopes
fig_num = 22;
% line_width = 2;
% axis_limits = [0 200 0 200];
% axis_stype = 'square';
% fcn_MapGen_OLD_plotPolytopes(polytopes,fig_num,'b',line_width,axis_limits,axis_stype);
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [0 0 1];
fillFormat = [];
h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>


% remove the edge polytopes that extend past the high and low points
fig_num = 23;
bounding_box = [0 0; 200 200];
trimmed_polytopes = ...
    fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);

%shrink polytopes to create space
fig_num = 24;
des_rad = 1; sigma_radius = 0.5; min_rad = 0.25;
shrunk_polytopes2=fcn_MapGen_polytopesShrinkToRadius(...
    trimmed_polytopes,des_rad,sigma_radius,min_rad,fig_num);




% generate error bubbles via fcn_MapGen_ugvSensorErrorBubble
[err] = fcn_MapGen_ugvSensorErrorBubble(shrunk_polytopes2, 0, 5);
assert(true);

% Convert error bounds into polytope structure
error_polytopes = shrunk_polytopes2; % Initialize the structure
for ii=1:length(shrunk_polytopes2)
    error_polytopes(ii).vertices = [err(ii).circ_x(err(ii).bubble)', err(ii).circ_y(err(ii).bubble)'];
end
error_polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(error_polytopes);

%verify
h_fig = figure('name','UGV Positioning Bubbles');
fig_num = h_fig.Number;
ax.ugv1=gca;
% hold on
% for ii=1:length(shrunk_polytopes2)
%     plot(shrunk_polytopes2(ii).vertices(:,1),shrunk_polytopes2(ii).vertices(:,2),'k')
%     plot(ax.ugv1,err(ii).circ_x(err(ii).bubble),err(ii).circ_y(err(ii).bubble),'r')
% end
% hold off

% fcn_MapGen_OLD_plotPolytopes(shrunk_polytopes2,fig_num,'k',line_width,axis_limits,axis_stype);
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [0 0 0];
fillFormat = [];
h_plot = fcn_MapGen_plotPolytopes(shrunk_polytopes2, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

% fcn_MapGen_OLD_plotPolytopes(error_polytopes,fig_num,'r',line_width,axis_limits,axis_stype);
plotFormat.LineWidth = 2;
plotFormat.MarkerSize = 10;
plotFormat.LineStyle = '-';
plotFormat.Color = [1 0 0];
fillFormat = [];
h_plot = fcn_MapGen_plotPolytopes(error_polytopes, (plotFormat), (fillFormat), (fig_num)); 

xlabel('X Distance [m]')
ylabel('Y Distance [m]')
title('UGV Positioning Bubbles')
legend('Real Object','Perceived Object')
grid on
axis equal

%%%%
% Testing work with contour plots, using function fcn_MapGen_ugvSensorError

%generating R and beta values for a grid of points
x = linspace(0,200,200);
y = linspace(100,-100,200);
[X,Y] = meshgrid(x,y);

R = sqrt(X.^2+Y.^2);  % perceived distance, nearly constant for flat object
beta = rad2deg(atan(Y./X));
kappa = zeros(size(R));

clear err;
[err.x, err.y, err.z] = fcn_MapGen_ugvSensorError({R, beta, kappa}, ...
    {0.08, 0.08, 0.08}, {0.03, 0.03, 0.4}, {0.02, -0.05});

err.R = sqrt(err.x.^2 + err.y.^2);

figure('name','Error in X')
contour(x,y,err.x,'ShowText','on')
% h=colorbar;
title('Error in X Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Error in Y')
contour(x,y,err.y,'ShowText','on')
% h=colorbar;
title('Error in Y Dimension [m]')
xlabel('X [m]')
ylabel('Y [m]')

figure('name','Total Error (Bubble Radius) [m]')
contour(x,y,err.R,'ShowText','on')
h=colorbar;
title('Total Error (Bubble Radius) [m]')
xlabel('X [m]')
ylabel('Y [m]')


