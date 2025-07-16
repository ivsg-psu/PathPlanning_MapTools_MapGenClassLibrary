function polyMapStats = fcn_MapGen_polytopesStatistics(polytopes, varargin)
% fcn_MapGen_polytopesStatistics calculates key statistics for the
% polytopes
%
% FORMAT:
%
% [polyMapStats] = fcn_MapGen_polytopesStatistics( polytopes, (fig_num))
%
% INPUTS:
%
%     POLYTOPES: polytopes to analyze
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     polyMapStats: a structure whose fields contain the polytope map's
%     statistics
%     each of the fields in the struct are as follows:
%     AABB is the axis aligned bounding box for the polytope field, bounding all polytopes
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_polytopeFindVertexAngles
%     fcn_MapGen_plotPolytopes
%     fcn_Path_findSensorHitOnWall
%
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopesStatistics
%
% This function was written on 2021_07_11 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2021_07_11 - S. Brennan
% -- first write of function
% 2022_06_17 - S. Harnett
% -- numerous features added to support length cost ratio prediction
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_10 by Sean Brennan
% -- changed fcn_MapGen_findIntersectionOfSegments to use
% fcn_Path_findSensorHitOnWall instead, as the Path function is much more
% tested/debugged and regularly updated
% -- updated some variable naming
% -- clean-up of comments
% 2025_07_10 by Sean Brennan
% -- updated header debugging and input area to fix global flags

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS");
    MATLABFLAG_MAPGEN_FLAG_DO_DEBUG = getenv("MATLABFLAG_MAPGEN_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_MAPGEN_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_MAPGEN_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end


%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (0==flag_max_speed)
    if flag_check_inputs
        % Are there the right number of inputs?
        if nargin < 1 || nargin > 2
            error('Incorrect number of input arguments')
        end

        % Check the polytopes input
        fcn_DebugTools_checkInputsToFunctions(...
            polytopes, 'polytopes');

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  2 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; 
        flag_do_plot = 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Npolys = length(polytopes);

% Loop through the polytopes, count the radii
Nradii = 0;
for ith_poly = 1:Npolys
    Nradii = Nradii+length(polytopes(ith_poly).radii);
end

Nverticies_per_poly = 20; % Maximum estimate
NrealVertices = 0; % Actual number of vertices
Ntotal_vertices = Npolys*Nverticies_per_poly;
all_walls_start = nan(Ntotal_vertices,3);
all_walls_end   = nan(Ntotal_vertices,3);
all_angles = nan(Nverticies_per_poly,Npolys);
all_lengths = nan(Nverticies_per_poly,Npolys);
all_side_count = zeros(Npolys,1);
all_max_radius = zeros(Npolys,1);
all_min_radius = zeros(Npolys,1);
all_radii = zeros(Nradii,1);
all_mean_radius = zeros(Npolys,1);
all_sharpness_ratios = zeros(Npolys,1);
all_areas = zeros(Npolys,1);

all_x_vertices = extractfield(polytopes,'xv');
all_y_vertices = extractfield(polytopes,'yv');
% Update the AABB to fit vertices (make it bigger/smaller as needed)
AABB(1) = min(all_x_vertices);
AABB(2) = min(all_y_vertices);
AABB(3) = max(all_x_vertices);
AABB(4) = max(all_y_vertices);

% Loop through the polytopes
Nradii_found = 0;
for ith_poly = 1:Npolys
    vertices = polytopes(ith_poly).vertices;
    Nvertices = length(vertices(:,1))-1;
    NrealVertices = NrealVertices + Nvertices;

    % Save the start/end wall locations created by the polytope walls
    offset_row = (ith_poly-1)*Nverticies_per_poly;
    rows_to_fill = (1:Nvertices) + offset_row;
    all_walls_start(rows_to_fill,1:2) = vertices(1:end-1,:);
    all_walls_start(rows_to_fill,3) = ith_poly;
    all_walls_end(rows_to_fill,1:2) = vertices(2:end,:);
    all_walls_end(rows_to_fill,3) = ith_poly;



    % Determine number of verticies
    Nangles = length(vertices(:,1))-1;

    % Calculate the angles, noting that the function gives outside angles,
    % so we have to subtract them from 180 degrees to get the interior
    % angles
    [angles] = fcn_MapGen_polytopeFindVertexAngles(vertices, -1);
    all_angles(1:Nangles,ith_poly)  = pi - angles;

    % Calculate the other stats
    all_lengths(1:Nangles,ith_poly) = polytopes(ith_poly).distances;
    all_side_count(ith_poly,1) = Nangles;
    all_max_radius(ith_poly,1) = polytopes(ith_poly).max_radius;
    all_min_radius(ith_poly,1) = polytopes(ith_poly).min_radius;
    all_min_radius(ith_poly,1) = polytopes(ith_poly).min_radius;
    all_mean_radius(ith_poly,1) = polytopes(ith_poly).mean_radius;

    % Add up radii
    Nthisradii = length(polytopes(ith_poly).radii);
    all_radii((Nradii_found+1):Nradii_found+Nthisradii,1) = polytopes(ith_poly).radii;
    Nradii_found = Nradii_found+Nthisradii;


    % note that sharpness is the ratio of max radius to min radius
    % this is not the same thing as aspect ratio which is AABB width to height
    all_sharpness_ratios(ith_poly,1) = ...
        polytopes(ith_poly).max_radius/polytopes(ith_poly).min_radius;
    all_areas(ith_poly,1) = polytopes(ith_poly).area;

    % % Plot the input polytopes in red
    % % fcn_MapGen_OLD_plotPolytopes(polytopes(ith_poly),fig_num,'r',2);
    % plotFormat.LineWidth = 2;
    % plotFormat.MarkerSize = 10;
    % plotFormat.LineStyle = '-';
    % plotFormat.Color = [1 0 0];
    % fillFormat = [];
    % h_plot = fcn_MapGen_plotPolytopes(polytopes(ith_poly), (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

end

% Remove the NaN values from all_walls variables
all_walls_start_no_nan = all_walls_start(~isnan(all_walls_start(:,1)),:);
all_walls_end_no_nan   = all_walls_end(~isnan(all_walls_end(:,1)),:);

% Now, reshape the angles into one column, remove the NaN values, and
% calculate means and standar deviations
angle_column = reshape(all_angles,Nverticies_per_poly*Npolys,1);
% angle_column = reshape(all_angles,1,[]); % this line is a useful alternative if the above line fails
angle_column_no_nan = angle_column(~isnan(angle_column));
average_vertex_angle = mean(angle_column_no_nan*180/pi,'omitmissing');
std_vertex_angle = std(angle_column_no_nan*180/pi,'omitmissing');

% Find the mean and std deviations of radii metrics
average_max_radius = mean(all_max_radius,'omitmissing');
average_min_radius = mean(all_min_radius,'omitmissing');
average_mean_radius = mean(all_mean_radius,'omitmissing');
average_radius = mean(all_radii,'omitmissing');
average_sharpness = mean(all_sharpness_ratios,'omitmissing');
std_max_radius = std(all_max_radius,'omitmissing');
% all_mean_radii = extractfield(polytopes,'true_mean_radius');
% average_mean_radius = nanmean(all_mean_radii);

% Determine the length properties related to sides of polytopes
length_column = reshape(all_lengths,Nverticies_per_poly*Npolys,1);
% length_column = reshape(all_lengths,1,[]); % this line is a useful alternative if the above fails
length_column_no_nan = length_column(~isnan(length_column));
total_perimeter = sum(length_column_no_nan);
average_side_length = mean(length_column_no_nan,'omitmissing');
std_side_length = std(length_column_no_nan,'omitmissing');
average_perimeter = total_perimeter/Npolys;

% Determine the area properties of the map
occupied_area = sum(all_areas);
total_area    = (AABB(3)-AABB(1))*(AABB(4)-AABB(2));

% Determine the density properties
point_density = length(polytopes)/total_area;
linear_density = point_density.^0.5;


polyMapStats.Npolys = Npolys;
polyMapStats.NtotalVertices = NrealVertices;
avg_r_D = average_max_radius*linear_density;

% AREA METRICS
polyMapStats.min_x = AABB(1);
polyMapStats.max_x = AABB(3);
polyMapStats.min_y = AABB(2);
polyMapStats.max_y = AABB(4);

polyMapStats.occupied_area = occupied_area;
polyMapStats.total_area = total_area;
polyMapStats.point_density = point_density;
polyMapStats.linear_density = linear_density;

polyMapStats.Nangles = length(angle_column_no_nan(:,1));
polyMapStats.average_vertex_angle = average_vertex_angle;
polyMapStats.std_vertex_angle = std_vertex_angle;
polyMapStats.angle_column_no_nan = angle_column_no_nan;

polyMapStats.average_max_radius = average_max_radius;
polyMapStats.average_min_radius = average_min_radius;
polyMapStats.average_mean_radius = average_mean_radius;
polyMapStats.average_radius = average_radius;
polyMapStats.all_average_radius = all_mean_radius;
polyMapStats.all_radii = all_radii;
polyMapStats.average_sharpness = average_sharpness;
polyMapStats.std_max_radius = std_max_radius;
polyMapStats.average_side_length = average_side_length;
polyMapStats.std_side_length = std_side_length;
polyMapStats.total_perimeter = total_perimeter;
polyMapStats.avg_r_D = avg_r_D;
polyMapStats.NtotalVertices = NrealVertices;
polyMapStats.average_perimeter = average_perimeter;


if flag_do_debug
    figure(fig_for_debug);
    clf;
    hold on
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);

    % plot the polytopes
    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_for_debug,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_for_debug)); %#ok<NASGU>

    % plot all vertices
    plot(all_walls_start_no_nan(:,1),all_walls_start_no_nan(:,2),'c.','Markersize',10);

    % plot the means in black
    temp = zeros(length(polytopes),2);
    for ith_poly = 1:length(polytopes)
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);

    % number the polytopes at seed points
    for ith_poly = 1:length(polytopes)
        text_location = polytopes(ith_poly).mean;
        text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    end
end


% Check to see how often a line through the polytopes hits the same
% polytopes around it (e.g. the effective size of polytopes)

range_y = AABB(4) - AABB(2);
random_y_values_to_try = rand(100,1)*range_y + AABB(2);
delta_y_range = linspace(-average_max_radius*2,average_max_radius*2,100)';
ratio_overlaps = 0*delta_y_range;

% Find some y values to test, randomly
good_indices = (random_y_values_to_try>(AABB(2)+average_max_radius*2)) & ...
    (random_y_values_to_try<(AABB(4)-average_max_radius*2));
random_y_values_to_try = random_y_values_to_try(good_indices);
Ntries = length(random_y_values_to_try);

% Loop through the y-values
for ith_try = 1:Ntries
    % Find overlap for this y-case
    this_y = random_y_values_to_try(ith_try);
    this_overlap = fcn_INTERNAL_findRatioOverlaps(...
        this_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan,...
        delta_y_range);

    % Average results
    ratio_overlaps = ratio_overlaps + (1/Ntries)*this_overlap;
end

%% Find the experimental linear density
y_search_range = linspace(...
    AABB(2)+average_max_radius*2,...
    AABB(4)-average_max_radius*2,...
    10)';
line_crossing_hits = fcn_INTERNAL_findLinearDensityStats(...
    y_search_range, AABB,all_walls_start_no_nan,all_walls_end_no_nan);
polyMapStats.linear_density_mean = mean(line_crossing_hits);
polyMapStats.linear_density_std = std(line_crossing_hits);


%% Plot results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plot

    fprintf(1,'\n\nSUMMARY STATISTICS: \n');
    fprintf(1,'\t Number of polytopes: %.0d\n',polyMapStats.Npolys);
    fprintf(1,'\t Number of vertices: %.0d\n',polyMapStats.NtotalVertices);

    fprintf(1,'\n\t AREA METRICS:\n');
    fprintf(1,'\t Min_x: %.2f\n',polyMapStats.min_x);
    fprintf(1,'\t Max_x: %.2f\n',polyMapStats.max_x);
    fprintf(1,'\t Min_y: %.2f\n',polyMapStats.min_y);
    fprintf(1,'\t Max_y: %.2f\n',polyMapStats.max_y);
    fprintf(1,'\t Occupied Area: %.2f\n',occupied_area);
    fprintf(1,'\t Total Area: %.2f\n',total_area);
    fprintf(1,'\t Point density: %.2f\n',point_density);
    fprintf(1,'\t Theoretical Linear density: %.2f\n',linear_density);
    fprintf(1,'\t Calculated Linear density: %.2f +/- %.4f\n',polyMapStats.linear_density_mean,2*polyMapStats.linear_density_std);

    fprintf(1,'\n\t ANGLE METRICS:\n');
    fprintf(1,'\t Total number of vertices: %.2f\n',polyMapStats.Nangles);
    fprintf(1,'\t Average vertex angle (deg): %.2f\n',average_vertex_angle);
    fprintf(1,'\t Std dev vertex angle (deg): %.2f\n',std_vertex_angle);

    fprintf(1,'\n\t LENGTH METRICS:\n');
    fprintf(1,'\t Average maximum radius: %.4f\n',average_max_radius);
    fprintf(1,'\t Std dev maximum radius: %.4f\n',std_max_radius);
    fprintf(1,'\t Average side length: %.4f\n',average_side_length);
    fprintf(1,'\t Std dev side length: %.4f\n',std_side_length);
    fprintf(1,'\t Total perimeter: %.4f\n',total_perimeter);



    figure(fig_num)

    subplot(2,3,1);
    axis equal
    grid on;
    grid minor;

    % Fill in the x and y data
    polytope_plot_data_x = [];
    polytope_plot_data_y = [];
    for polys = 1:size(polytopes,2) % plot each polytope
        polytope_plot_data_x = [polytope_plot_data_x; polytopes(polys).vertices(:,1); nan]; %#ok<AGROW>
        polytope_plot_data_y = [polytope_plot_data_y; polytopes(polys).vertices(:,2); nan]; %#ok<AGROW>
    end
    % Plot polytopes
    plot(polytope_plot_data_x,polytope_plot_data_y,'-','Linewidth',2)
    title('Polytopes being analyzed');


    subplot(2,3,2);
    histogram(angle_column_no_nan*180/pi,100);
    xlabel('Angles (deg)');
    title(sprintf('Histogram of angles. Mean: %.2f, StdDev: %.2f',average_vertex_angle,std_vertex_angle));

    subplot(2,3,3);
    histogram(all_side_count,100);
    xlabel('N vertices (count)');
    title(sprintf('Histogram of side count. Mean: %.4f, StdDev: %.4f',mean(all_side_count,'omitmissing'),std(all_side_count,'omitmissing')));

    subplot(2,3,4);
    histogram(length_column_no_nan,100);
    xlabel('Side Lengths');
    title(sprintf('Histogram of side length. Mean: %.4f, StdDev: %.4f',average_side_length,std_side_length));

    subplot(2,3,5);
    histogram(all_max_radius,100);
    xlabel('Max radius');
    title(sprintf('Histogram of max radius. Mean: %.4f, StdDev: %.4f',average_max_radius,std_max_radius));

    subplot(2,3,6);
    hold on;
    plot(delta_y_range,ratio_overlaps,'k-');
    xline(-average_max_radius,'r-')
    xline(average_max_radius,'r-')
    xlabel('Range of deviation');
    title('% Similarity of polytopes hit versus range of deviation');
    legend('Similarity','+/- ave max radius');

    figure;
    box on;
    subplot(2,2,1);
    histogram(all_radii);
    title(sprintf("Distribution of All Polytope Radii.\nMean: %.4f, Median: %.4f, Mode: %.4f",mean(all_radii),mode(all_radii),mode(all_radii)));

    subplot(2,2,2);
    box on;
    histogram(all_min_radius);
    title(sprintf("Distribution of Min Polytope Radii.\nMean: %.4f, Median: %.4f, Mode: %.4f",mean(all_min_radius),mode(all_min_radius),mode(all_min_radius)));

    subplot(2,2,3);
    box on;
    histogram(all_max_radius);
    title(sprintf("Distribution of Max Polytope Radii.\nMean: %.4f, Median: %.4f, Mode: %.4f",mean(all_max_radius),mode(all_max_radius),mode(all_max_radius)));

    subplot(2,2,4);
    box on;
    histogram(all_mean_radius);
    title(sprintf("Distribution of Mean Radius per Polytope.\nMean: %.4f, Median: %.4f, Mode: %.4f",mean(all_mean_radius),mode(all_mean_radius),mode(all_mean_radius)));
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function line_crossing_hits = fcn_INTERNAL_findLinearDensityStats(...
    y_search_range, AABB,walls_start,walls_end)

flag_do_debug = 0;
fig_for_debug = 8456;  %#ok<NASGU>

if flag_do_debug
    figure(fig_for_debug);
    clf
    grid on;
    hold on;

    % Fill in the x and y data to plot the polytopes
    Nwalls = length(walls_start(:,1));
    polytope_plot_data_x = zeros(3*length(walls_start(:,1)),1);
    polytope_plot_data_y = zeros(3*length(walls_start(:,1)),1);
    for ith_wall = 1:Nwalls % plot each wall
        row_range = (ith_wall-1)*3 + (1:3);
        polytope_plot_data_x(row_range,1) = [walls_start(ith_wall,1); walls_end(ith_wall,1); nan];
        polytope_plot_data_y(row_range,1) = [walls_start(ith_wall,2); walls_end(ith_wall,2); nan];
    end

    % Plot polytopes
    plot(polytope_plot_data_x,polytope_plot_data_y,'-','Linewidth',2)
    title('Polytopes being analyzed for linear density statistics');

end

line_crossing_hits = zeros(length(y_search_range(:,1)),1);
for ith_y = 1:length(y_search_range)
    this_y = y_search_range(ith_y);
    polytopes_hit = fcn_INTERNAL_findPolysCrossedAtY(...
        this_y,AABB,walls_start,walls_end);
    line_crossing_hits(ith_y,1) = length(polytopes_hit(:,1));

    if flag_do_debug
        figure(fig_for_debug);
        yline(this_y,'r-');
        text(AABB(3),this_y,sprintf('%.0d hits',line_crossing_hits(ith_y,1)));
    end
end


end

function polytopes_hit = fcn_INTERNAL_findPolysCrossedAtY(...
    y_value,AABB,walls_start,walls_end)



width = AABB(3)-AABB(1);
bisector_start = [AABB(1)-0.1*width, y_value];
bisector_end = [AABB(3)+0.1*width,  y_value];

[~,~,walls_that_were_hit] = ...
        fcn_Path_findSensorHitOnWall(...
        walls_start(:,1:2),...   % wall start
        walls_end(:,1:2),...     % wall end
        bisector_start,...       % sensor_vector_start
        bisector_end,...         % sensor_vector_end
        (1), ...                 % (flag_search_return_type) -- 1 means ALL hits of any results,
        (0), ...                 % (flag_search_range_type)  -- 0 means only if overlapping wall/sensor, ...
        ([]),...                 % (tolerance) -- default is eps * 1000,
        (-1));                   % (fig_num) -- -1 means to use "fast mode")

polytopes_hit = unique(walls_start(walls_that_were_hit,3));
end

function ratio_overlaps = fcn_INTERNAL_findRatioOverlaps(...
    midpoint_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan,...
    delta_y_range)

flag_do_debug = 0;
fig_for_debug = 999; %#ok<NASGU>

original_polytopes_hit = fcn_INTERNAL_findPolysCrossedAtY(...
    midpoint_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan);
original_num_hit = length(original_polytopes_hit);

if flag_do_debug
    figure(fig_for_debug);
    yline(midpoint_y,...
        'k-','Linewidth',3);
end


ratio_overlaps = zeros(length(delta_y_range(:,1)),1);
for ith_y = 1:length(delta_y_range)
    this_y = delta_y_range(ith_y)+midpoint_y;
    new_polytopes_hit = fcn_INTERNAL_findPolysCrossedAtY(...
        this_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan);
    number_overlap = sum(ismember(new_polytopes_hit,original_polytopes_hit));
    ratio_overlaps(ith_y) = number_overlap/original_num_hit;
end
end
