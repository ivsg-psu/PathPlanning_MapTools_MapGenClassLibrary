function poly_map_stats = fcn_MapGen_polytopesStatistics(...
    polytopes,...
    varargin)
% fcn_MapGen_polytopesStatistics calculates key statistics for the
% polytopes
%
% FORMAT:
%
% [poly_map_stats] = ...
%     fcn_MapGen_polytopesStatistics(...
%     polytopes,...
%     (fig_num))
%
% INPUTS:
%
%     POLYTOPES: polytopes to analyze
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
% OUTPUTS:
%
%     poly_map_stats: a structure whose fields contain the polytope map's
%     statistics
%     each of the fields in the struct are as follows:
%     AABB is the axis aligned bounding box for the polytope field, bounding all polytopes
%
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_polytopeFindVertexAngles
%     fcn_MapGen_plotPolytopes
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

% TO DO
% -- TBD

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 9453;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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

if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');

end


% Does user want to show the plots?
if  2 == nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
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
all_radii = [];
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
    [angles] = fcn_MapGen_polytopeFindVertexAngles(vertices);
    all_angles(1:Nangles,ith_poly)  = pi - angles;

    % Calculate the other stats
    all_lengths(1:Nangles,ith_poly) = polytopes(ith_poly).distances;
    all_side_count(ith_poly,1) = Nangles;
    all_max_radius(ith_poly,1) = polytopes(ith_poly).max_radius;
    all_min_radius(ith_poly,1) = polytopes(ith_poly).min_radius;
    all_min_radius(ith_poly,1) = polytopes(ith_poly).min_radius;
    all_mean_radius(ith_poly,1) = polytopes(ith_poly).mean_radius;
    all_radii = [all_radii; polytopes(ith_poly).radii];
    % note that sharpness is the ratio of max radius to min radius
    % this is not the same thing as aspect ratio which is AABB width to height
    all_sharpness_ratios(ith_poly,1) = ...
        polytopes(ith_poly).max_radius/polytopes(ith_poly).min_radius;
    all_areas(ith_poly,1) = polytopes(ith_poly).area;

    % Plot the input polytopes in red
    %fcn_MapGen_plotPolytopes(polytopes(ith_poly),fig_num,'r',2);

end

% Remove the NaN values from all_walls variables
all_walls_start_no_nan = all_walls_start(~isnan(all_walls_start(:,1)),:);
all_walls_end_no_nan   = all_walls_end(~isnan(all_walls_end(:,1)),:);

% Now, reshape the angles into one column, remove the NaN values, and
% calculate means and standar deviations
% angle_column = reshape(all_angles,Nverticies_per_poly*Npolys,1);
angle_column = all_angles;
angle_column_no_nan = all_angles;
% angle_column_no_nan = angle_column(~isnan(angle_column));
average_vertex_angle = nanmean(angle_column_no_nan*180/pi);
std_vertex_angle = nanstd(angle_column_no_nan*180/pi);

% Find the mean and std deviations of radii metrics
average_max_radius = nanmean(all_max_radius);
average_min_radius = nanmean(all_min_radius);
average_mean_radius = nanmean(all_mean_radius);
average_radius = nanmean(all_radii);
average_sharpness = nanmean(all_sharpness_ratios);
std_max_radius = nanstd(all_max_radius);
% all_mean_radii = extractfield(polytopes,'true_mean_radius');
% average_mean_radius = nanmean(all_mean_radii);

% Determine the length properties related to sides of polytopes
% length_column = reshape(all_lengths,Nverticies_per_poly*Npolys,1);
length_column = all_lengths;
length_column_no_nan = all_lengths;
% length_column_no_nan = length_column(~isnan(length_column));
total_perimeter = sum(length_column_no_nan,'omitnan');
average_side_length = nanmean(length_column_no_nan);
std_side_length = nanstd(length_column_no_nan);
average_perimeter = total_perimeter/Npolys;

% Determine the area properties of the map
occupied_area = sum(all_areas);
total_area    = (AABB(3)-AABB(1))*(AABB(4)-AABB(2));

% Determine the density properties
point_density = length(polytopes)/total_area;
linear_density = point_density.^0.5;


poly_map_stats.Npolys = Npolys;
poly_map_stats.NtotalVertices = NrealVertices;
avg_r_D = average_max_radius*linear_density;

% AREA METRICS
poly_map_stats.min_x = AABB(1);
poly_map_stats.max_x = AABB(3);
poly_map_stats.min_y = AABB(2);
poly_map_stats.max_y = AABB(4);

poly_map_stats.occupied_area = occupied_area;
poly_map_stats.total_area = total_area;
poly_map_stats.point_density = point_density;
poly_map_stats.linear_density = linear_density;

poly_map_stats.Nangles = length(angle_column_no_nan(:,1));
poly_map_stats.average_vertex_angle = average_vertex_angle;
poly_map_stats.std_vertex_angle = std_vertex_angle;
poly_map_stats.angle_column_no_nan = angle_column_no_nan;

poly_map_stats.average_max_radius = average_max_radius;
poly_map_stats.average_min_radius = average_min_radius;
poly_map_stats.average_mean_radius = average_mean_radius;
poly_map_stats.average_radius = average_radius;
poly_map_stats.all_average_radius = all_mean_radius;
poly_map_stats.all_radii = all_radii;
poly_map_stats.average_sharpness = average_sharpness;
poly_map_stats.std_max_radius = std_max_radius;
poly_map_stats.average_side_length = average_side_length;
poly_map_stats.std_side_length = std_side_length;
poly_map_stats.total_perimeter = total_perimeter;
poly_map_stats.avg_r_D = avg_r_D;
poly_map_stats.NtotalVertices = NrealVertices;
poly_map_stats.average_perimeter = average_perimeter;


if flag_do_debug
    figure(fig_for_debug);
    clf;
    hold on
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);

    % plot the polytopes
    fcn_MapGen_plotPolytopes(polytopes,fig_for_debug,'b',2);

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
    this_overlap = INTERNAL_fcn_findRatioOverlaps(...
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
line_crossing_hits = INTERNAL_fcn_findLinearDensityStats(...
    y_search_range, AABB,all_walls_start_no_nan,all_walls_end_no_nan);
poly_map_stats.linear_density_mean = mean(line_crossing_hits);
poly_map_stats.linear_density_std = std(line_crossing_hits);


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
    fprintf(1,'\t Number of polytopes: %.0d\n',poly_map_stats.Npolys);
    fprintf(1,'\t Number of vertices: %.0d\n',poly_map_stats.NtotalVertices);

    fprintf(1,'\n\t AREA METRICS:\n');
    fprintf(1,'\t Min_x: %.2f\n',poly_map_stats.min_x);
    fprintf(1,'\t Max_x: %.2f\n',poly_map_stats.max_x);
    fprintf(1,'\t Min_y: %.2f\n',poly_map_stats.min_y);
    fprintf(1,'\t Max_y: %.2f\n',poly_map_stats.max_y);
    fprintf(1,'\t Occupied Area: %.2f\n',occupied_area);
    fprintf(1,'\t Total Area: %.2f\n',total_area);
    fprintf(1,'\t Point density: %.2f\n',point_density);
    fprintf(1,'\t Theoretical Linear density: %.2f\n',linear_density);
    fprintf(1,'\t Calculated Linear density: %.2f +/- %.4f\n',poly_map_stats.linear_density_mean,2*poly_map_stats.linear_density_std);

    fprintf(1,'\n\t ANGLE METRICS:\n');
    fprintf(1,'\t Total number of vertices: %.2f\n',poly_map_stats.Nangles);
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
    title(sprintf('Histogram of side count. Mean: %.4f, StdDev: %.4f',nanmean(all_side_count),nanstd(all_side_count)));

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

function line_crossing_hits = INTERNAL_fcn_findLinearDensityStats(...
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
    polytopes_hit = INTERNAL_fcn_findPolysCrossedAtY(...
        this_y,AABB,walls_start,walls_end);
    line_crossing_hits(ith_y,1) = length(polytopes_hit(:,1));

    if flag_do_debug
        figure(fig_for_debug);
        yline(this_y,'r-');
        text(AABB(3),this_y,sprintf('%.0d hits',line_crossing_hits(ith_y,1)));
    end
end


end

function polytopes_hit = INTERNAL_fcn_findPolysCrossedAtY(...
    y_value,AABB,walls_start,walls_end)



width = AABB(3)-AABB(1);
bisector_start = [AABB(1)-0.1*width, y_value];
bisector_end = [AABB(3)+0.1*width,  y_value];

[~,~,walls_that_were_hit] = ...
    fcn_MapGen_findIntersectionOfSegments(...
    walls_start(:,1:2),...
    walls_end(:,1:2),...
    bisector_start,...
    bisector_end,...
    2);
polytopes_hit = unique(walls_start(walls_that_were_hit,3));
end

function ratio_overlaps = INTERNAL_fcn_findRatioOverlaps(...
    midpoint_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan,...
    delta_y_range)

flag_do_debug = 0;
fig_for_debug = 999; %#ok<NASGU>

original_polytopes_hit = INTERNAL_fcn_findPolysCrossedAtY(...
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
    new_polytopes_hit = INTERNAL_fcn_findPolysCrossedAtY(...
        this_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan);
    number_overlap = sum(ismember(new_polytopes_hit,original_polytopes_hit));
    ratio_overlaps(ith_y) = number_overlap/original_num_hit;
end
end
