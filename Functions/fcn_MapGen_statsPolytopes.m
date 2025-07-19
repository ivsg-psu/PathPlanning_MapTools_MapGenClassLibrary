function polyMapStats = fcn_MapGen_statsPolytopes(polytopes, varargin)
% fcn_MapGen_statsPolytopes calculates key statistics for the
% polytopes
%
% FORMAT:
%
%     polyMapStats = fcn_MapGen_statsPolytopes( polytopes, (fig_num))
%
% INPUTS:
%
%     POLYTOPES: polytopes to analyze
%
%     (optional inputs)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     polyMapStats: a structure whose fields contain the polytope map's
%     statistics. A summary copy is below:
%         SUMMARY STATISTICS: 
%         	 Number of polytopes: 21
%         	 Number of vertices: 110
%         
%         	 AREA METRICS:
%         	 Min_x: 0.00
%         	 Max_x: 1.00
%         	 Min_y: 0.00
%         	 Max_y: 1.00
%         	 Occupied Area: 1.00
%         	 Total Area: 1.00
%         	 Point density: 21.00
%         	 Theoretical Linear density: 4.58
%         
%         	 ANGLE METRICS:
%         	 Total number of vertices: 110.00
%         	 Average vertex angle (deg): 111.27
%         	 Std dev vertex angle (deg): 32.00
%         
%         	 LENGTH METRICS:
%         	 Average maximum radius: 0.2030
%         	 Std dev maximum radius: 0.0354
%         	 Average side length: 0.1788
%         	 Std dev side length: 0.1141
%         	 Total perimeter: 19.6705
%         
%         	 POLYTOPE ENCOUNTER METRICS:
%         	 Max similar path width (non-zero width of ratio_overlaps): 0.6043
%         	 Theoretical Linear density: 4.58
%         	 Experimental_linear_density mean/2sigma: 5.9000 +/- 1.4757
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
% For additional examples, see: script_test_fcn_MapGen_statsPolytopes
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
% 2025_07_17 by Sean Brennan
% -- standardized Debugging and Input checks area, Inputs area
% -- made codes use MAX_NARGIN definition at top of code, narginchk
% -- made plotting flag_do_plots and code consistent across all functions
% -- cleanup of header docstrings
% -- added debugging visualization to see what is happening

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
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
    debug_fig_num = 999978; 
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
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(1,MAX_NARGIN);

        % Check the polytopes input
        fcn_DebugTools_checkInputsToFunctions(polytopes, 'polytopes');

    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
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
if flag_do_debug
    % Prep the figure
    figure(debug_fig_num);
    clf;
    hold on;

    % Plot the input polytopes in black
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (debug_fig_num)); 
    set(h_plot,'DisplayName','Input: polytopes');    

    % Set up to plot individual polytopes
    h_poly = fcn_MapGen_plotPolytopes(polytopes(1), (plotFormat), (fillFormat), (debug_fig_num));
    set(h_poly,'DisplayName','current poly');    

    legend('Interpreter','none','Location','best');

end

Npolys = length(polytopes);

% Set up matrices, all_angles and all_lengths, that store in each row, the
% angles and lengths for each polytope, where each polytope is a column
Nverticies_per_poly = 20; % Maximum estimate
NrealVertices   = 0; % Actual number of vertices
Ntotal_vertices = Npolys*Nverticies_per_poly;
all_walls_start = nan(Ntotal_vertices,3);
all_walls_end   = nan(Ntotal_vertices,3);
all_angles      = nan(Nverticies_per_poly,Npolys);
all_lengths     = nan(Nverticies_per_poly,Npolys);
all_side_count  = zeros(Npolys,1);

% Fill in values that are easily extracted
all_radii = extractfield(polytopes,'radii');
all_x_vertices  = extractfield(polytopes,'xv');
all_y_vertices  = extractfield(polytopes,'yv');
all_max_radius  = extractfield(polytopes,'max_radius');
all_min_radius  = extractfield(polytopes,'min_radius');
all_mean_radius = extractfield(polytopes,'mean_radius');
all_areas       = extractfield(polytopes,'area');
all_sharpness_ratios = all_max_radius./all_min_radius;

% Update the AABB to fit vertices (make it bigger/smaller as needed)
AABB(1) = min(all_x_vertices);
AABB(2) = min(all_y_vertices);
AABB(3) = max(all_x_vertices);
AABB(4) = max(all_y_vertices);

% Loop through the polytopes to fill in:
% - angles in each polytope
% - where side "walls" start/stop
% - number of sides in the polytope
allPlottingVertices = [];
for ith_poly = 1:Npolys
    % Grab this polytope
    thisPoly = polytopes(ith_poly);

    % Show what poly we are looking at by plotting it in red
    if flag_do_debug
        % Remove previous poly from legend
        h_poly.('HandleVisibility') = 'off'; % Take previous poly off legend
        h_poly.Color = 0.8*[1 1 1]; % Dim previous poly

        % Plot the current polytope in red
        plotFormat.LineWidth = 2;
        plotFormat.MarkerSize = 10;
        plotFormat.LineStyle = '-';
        plotFormat.Color = [1 0 0];
        plotFormat.DisplayName = 'current poly';
        fillFormat = [];
        h_poly = fcn_MapGen_plotPolytopes(polytopes(ith_poly), (plotFormat), (fillFormat), (debug_fig_num));
        pause(0.01);
    end
    
    % Grab the vertices
    vertices = thisPoly.vertices;
    Nvertices = length(vertices(:,1))-1;
    NrealVertices = NrealVertices + Nvertices;

    % Save for plotting?
    if 1==flag_do_plots
        allPlottingVertices = [allPlottingVertices; vertices; nan nan]; %#ok<AGROW>
    end

    % Save the start/end wall locations created by the polytope walls
    offset_row = (ith_poly-1)*Nverticies_per_poly;
    rows_to_fill = (1:Nvertices) + offset_row;
    all_walls_start(rows_to_fill,1:2) = vertices(1:end-1,:);
    all_walls_start(rows_to_fill,3) = ith_poly;
    all_walls_end(rows_to_fill,1:2) = vertices(2:end,:);
    all_walls_end(rows_to_fill,3) = ith_poly;

    % Determine number of sides. This is # vertices minus 1
    Nsides = length(vertices(:,1))-1;
    all_side_count(ith_poly,1) = Nsides;

    % Calculate the angles, noting that the function gives outside angles,
    % so we have to subtract them from 180 degrees to get the interior
    % angles
    [angles] = fcn_MapGen_polytopeFindVertexAngles(vertices, -1);
    all_angles(1:Nsides,ith_poly)  = pi - angles;

    % Save lengths for use in perimeter calculations later
    all_lengths(1:Nsides,ith_poly) = thisPoly.distances;
end

% Remove the NaN values from all_walls variables
all_walls_start_no_nan = all_walls_start(~isnan(all_walls_start(:,1)),:);
all_walls_end_no_nan   = all_walls_end(~isnan(all_walls_end(:,1)),:);

% Now, reshape the angles into one column, remove the NaN values, and
% calculate means and standard deviations
angle_column = reshape(all_angles,Nverticies_per_poly*Npolys,1);
% angle_column = reshape(all_angles,1,[]); % this line is a useful alternative if the above line fails
angle_column_no_nan = angle_column(~isnan(angle_column));
average_vertex_angle = mean(angle_column_no_nan*180/pi,'omitmissing');
std_vertex_angle = std(angle_column_no_nan*180/pi,'omitmissing');

% Find the mean and std deviations of radii metrics
max_max_radius = max(all_max_radius,[],'all','omitmissing');
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
point_density = Npolys/total_area;
linear_density = point_density.^0.5;
Nangles = length(angle_column_no_nan(:,1));

polyMapStats.Npolys = Npolys;
polyMapStats.NtotalVertices = NrealVertices;
avg_r_D = average_max_radius*linear_density;


%% Find ratio_overlaps
% Check to see how often a line through the polytopes hits the same
% polytopes around it (e.g. the effective size of polytopes)

% Sample the map randomly, NrandomRepeats times, at different y values in
% the map. At each location, find the ratio_overlaps value and add these up
% to generate the total ratio_overlaps.
NrandomRepeats = 10;
minY = AABB(2)+2*max_max_radius; % Have to go by 2 times because the search delta_y_range goes this far below and above.
maxY = AABB(4)-2*max_max_radius;
range_y = maxY - minY;
random_y_values_to_try = rand(NrandomRepeats,1)*range_y + minY;

delta_y_range = linspace(-max_max_radius*2,max_max_radius*2,31)';
ratio_overlaps = 0*delta_y_range;

% Loop through the y-values
for ith_try = 1:NrandomRepeats
    % Find overlap for this y-case
    this_y = random_y_values_to_try(ith_try);
    this_overlap = fcn_INTERNAL_findRatioOverlaps(...
        this_y, AABB,all_walls_start_no_nan,all_walls_end_no_nan,...
        delta_y_range);

    % Average results
    ratio_overlaps = ratio_overlaps + (1/NrandomRepeats)*this_overlap;
    if this_overlap(16)~=1 || this_overlap(1)~=0 || this_overlap(end)~=0
        error('Unexpected results found in overlap calculations between adjacent paths.');
    end
end
ratio_overlaps_ySampling = delta_y_range;

% For debugging
if 1==0
    figure(47474);
    clf;
    plot(ratio_overlaps_ySampling,ratio_overlaps)
end

firstIndex = find(ratio_overlaps~=0,1,'first');
lastIndex = find(ratio_overlaps~=0,1,'last');
maxSimilarPathWidth = delta_y_range(lastIndex) - delta_y_range(firstIndex);

%% Find the experimental linear density
% Samples at 10 points in y direction, evenly spaced

y_search_range = linspace(...
    AABB(2)+average_max_radius*2,...
    AABB(4)-average_max_radius*2,...
    10)';
line_crossing_hits = fcn_INTERNAL_findLinearDensityStats(...
    y_search_range, AABB, all_walls_start_no_nan,all_walls_end_no_nan);

experimental_linear_density = line_crossing_hits;
experimental_linear_density_ySampling = y_search_range;
experimental_linear_density_mean = mean(line_crossing_hits);
experimental_linear_density_std = std(line_crossing_hits);


%% Save metrics
polyMapStats.AABB = AABB; % The axis-aligned bounding box for the polytopes
polyMapStats.min_x = AABB(1); % Minimum of all x values
polyMapStats.max_x = AABB(3); % Maximum of all x values
polyMapStats.min_y = AABB(2); % Minimum of all y values
polyMapStats.max_y = AABB(4); % Maximum of all x values
polyMapStats.occupied_area = occupied_area;   % Sum of all areas of all polytopes
polyMapStats.total_area = total_area;         % Total area within the AABB
polyMapStats.point_density = point_density;   % Number of polytopes per unit area in AABB
polyMapStats.linear_density = linear_density; % Square root of point density
polyMapStats.Nangles = Nangles; % Number of angles within polytopes
polyMapStats.average_vertex_angle = average_vertex_angle; % Average vertex angle (interior angle)
polyMapStats.std_vertex_angle = std_vertex_angle;  % Standard deviation of vertex angle
polyMapStats.angle_column_no_nan = angle_column_no_nan; % Column of all angles
polyMapStats.average_max_radius = average_max_radius; % Mean of max radius of each poly
polyMapStats.average_min_radius = average_min_radius; % Mean of min radius of each poly
polyMapStats.average_mean_radius = average_mean_radius; % Mean of mean radius of each poly
polyMapStats.average_radius = average_radius; % Average radius across all polys
polyMapStats.all_average_radius = all_mean_radius; % A column of all mean radii from each poly
polyMapStats.all_radii = all_radii; % A column of all radii from each poly
polyMapStats.average_sharpness = average_sharpness; % Average sharpness 
polyMapStats.std_max_radius = std_max_radius; % Std in max radius
polyMapStats.average_side_length = average_side_length; % Average side length
polyMapStats.std_side_length = std_side_length; % Std in side length
polyMapStats.total_perimeter = total_perimeter; % Total perimeter, e.g. sum of all side lengths
polyMapStats.avg_r_D = avg_r_D; % average_max_radius*linear_density (scalar)
polyMapStats.NtotalVertices = NrealVertices; % Total number of vertices in polys
polyMapStats.average_perimeter = average_perimeter; % Average perimeter of all polys
polyMapStats.ratio_overlaps = ratio_overlaps; % Measures similarity of paths to each other. 
polyMapStats.ratio_overlaps_ySampling = ratio_overlaps_ySampling; % The y-spacing from mean where ratio_overlaps is calculated
polyMapStats.experimental_linear_density = experimental_linear_density; % Linear density raw data
polyMapStats.experimental_linear_density_ySampling = experimental_linear_density_ySampling; % The y-spacing values for linear density
polyMapStats.experimental_linear_density_mean = experimental_linear_density_mean; % Ave number of polys encountered in straight-line pass across map in y-direction
polyMapStats.experimental_linear_density_std = experimental_linear_density_std; % Std deviation in poly count, of polys encountered in straight-line pass across map in y-direction
polyMapStats.maxSimilarPathWidth = maxSimilarPathWidth; % Maximum width wherein one y-traversal is similar to another.

if flag_do_debug
    figure(fig_for_debug);
    clf;
    hold on
    scale = max(AABB,[],'all') - min(AABB,[],'all');
    new_axis = [AABB(1)-scale/2 AABB(3)+scale/2 AABB(2)-scale/2 AABB(4)+scale/2];
    axis(new_axis);

    % plot the polytopes
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_for_debug));

    % plot all vertices
    plot(all_walls_start_no_nan(:,1),all_walls_start_no_nan(:,2),'c.','Markersize',10);

    % plot the means in black
    temp = zeros(Npolys,2);
    for ith_poly = 1:Npolys
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);

    % number the polytopes at seed points
    for ith_poly = 1:Npolys
        text_location = polytopes(ith_poly).mean;
        text(text_location(1,1),text_location(1,2),sprintf('%.0d',ith_poly));
    end
end

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

if flag_do_plots

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
    fprintf(1,'\t Theoretical Linear density: %.2f\n',polyMapStats.linear_density);

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

    fprintf(1,'\n\t POLYTOPE ENCOUNTER METRICS:\n');
    fprintf(1,'\t Max similar path width (non-zero width of ratio_overlaps): %.4f\n',maxSimilarPathWidth);
    fprintf(1,'\t Theoretical Linear density: %.2f\n',polyMapStats.linear_density);
    fprintf(1,'\t Experimental_linear_density mean/2sigma: %.4f +/- %.4f\n',experimental_linear_density_mean,2*experimental_linear_density_std);

    figure(fig_num)

    subplot(2,3,1);
    axis equal
    grid on;
    grid minor;

    % Plot polytopes
    plot(allPlottingVertices(:,1),allPlottingVertices(:,2),'-','Linewidth',2)
    title('Polytopes being analyzed');

    subplot(2,3,2);
    histogram(angle_column_no_nan*180/pi,36);
    ave_value = average_vertex_angle;
    std_value = std_vertex_angle;
    xline(ave_value,'g-')
    xline(ave_value-2*std_value,'r-')
    xline(ave_value+2*std_value,'r-')
    text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
    text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
    text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
    xlabel('Angles (deg)');
    title('Histogram of angles');

    subplot(2,3,3);
    histogram(all_side_count,10);
    ave_value = mean(all_side_count,'omitmissing');
    std_value = std(all_side_count,'omitmissing');
    xline(ave_value,'g-')
    xline(ave_value-2*std_value,'r-')
    xline(ave_value+2*std_value,'r-')
    text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
    text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
    text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
    xlabel('N vertices (count)');
    title('Histogram of side count');

    subplot(2,3,4);
    histogram(length_column_no_nan,20);
    ave_value = average_side_length;
    std_value = std_side_length;
    xline(ave_value,'g-')
    xline(ave_value-2*std_value,'r-')
    xline(ave_value+2*std_value,'r-')
    text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
    text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
    text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
    xlabel('Side Lengths');
    title('Histogram of side length');

    subplot(2,3,5);
    histogram(all_max_radius,20);
    ave_value = average_max_radius;
    std_value = std_max_radius;
    xline(ave_value,'g-')
    xline(ave_value-2*std_value,'r-')
    xline(ave_value+2*std_value,'r-')
    text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
    text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
    text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
    xlabel('Max radius [m]');
    title('Histogram of max radius'); 

    subplot(2,3,6);
    hold on;
    plot(delta_y_range,ratio_overlaps,'k-');
    xline(-average_max_radius,'g-')
    xline(average_max_radius,'g-')
    xline(-maxSimilarPathWidth/2,'r-')
    xline(maxSimilarPathWidth/2,'r-')
    text(-average_max_radius,0,'- Ave. Max radius','Rotation',90,'Color',[0 1 0]);
    text(average_max_radius,0,'+ Ave. Max radius','Rotation',90,'Color',[0 1 0]);
    text(-maxSimilarPathWidth/2,0,'-maxSimilarPathWidth/2,','Rotation',90,'Color',[1 0 0]);
    text(maxSimilarPathWidth/2,0,'maxSimilarPathWidth/2,','Rotation',90,'Color',[1 0 0]);
    xlabel('Range of deviation [m]');
    title('% Similarity of polytopes hit versus range of deviation');

    if flag_do_debug
        figure;
        box on;
        subplot(2,2,1);
        histogram(all_radii);
        ave_value = mean(all_radii);
        std_value = std(all_radii);
        xline(ave_value,'g-')
        xline(ave_value-2*std_value,'r-')
        xline(ave_value+2*std_value,'r-')
        text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
        text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
        text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
        title('Distribution of All Polytope Radii');

        subplot(2,2,2);
        box on;
        histogram(all_min_radius);
        ave_value = mean(all_min_radius);
        std_value = std(all_min_radius);
        xline(ave_value,'g-')
        xline(ave_value-2*std_value,'r-')
        xline(ave_value+2*std_value,'r-')
        text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
        text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
        text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
        title('Distribution of Min Polytope Radii');

        subplot(2,2,3);
        box on;
        histogram(all_max_radius);
        ave_value = mean(all_max_radius);
        std_value = std(all_max_radius);
        xline(ave_value,'g-')
        xline(ave_value-2*std_value,'r-')
        xline(ave_value+2*std_value,'r-')
        text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
        text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
        text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
        title('Distribution of Max Polytope Radii');

        subplot(2,2,4);
        box on;
        histogram(all_mean_radius);
        ave_value = mean(all_mean_radius);
        std_value = std(all_mean_radius);
        xline(ave_value,'g-')
        xline(ave_value-2*std_value,'r-')
        xline(ave_value+2*std_value,'r-')
        text(ave_value,0,'Mean','Rotation',90,'Color',[0 1 0]);
        text(ave_value+2*std_value,0,'+2sigma','Rotation',90,'Color',[1 0 0]);
        text(ave_value-2*std_value,0,'-2sigma','Rotation',90,'Color',[1 0 0]);
        title('Distribution of Mean Radius per Polytope');
    end
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
% Finds how many polytopes are encountered for a straight line "ray" in
% y-direction across polytopes. Results are recorded for each value in
% y_search_range

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

%% fcn_INTERNAL_findPolysCrossedAtY
function polytopes_hit = fcn_INTERNAL_findPolysCrossedAtY(...
    y_value,AABB,walls_start,walls_end)
% Finds the polytopes that are passed through if a projection is made at a
% y-value offset. Saves the IDs of the polytopes that were hit.

width = AABB(3)-AABB(1);

% Start the bisector just a little bit outside the AABB in both the low and
% high x directions.
bisector_start = [AABB(1)-0.1*width, y_value];
bisector_end = [AABB(3)+0.1*width,  y_value];

% Find which walls were hit
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

% Keep only unique values. The polytope ID is saved in the 3rd column of
% the walls_start array.
if ~all(isnan(walls_that_were_hit))
    polytopes_hit = unique(walls_start(walls_that_were_hit,3));
else
    polytopes_hit = 0;
end

end % Ends fcn_INTERNAL_findPolysCrossedAtY

%% fcn_INTERNAL_findRatioOverlaps
function ratio_overlaps = fcn_INTERNAL_findRatioOverlaps(...
    midpoint_y,AABB,all_walls_start_no_nan,all_walls_end_no_nan,...
    delta_y_range)
% Given a midpoint_y value and the AABB, finds how many polytopes are hit
% at the midpoint by passing a "ray" through the AABB at the y value and
% recording each polytope index. It then searches through a delta_y_range
% and again shoots a ray value through, counting how many of the midpoint
% polytope list are still hit. This result is then saved as a ratio. Thus,
% a 100% ratio means that the ray went through exactly the same polys as
% the midpoint, whereas a 0% ratio means the ray did not encounter any
% hits. The ratio is therefore a good metric of how similar one path is to
% another adjacent path, if they are both going in the same direction but
% offset from each other. For maps to be reused at different areas, their
% ratios at different areas should be zero; otherwise, the results will be
% correlated to each other.

flag_do_debug = 0;
fig_for_debug = 999; %#ok<NASGU>

% Find which polytopes are hit at the midpoint y-value
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
        this_y,AABB, all_walls_start_no_nan, all_walls_end_no_nan);
    number_overlap = sum(ismember(new_polytopes_hit,original_polytopes_hit));
    ratio_overlaps(ith_y) = number_overlap/original_num_hit;
end
end % Ends fcn_INTERNAL_findRatioOverlaps
