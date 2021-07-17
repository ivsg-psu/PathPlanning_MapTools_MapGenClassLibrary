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
%
% DEPENDENCIES:
% 
%     fcn_MapGen_checkInputsToFunctions
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
Nangles = Npolys*Nverticies_per_poly;
all_angles = nan(Nverticies_per_poly,Npolys);
all_lengths = nan(Nverticies_per_poly,Npolys);
all_side_count = zeros(Npolys,1);
all_max_radius = zeros(Npolys,1);
all_areas = zeros(Npolys,1);
AABB = [inf inf -inf -inf];

% Loop through the polytopes
for ith_poly = 1:Npolys
    vertices = polytopes(ith_poly).vertices;

    % Update the AABB
    AABB(1) = min(AABB(1),min(vertices(:,1)));
    AABB(2) = min(AABB(2),min(vertices(:,2)));
    AABB(3) = max(AABB(3),max(vertices(:,1)));
    AABB(4) = max(AABB(4),max(vertices(:,2)));

    
    Nangles = length(vertices(:,1))-1;
    calculation_vertices = [vertices; vertices(2,:)];
    
    starting_vector = calculation_vertices(2:end-1,:)- calculation_vertices(1:end-2,:);
    ending_vector = calculation_vertices(3:end,:)- calculation_vertices(2:end-1,:);
    
    unit_starting_vector = starting_vector./(sum(starting_vector.^2,2).^0.5);
    unit_ending_vector   = ending_vector./(sum(ending_vector.^2,2).^0.5);
   
    cross_products = cross([unit_starting_vector zeros(Nangles,1)],[unit_ending_vector zeros(Nangles,1)]);
    cross_result = cross_products(:,3);
    
    angles_cross = asin(cross_result);
    dot_products = dot(unit_starting_vector,unit_ending_vector,2);
    angles_dot = acos(dot_products);
    
    all_angles(1:Nangles,ith_poly)  = angles_dot.*sign(angles_cross);
    all_lengths(1:Nangles,ith_poly) = polytopes(ith_poly).distances;
    all_side_count(ith_poly,1) = Nangles;
    all_max_radius(ith_poly,1) = polytopes(ith_poly).max_radius;
    all_areas(ith_poly,1) = polytopes(ith_poly).area;
        
    % Plot the input polytopes in red
    %fcn_MapGen_plotPolytopes(polytopes(ith_poly),fig_num,'r',2);
    
end
angle_column = reshape(all_angles,Nverticies_per_poly*Npolys,1);
angle_column_no_nan = angle_column(~isnan(angle_column));
average_vertex_angle = nanmean(angle_column_no_nan*180/pi);
std_vertex_angle = nanstd(angle_column_no_nan*180/pi);

average_max_radius = nanmean(all_max_radius);
std_max_radius = nanstd(all_max_radius);

length_column = reshape(all_lengths,Nverticies_per_poly*Npolys,1);
length_column_no_nan = length_column(~isnan(length_column));
total_perimeter = sum(length_column_no_nan);

average_side_length = nanmean(length_column_no_nan);
std_side_length = nanstd(length_column_no_nan);

occupied_area = sum(all_areas);
total_area    = (AABB(3)-AABB(1))*(AABB(4)-AABB(2));
unoccupied_area = total_area-occupied_area;
unoccupancy_ratio = (total_area - occupied_area)/total_area;

point_density = length(polytopes)/total_area;
linear_density = point_density.^0.5;

average_gap_size_G_bar = (unoccupancy_ratio/point_density)^0.5; % See Eq. (4.24 in Seth Tau's thesis)
perimeter_gap_size = 2*unoccupied_area/(total_perimeter);

% Fill in results
poly_map_stats.occupied_area = occupied_area;
poly_map_stats.total_area = total_area;
poly_map_stats.unoccupied_area = unoccupied_area;
poly_map_stats.unoccupancy_ratio = unoccupancy_ratio;
poly_map_stats.point_density = point_density;
poly_map_stats.linear_density = linear_density;
poly_map_stats.average_gap_size_G_bar = average_gap_size_G_bar;
poly_map_stats.perimeter_gap_size = perimeter_gap_size;
poly_map_stats.average_vertex_angle = average_vertex_angle;
poly_map_stats.std_vertex_angle = std_vertex_angle;
poly_map_stats.average_max_radius = average_max_radius;
poly_map_stats.std_max_radius = std_max_radius;
poly_map_stats.average_side_length = average_side_length;
poly_map_stats.std_side_length = std_side_length;
poly_map_stats.total_perimeter = total_perimeter;


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
    fprintf(1,'\tOccupied Area: %.2f\n',occupied_area);
    fprintf(1,'\tTotal Area: %.2f\n',total_area);
    fprintf(1,'\tUnoccupied Area: %.2f\n',unoccupied_area);
    fprintf(1,'\tUnoccupancy ratio: %.2f\n',unoccupancy_ratio);
    fprintf(1,'\tPoint density: %.2f\n',point_density);
    fprintf(1,'\tLinear density: %.2f\n',linear_density);
    fprintf(1,'\tAverage gap size, G-bar: %.5f\n',average_gap_size_G_bar);
    fprintf(1,'\tPerimeter gap size: %.5f\n',perimeter_gap_size);
    fprintf(1,'\tAverage vertex angle (deg): %.2f\n',average_vertex_angle);
    fprintf(1,'\tStd dev vertex angle (deg): %.2f\n',std_vertex_angle);
    fprintf(1,'\tAverage maximum radius: %.4f\n',average_max_radius);
    fprintf(1,'\tStd dev maximum radius: %.4f\n',std_max_radius);
    fprintf(1,'\tAverage side length: %.4f\n',average_side_length);
    fprintf(1,'\tStd dev side length: %.4f\n',std_side_length);
    fprintf(1,'\tTotal perimeter: %.4f\n',total_perimeter);
    
    
    
    figure(fig_num)
    
    subplot(2,3,1);   
    % Fill in the x and y data
    polytope_plot_data_x = [];
    polytope_plot_data_y = [];
    for polys = 1:size(polytopes,2) % plot each polytope
        polytope_plot_data_x = [polytope_plot_data_x; polytopes(polys).vertices(:,1); nan]; %#ok<AGROW>
        polytope_plot_data_y = [polytope_plot_data_y; polytopes(polys).vertices(:,2); nan]; %#ok<AGROW>
    end
    % Plot polytopes
    plot(polytope_plot_data_x,polytope_plot_data_y,'-','Linewidth',2)
    
    
    
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


