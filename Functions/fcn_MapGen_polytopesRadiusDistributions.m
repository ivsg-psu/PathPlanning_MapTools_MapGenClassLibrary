function poly_size_stats = fcn_MapGen_polytopesRadiusDistributions(polytopes, varargin)
% fcn_MapGen_polytopesRadiusDistributions uses statistics and calculaus to get
% estimates of polytope radius distribution, expected radius, and expected
% width, independent of traversal location and direction
% by rotating polytopes about their centroids to produce a radial probability of
% occupation
%
% FORMAT:
% poly_size_stats = fcn_MapGen_polytopesRadiusDistributions(polytopes, (fig_num))
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
%     poly_size_stats: a structure whose fields contain radii distributions and
%     expected values for polytopes in the field
%
% DEPENDENCIES:
%
%     fcn_MapGen_polytopesStatistics
%
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopesRadiusDistributions
%
% This function was written in 07_2022 by Steve Harnett
% Questions or comments? sjharnett@psu.edu
%
% Revision History:
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(1,2);

        % % Check the polytopes input, make sure it is 'polytopes' type
        % fcn_DebugTools_checkInputsToFunctions(...
        %     polytopes, 'polytopes');
        % 
        % 
        % % Check the exp_dist input, make sure it is 'positive_column_of_numbers' type
        % fcn_DebugTools_checkInputsToFunctions(...
        %     exp_dist, 'positive_1column_of_numbers',1);

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (2 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; %#ok<NASGU>
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
fig_idx = 1;
% r vs theta curve for all polys (vectors per poly)
poly_size_stats.r_of_theta_all_polys = [];
% probability of occupation for all polys (vectors per poly)
poly_size_stats.N_all_polys = [];
poly_size_stats.bin_edges_all_polys = [];
% expected radius values for all polys (scalar per poly)
poly_size_stats.expected_radii = [];
% effective depths as a function of offset for all polys (vectors per poly)
poly_size_stats.offsets = [];
poly_size_stats.effective_depths = [];

poly_size_stats.effective_depth_scalars = [];

poly_size_stats.d_eff_i_at_middle_o = [];

poly_size_stats.a_d_eff_i_eval_at_o_i = [];
poly_size_stats.a_d_eff_i_eval_at_o_avg = [];

% max radius of all polys in the field
max_radius_field = max(extractfield(polytopes,'max_radius'));
field_stats = fcn_MapGen_polytopesStatistics(polytopes);
field_avg_r_D = field_stats.avg_r_D;
N_int = field_stats.linear_density_mean;
D_L = field_avg_r_D;
D_L_i_avg = D_L/N_int*1000;
o_avg = mean(extractfield(polytopes,'max_radius'))*1000-D_L_i_avg;
if o_avg < 0
    o_avg = 0;
end

for ith_poly = 1:length(polytopes)

    %% find average radius from calculus
    % gather information for this polytope
    poly = polytopes(ith_poly);
    centroid = poly.mean;
    if poly.area == 0
        continue
    end
    n_sides = length(poly.xv);
    delta_theta = 1;
    % translate polytope so centroid is at (0,0)
    vertices_recentered = [poly.vertices(:,1)-centroid(1),poly.vertices(:,2)-centroid(2)];
    if flag_do_plot
        % plot original polytope
        figure(12345321)
        clf;
        plot(1000*poly.vertices(:,1),1000*poly.vertices(:,2))
        box on
        hold on
        plot(centroid(1),centroid(2),'kd')
        figure(12345322)
        clf;
        % plot recentered sides
        plot(1000*vertices_recentered(:,1),1000*vertices_recentered(:,2))
        hold on
        box on
        % plot recentered centroid
        plot(0,0,'kd')
    end
    % initialize array for storing r(theta) function for polytope sides
    r_of_theta_all_sides = [];
    % althought we know the range is 0 to 360 we don't know where it staarts and ends
    theta_range_ordered = [];
    % loop through all sides on this polytope
    for ith_side = 1:n_sides
        % need to find the equation for a line representing this side
        x1 = vertices_recentered(ith_side,1);
        x2 = vertices_recentered(ith_side+1,1);
        y1 = vertices_recentered(ith_side,2);
        y2 = vertices_recentered(ith_side+1,2);
        % find slope of this line
        m = (y2-y1)/(x2-x1);
        % find the angle from the centroid to the first and second vertecies
        theta1 = atan2d(y1,x1);
        theta2 = atan2d(y2,x2);
        if theta1 < 0
            theta1 = theta1 + 360;
        end
        if theta2 < 0
            theta2 = theta2 + 360;
        end
        % check to see if this is the point crossing 0 degrees
        if theta1 > theta2
            % if it is we want to go to theta1 to 0, then from 0 to theta2
            theta_range1 = ceil(theta1):delta_theta:359;
            theta_range2 = 0:delta_theta:floor(theta2);
            theta_range = [theta_range1, theta_range2];
        else
            theta_range = ceil(theta1):delta_theta:floor(theta2);
        end
        % convert our y(x) formula for a line to r(theta), then find r for the appropriate theta
        r_of_theta = (-1*m*x2+y2)./(sind(theta_range)-m*cosd(theta_range));
        % add the r(theta) data for this side to the vector for the whole shape
        r_of_theta_all_sides = [r_of_theta_all_sides, r_of_theta];
        theta_range_ordered = [theta_range_ordered, theta_range];
        if flag_do_plot
            figure(12345322)
            % plot vertices of interest
            plot(1000*x1,1000*y1,'ro')
            plot(1000*x2,1000*y2,'ro')
            % plot polytope to check that polar form of line forming side matches rectangular
            plot(1000*r_of_theta.*cosd(theta_range),1000*r_of_theta.*sind(theta_range),'kx')
            box on
            xlabel('x [m]')
            ylabel('y [m]')
            title('polytope, plotted from Cartesian and Polar coordinates')
            legend('Cartesian','centroid','','','polar')
        end
    end
    if flag_do_plot
        figure(12345323)
        clf;
        % plot r(theta) for whole polytope on rectangular axes
        box on
        plot(theta_range_ordered,1000*r_of_theta_all_sides,'kx')
        xlabel('\theta [deg]')
        ylabel('R (\theta) [m]')
        title('Radius, R, as a function of angle from centroid to wall, \theta, for one polytope')
    end
    % the average distance from centroid to sides is the average value of r(theta) from 0 to 360
    % order the r array by theta
    r_and_theta_all_sides = [theta_range_ordered; r_of_theta_all_sides];
    [~,idx] = sort(r_and_theta_all_sides(1,:));
    r_and_theta_all_sides = r_and_theta_all_sides(:,idx);
    % now we can store the r-curve for this all polys, knowing they all are ordered
    % from 0 degrees to 359 with 1 degree increments
%     poly_size_stats.r_of_theta_all_polys = [poly_size_stats.r_of_theta_all_polys; r_and_theta_all_sides(2,:)];
    if flag_do_plot
        figure(12345324)
        clf;
        % plot histogram of single polytope
        h = histogram(1000*r_of_theta_all_sides)
        box on
        xlabel('radius value, R [m]')
        ylabel('count')
        title(sprintf('Histogram of Polytope Radius Values,\n measured from a single polytope'))
        h.BinWidth = 0.2;

        % plot PDF of single polytope
        figure(12345325)
        clf;
        h2 = histogram(1000*r_of_theta_all_sides,'Normalization','pdf');
        h2.BinWidth = 0.2;
        box on
        title(sprintf('Empirical PDF of Polytope Radius Values,\n measured from a single polytope'))
        xlabel('radius value, R [m]')
        ylabel('P(R_0 = R)')

        % plot CDF of single polytope
        figure(123453255)
        clf;
        h3 = histogram(1000*r_of_theta_all_sides,'Normalization','cdf');
        h3.BinWidth = 0.2;
        box on
        title(sprintf('Empirical CDF of Polytope Radius Values,\n measured from a single polytope'))
        xlabel('radius value, R [m]')
        ylabel('P(R_0 \geq R)')
    end
    % instead of probability of radius having that value or less
    % max bin should be factor of safety of 1.2 over max_radius in m
    bin_edges = 0:0.1:1000*max_radius_field*1.2;
    if max_radius_field == 0
        sprintf('field average max radius is 0')
        break
    end
    % plot 1-CDF of single polytope (probability of occupation
    [N,edges] = histcounts(1000*r_of_theta_all_sides,bin_edges,'Normalization','cdf');
    bin_center = (edges(1:end-1)+edges(2:end))/2;
    poly_size_stats.N_all_polys = [poly_size_stats.N_all_polys; 1-N];
    poly_size_stats.bin_edges_all_polys = [poly_size_stats.bin_edges_all_polys; bin_center];
    if flag_do_plot
        % plot 1-CDF of a single polytope - this is probability of occupation
        figure(123453256)
        clf;
        plot(bin_center,1-N)
        box on
        title(sprintf('Probability of occupation at a given radius\n measured from a single polytope'))
        xlabel('radius value, R [m]')
        ylabel('P that space is occupied by obstacle at R')
    end
    % find average depth from integrating P(occupied)
    expected_radius = 0;
    % for number of bins...
    for i = 1:1:(length(bin_center)-1)
        % expected value is P(R)*delta_R
        delta_R = bin_center(i+1)-bin_center(i);
        % add to total length the current length interval, weighted by
        % its probability
        expected_radius = ...
            expected_radius + ...
            (1-N(i))*delta_R;
    end
    poly_size_stats.expected_radii = [poly_size_stats.expected_radii, expected_radius];
    if flag_do_plot
        % show expected radius on probability of occupation plot
        figure(123453256)
        hold on;
        box on
        plot([expected_radius, expected_radius],[-0.1,1.1],'-k')
        leg2 = sprintf('expected value of R = %.2f',expected_radius);
        legend('probability of occupation',leg2)
        xlim([0,1000*max_radius_field*1.1])
        ylim([-0.05,1.05])
    end
    o_this_poly = max(r_of_theta_all_sides)*1000-D_L_i_avg;
    if o_this_poly < 0
        o_this_poly = 0;
    end
    %% find effective depth for d at some offset, o
    % create offset range, to R_max
    o = 0:0.2:1000*max_radius_field*1.2; % TODO
    effective_depths = [];
    % loop through all possible offets (i.e. strike positions)
    for j = 1:1:length(o)
        offset = o(j);
        % initialize effective depth
        effective_depth = 0;
        % integrate over all probabilities
        delta_d = 0.1;
        d_max = sqrt(max(1000*max(r_of_theta_all_sides))^2-offset^2); % TODO
        for i = 0:delta_d:d_max
            % expected value is P(R)*delta_R
            % delta_R = bin_center(i+1)-bin_center(i);
            current_d = i;%delta_d*(i-1);
            R_at_current_d = sqrt(current_d^2+offset^2);
            % the r we're at may not be have an associated probability
            % so let's find the closest R we have prob for
            [~,closest_index] = min(abs(bin_center-R_at_current_d));
            % now find prob for this R
            prob_occ_at_R = 1-N(closest_index);
            effective_depth = ...
                effective_depth + ...
                prob_occ_at_R*delta_d;
        end
        effective_depths = [effective_depths,effective_depth];
    end

    % find the closest o to the estiamted o for this poly
    [~,closest_index] = min(abs(o-o_this_poly));
    d_eff_i_eval_at_o_i = effective_depths(closest_index);
    poly_size_stats.a_d_eff_i_eval_at_o_i = [poly_size_stats.a_d_eff_i_eval_at_o_i, d_eff_i_eval_at_o_i];

%     middle_o_this_poly = o(floor(length(effective_depths(effective_depths>0.0001))/2)); % TODO
%     [~,closest_index] = min(abs(o-middle_o_this_poly));
%     poly_size_stats.d_eff_i_at_middle_o = [poly_size_stats.d_eff_i_at_middle_o, effective_depths(closest_index)];
    % find the closest o to the estimate o for all polys
    [~,closest_index] = min(abs(o-o_avg));
    d_eff_i_eval_at_o_avg = effective_depths(closest_index);
    poly_size_stats.a_d_eff_i_eval_at_o_avg = [poly_size_stats.a_d_eff_i_eval_at_o_avg, d_eff_i_eval_at_o_avg];

    poly_size_stats.offsets = [poly_size_stats.offsets; o];
    poly_size_stats.effective_depths = [poly_size_stats.effective_depths;effective_depths];
    poly_size_stats.effective_depth_scalars = [poly_size_stats.effective_depth_scalars,mean(effective_depths(effective_depths>0.0001))];
    if flag_do_plot
        figure(12345327)
        hold on
        box on
        plot(o,effective_depths)
        xlabel('offset from centroid, o [m]')
        ylabel('effective depth, d_{eff} [m]')
        title('Effective depth of a polytope as a function of offset, d_{eff}(o)')
    end
% end loop through all polys
end
mean_d_eff = mean(poly_size_stats.effective_depths,1);
poly_size_stats.mean_d_eff = mean_d_eff;

% find the closest o to the estimate o for all polys
[~,closest_index] = min(abs(o-o_avg));
d_eff_avg_eval_at_o_avg = mean_d_eff(closest_index);
poly_size_stats.d_eff_avg_eval_at_o_avg = d_eff_avg_eval_at_o_avg;

if flag_do_plot
    plot(o,mean_d_eff,'LineWidth',3,'Color','k')
    mean(poly_size_stats.effective_depths(5,:),1);
    num_curves = [1,5,20,50,100,500];
    for i = 1:1:length(num_curves)
        figure(12345329)
        clf
        hold on
        box on
        avg_curve = mean(poly_size_stats.effective_depths(1:num_curves(i),:),1);
        if flag_do_plot
            plot(o,avg_curve,'LineWidth',3,'Color','k')
            xlabel('offset from centroid, o [m]')
            ylabel('effective depth, d_{eff} [m]')
            title(sprintf('Effective depth of a polytope as a function of offset, d_{eff}(o) from %.0d curves',num_curves(i)))
            savefig(sprintf('d_eff_%d',i))
        end
    end
end
% find average of average scalar value of average d_eff curve
poly_size_stats.mean_d_eff_scalar = mean(mean_d_eff(mean_d_eff>0.0001));
% middle_o_avg_curve = o(floor(length(mean_d_eff(mean_d_eff>0.0001))/2));
% [~,closest_index] = min(abs(o-middle_o_avg_curve));
% poly_size_stats.mean_d_eff_at_middle_o = mean_d_eff(closest_index);
if flag_do_plot
    figure(12345327)
    hold on
    plot([0,max(o)],[poly_size_stats.mean_d_eff_scalar,poly_size_stats.mean_d_eff_scalar],'k-')
end

%% Plot the results (for debugging)?
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
    figure(fig_num)
    clf;
    % 
    % LineWidth = 2;
    % fcn_MapGen_plotPolytopes(polytopes,fig_num,'r-',LineWidth);
    % fcn_MapGen_plotPolytopes(exp_polytopes,fig_num,'b-',LineWidth,'square');
    % legend('Original','Expanded')
    % box on
    % xlabel('X Position')
    % ylabel('Y Position')

end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends main function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§


