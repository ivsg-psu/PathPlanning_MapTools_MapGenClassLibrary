function [ ...
filled_polytopes ...
] = ...
fcn_MapGen_fillPolytopeFieldsFromVertices( ...
polytopes, ...
varargin...
)
% fcn_MapGen_fillPolytopeFieldsFromVertices
% Given a polytoope structure array where the vertices field is filled,
% calculates the values for all the other fields.
%
%
%
% FORMAT:
%
%    [ ...
%    filled_polytopes ...
%    ] = ...
%    fcn_MapGen_fillPolytopeFieldsFromVertices( ...
%    polytopes, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     polytopes: an individual structure or structure array of 'polytopes'
%     type that stores the polytopes to be filled
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
%
%
% OUTPUTS:
%
%     filled_polytopes: the polytopes array with all fields completed
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_polytopeCentroidAndArea
%     fcn_MapGen_checkInputsToFunctions
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_fillPolytopeFieldsFromVertices
% for a full test suite.
%
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of function

%
% TO DO:
%
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 1;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 3;
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


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it has vertices
    if ~isfield(polytopes,'vertices')
        error('Field of vertices was not found');
    end

    % Check the vertices input to have 4 or more rows, 2 columns
    %     fcn_MapGen_checkInputsToFunctions(...
    %         polytopes.vertices, '2column_of_numbers',[4 5]);


end

% Does user want to show the plots?
if  2== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
        flag_do_plot = 1;
    end
end

%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

% Initialize variables
filled_polytopes = polytopes;
num_poly = length(polytopes);

fig_idx = 1;

for ith_poly = 1:num_poly % pull each polytope
    % adjust polytopes
    filled_polytopes(ith_poly).xv        = (polytopes(ith_poly).vertices(1:end-1,1)');
    filled_polytopes(ith_poly).yv        = (polytopes(ith_poly).vertices(1:end-1,2)');
    filled_polytopes(ith_poly).distances = ...
        sum((polytopes(ith_poly).vertices(1:end-1,:) - ...
        polytopes(ith_poly).vertices(2:end,:)).^2,2).^0.5;

    % Calculate the mean and area
    [filled_polytopes(ith_poly).mean,filled_polytopes(ith_poly).area] = ...
        fcn_MapGen_polytopeCentroidAndArea(polytopes(ith_poly).vertices);

    % Find max radius
    radii = sum(...
        (filled_polytopes(ith_poly).vertices(1:end-1,:) - ...
        ones(length(filled_polytopes(ith_poly).xv),1)*filled_polytopes(ith_poly).mean).^2,2).^0.5;
    filled_polytopes(ith_poly).max_radius = ...
        max(radii);
    filled_polytopes(ith_poly).min_radius = ...
        min(radii);
    filled_polytopes(ith_poly).mean_radius = ...
        mean(radii);
    filled_polytopes(ith_poly).radii = radii;
    filled_polytopes(ith_poly).cost = rand;

    %% find average radius from calculus
    % gather information for this polytope
    poly = filled_polytopes(ith_poly);
    centroid = poly.mean;
    n_sides = length(poly.xv);
    % translate polytope so centroid is at (0,0)
    vertices_recentered = [poly.vertices(:,1)-centroid(1),poly.vertices(:,2)-centroid(2)];
    if flag_do_plot
        % plot original polytope
        figure(12345321)
        clf;
        plot(1000*poly.vertices(:,1),1000*poly.vertices(:,2))
        hold on
        plot(centroid(1),centroid(2),'kd')
        figure(12345322)
        clf;
       % plot recentered sides
        plot(1000*vertices_recentered(:,1),1000*vertices_recentered(:,2))
        hold on
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
            theta_range1 = ceil(theta1):0.25:359;
            theta_range2 = 0:0.25:floor(theta2);
            theta_range = [theta_range1, theta_range2];
        else
            theta_range = ceil(theta1):0.25:floor(theta2);
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
            xlabel('x [m]')
            ylabel('y [m]')
            title('polytope, plotted from Cartesian and Polar coordinates')
            legend('Cartesian','centroid','','','polar')
            savefig(sprintf('poly_rect_and_polar_%d',fig_idx))
        end
    end
    if flag_do_plot
        figure(12345323)
        clf;
        % plot r(theta) for whole polytope on rectangular axes
        plot(theta_range_ordered,1000*r_of_theta_all_sides,'kx')
        xlabel('\theta [deg]')
        ylabel('R (\theta) [m]')
        title('Radius, R, as a function of angle from centroid to wall, \theta, for one polytope')
        savefig(sprintf('r_of_theta_%d',fig_idx))
    end
    % the average distance from centroid to sides is the average value of r(theta) from 0 to 360
    filled_polytopes(ith_poly).true_mean_radius = mean(r_of_theta_all_sides);
    filled_polytopes(ith_poly).radii_dist = r_of_theta_all_sides;
    filled_polytopes(ith_poly).theta_dist = r_of_theta_all_sides;
    if flag_do_plot
        figure(12345324)
        clf;
        % plot histogram (pdf) of single polytope
        h = histogram(1000*r_of_theta_all_sides)
        xlabel('radius value, R [m]')
        ylabel('count')
        title(sprintf('Histogram of Polytope Radius Values,\n measured from a single polytope'))
        savefig(sprintf('r_hist_%d',fig_idx))
        h.BinWidth = 0.2;
        savefig(sprintf('r_hist_narrow_%d',fig_idx))
        figure(12345325)
        clf;
        h2 = histogram(1000*r_of_theta_all_sides,'Normalization','cdf');
        [N,edges] = histcounts(1000*r_of_theta_all_sides,'Normalization','cdf')
        bin_center = (edges(1:end-1)+edges(2:end))/2;
        [N,edges] = histcounts(1000*r_of_theta_all_sides,'Normalization','pdf')
        h2.BinWidth = 0.2;
        bin_edges = 0:0.1:20;
        [N,edges] = histcounts(1000*r_of_theta_all_sides,bin_edges,'Normalization','cdf')
        bin_center = (edges(1:end-1)+edges(2:end))/2;
        plot(bin_center,1-N)
        title(sprintf('Empirical CDF of Polytope Radius Values,\n measured from a single polytope'))
        xlabel('radius value, R [m]')
        ylabel('P(R_0<=R)')
        savefig(sprintf('r_cdf_narrow_%d',fig_idx))
        figure(12345326)
        clf;
        h3 = histogram(1000*r_of_theta_all_sides,'Normalization','pdf');
        h3.BinWidth = 0.2;
        title(sprintf('Empirical PDF of Polytope Radius Values,\n measured from a single polytope'))
        xlabel('radius value, R [m]')
        ylabel('P(R_0=R)')
        savefig(sprintf('r_pdf_narrow_%d',fig_idx))
        fig_idx = fig_idx + 1;
    end
end

r_of_theta_all_polytopes = extractfield(filled_polytopes,'radii_dist');
if flag_do_plot
    % plot histogram (pdf) of all polytopes
    figure(12345327)
    clf;
    histogram(1000*r_of_theta_all_polytopes)
    savefig(sprintf('all_r_hist'))
end
save('workspace')

%§
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
    figure(fig_num);
    hold on

    % plot the polytopes
    fcn_MapGen_plotPolytopes(filled_polytopes,fig_num,'b',2);

    % plot the means in black
    temp = zeros(length(filled_polytopes),2);
    for ith_poly = 1:length(filled_polytopes)
        temp(ith_poly,:) = filled_polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);
end % Ends the flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end


end % Ends the function

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



