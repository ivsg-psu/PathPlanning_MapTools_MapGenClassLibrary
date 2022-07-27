function poly_size_stats = fcn_MapGen_polytopesRadiusDistributions(polytopes)
% fcn_MapGen_polytopesRadiusDistributions uses statistics and calculaus to get
% estimates of polytope radius distribution, expected radius, and expected
% width, independent of traversal location and direction
%
% FORMAT:
% poly_size_stats = fcn_MapGen_polytopesRadiusDistributions(polytopes)
%
% INPUTS:
%
%     POLYTOPES: polytopes to analyze
%
% OUTPUTS:
%
%     poly_size_stats: a structure whose fields contain radii distributions and
%     expected values for polytopes in the field
%
% DEPENDENCIES:
%
%     none from IVSG libraries, only Matlab built-in depencies
%
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopesRadiusDistributions
%
% This function was written in 07_2022 by Steve Harnett
% Questions or comments? sjharnett@psu.edu
%
% Revision History:
%
% TO DO
%
%% Debugging and Input checks
flag_check_inputs = 0; % Set equal to 1 to check the input arguments
flag_do_plot = 1;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 9453;
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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
poly_size_stats.r_of_theta_all_polys = [];

for ith_poly = 1:length(polytopes)

    %% find average radius from calculus
    % gather information for this polytope
    poly = polytopes(ith_poly);
    centroid = poly.mean;
    n_sides = length(poly.xv);
    delta_theta = 1;
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
    % order the r array by theta
    r_and_theta_all_sides = [theta_range_ordered; r_of_theta_all_sides];
    [~,idx] = sort(r_and_theta_all_sides(1,:));
    r_and_theta_all_sides = r_and_theta_all_sides(:,idx);
    poly_size_stats.r_of_theta_all_polys = [poly_size_stats.r_of_theta_all_polys; r_and_theta_all_sides(2,:)];
    if flag_do_plot
        figure(12345324)
        clf;
        % plot histogram of single polytope
        h = histogram(1000*r_of_theta_all_sides)
        xlabel('radius value, R [m]')
        ylabel('count')
        title(sprintf('Histogram of Polytope Radius Values,\n measured from a single polytope'))
        savefig(sprintf('r_hist_%d',fig_idx))
        h.BinWidth = 0.2;
        savefig(sprintf('r_hist_narrow_%d',fig_idx))
        % plot CDF of single polytope
        figure(12345325)
        clf;
        h2 = histogram(1000*r_of_theta_all_sides,'Normalization','cdf');
        [N,edges] = histcounts(1000*r_of_theta_all_sides,'Normalization','cdf')
        bin_center = (edges(1:end-1)+edges(2:end))/2;
        [N,edges] = histcounts(1000*r_of_theta_all_sides,'Normalization','pdf')
        h2.BinWidth = 0.2;
        bin_edges = 0:0.1:20;
        % plot 1-CDF of single polytope (probability of occupation
        [N,edges] = histcounts(1000*r_of_theta_all_sides,bin_edges,'Normalization','cdf')
        bin_center = (edges(1:end-1)+edges(2:end))/2;
        plot(bin_center,1-N)
        % TODO @sjharnett
        % all_prob_occ = [all_prob_occ; 1-N];
        title(sprintf('Empirical CDF of Polytope Radius Values,\n measured from a single polytope'))
        xlabel('radius value, R [m]')
        ylabel('P(R_0 \geq R)')
        savefig(sprintf('r_cdf_narrow_%d',fig_idx))
        % find average depth from average of P(solid)
        expected_radius = 0;
        for i = 1:1:(length(bin_center)-1)
            % expected value is P(R)*delta_R
            delta_R = bin_center(i+1)-bin_center(i);
            expected_radius = ...
                expected_radius + ...
                (1-N(i))*delta_R;
        end
        hold on;
        plot([expected_radius, expected_radius],[-0.1,1.1],'-k')
        legend('probability of occupation','expected value of R')
        xlim([0,20])
        ylim([-0.05,1.05])
        hold on;
        plot([expected_radius, expected_radius],[-0.1,1.1],'-k')
        legend('probability of occupation','expected value of R')
        %% find effective depth for d at some offset, o
        % create offset range, to R_max
        % TODO @sjharnett see if you can make o range from 0 to a giant
        % number so it will align for all curves
        o = 0:0.2:max(1000*r_of_theta_all_sides);
        effective_depths = [];
        % loop through all possible offets (i.e. strike positions)
        for j = 1:1:length(o)
            offset = o(j);
            % initialize effective depth
            effective_depth = 0;
            % integrate over all probabilities
            % TODO @sjharnett should this be looping from 1 to d(o)?
            % want to go from 0 to the maximum d(o) which is r_max
            delta_d = 0.1;
            d_max = sqrt(max(1000*r_of_theta_all_sides)^2-offset^2);
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
            % TODO @sjharnett store this in all effective depth curves for
            % all shapes
        end
        figure(12345327)
        hold on
        plot(o,effective_depths)
        xlabel('offset from centroid, o [m]')
        ylabel('effective depth, d_{eff} [m]')
        title('Effective depth of a polytope as a function of offset, d(o)')
        % plot PDF of single polytope
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

%# r_of_theta_all_polytopes = extractfield(filled_polytopes,'radii_dist');
%# if flag_do_plot
%#     % plot histogram (pdf) of all polytopes
%#     figure(12345327)
%#     clf;
%#     histogram(1000*r_of_theta_all_polytopes)
%#     savefig(sprintf('all_r_hist'))
%# end
% TODO @sjharnett take average of effective depth curve
% also take average of CDF curve
% write code to save figure every few iterations (e.g. avg of 1 shape, 2
% shapes, 5, 10, 20)
% find average of average d_eff curve
% build this expected value of avg d_eff curve into predictions
save('workspace')
