function [shrunk_polytopes,mu_final,sigma_final] = ...
    fcn_MapGen_polytopesShrinkToRadius(...
    polytopes,...
    des_radius,...
    sigma_radius,...
    min_rad,...
    varargin)
% fcn_MapGen_polytopesShrinkToRadius shrinks the polytopes to achieve the
% desired mean radius and specified variance
%
% FORMAT:
% 
% [shrunk_polytopes,mu_final,sigma_final] = ...
%     fcn_MapGen_polytopesShrinkToRadius(...
%     polytopes,...
%     des_radius,...
%     sigma_radius,...
%     min_rad,...
%     (fig_num))
%
% INPUTS:
%
%     POLYTOPES: original polytopes with same fields as shrunk_polytopes
%
%     DES_RAD: desired average max radius   
%
%     SIGMA_RADIUS: desired variance in the radii 
%
%     MIN_RAD: minimum acceptable radius
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
%     SHRUNK_POLYTOPES: a 1-by-n seven field structure of shrunken polytopes, 
%     where n <= number of polytopes with fields:
%       vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
%       xv: a 1-by-m vector of vertice x-coordinates
%       yv: a 1-by-m vector of vertice y-coordinates
%       distances: a 1-by-m vector of perimeter distances from one point to the
%         next point, distances(i) = distance from vertices(i) to vertices(i+1)
%       mean: centroid xy coordinate of the polytope
%       area: area of the polytope
%       max_radius: distance from the mean to the farthest vertex
%
%     MU_FINAL: final average maximum radius achieved
%
%     SIGMA_FINAL: final variance achieved
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_plotPolytopes
%     fcn_MapGen_polytopeShrinkToRadius
% 
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopesShrinkToRadius
%
% This function was written on 2019_08_29 by Seth Tau
% Questions or comments? sat5340@psu.edu 
%

% Revision History:
% 2021-06-08 - S. Brennan
% -- revised function to prep for MapGen class 
% -- added plotting option
% -- added comments, added debugging option
% 2021-06-12
% -- added minimum radius check on inputs, throws warning if min radius
% less than zero.
% -- made warnings more clear on the truncations.
% 2021-06-16
% -- added more comments via Seth's inputs
% 2021-07-06
% -- added positive input checking to fcn_MapGen_polytopesShrinkToRadius
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions

% TO DO
% -- Vectorize the for loop if possible
% -- check inputs are positive numbers where appropriate (e.g. make a
% "positive number" check

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
    if flag_check_inputs
        % Are there the right number of inputs?
        if nargin < 4 || nargin > 5
            error('Incorrect number of input arguments')
        end

        % Check the polytopes input
        fcn_DebugTools_checkInputsToFunctions(...
            polytopes, 'polytopes');

        % Check the des_radius input
        fcn_DebugTools_checkInputsToFunctions(...
            des_radius, 'positive_1column_of_numbers',1);

        % Check the sigma_radius input
        fcn_DebugTools_checkInputsToFunctions(...
            sigma_radius, '1column_of_numbers',1);

        % Check the min_rad input
        fcn_DebugTools_checkInputsToFunctions(...
            min_rad, 'positive_1column_of_numbers',1);

    end
end
    

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  5 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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


%% find current distribution
old_max_radii = [polytopes.max_radius]';
Nradii = length(old_max_radii);
old_r_mu = mean(old_max_radii);
old_r_sigma = std(old_max_radii);

if flag_do_debug
    fprintf(1,'Target distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',des_radius);
    fprintf(1,'\tStd dev: %.4f\n',sigma_radius);
    
    fprintf(1,'Input distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',old_r_mu);
    fprintf(1,'\tStd dev: %.4f\n',old_r_sigma);
    
    fcn_MapGen_plotPolytopes(polytopes,fig_for_debug,'b',2);
    
    figure(fig_for_debug+1);
    histogram(old_max_radii,20)
    title(sprintf('Histogram of input radii. Mean: %.4f, Std-dev: %.4f. Targets are: %.4f and %.4f',...
        old_r_mu,old_r_sigma,...
        des_radius,...
        sigma_radius));
end


if old_r_mu < des_radius
    error('cannot achieve the desired radius by shrinking because average radius is already smaller than desired radius')
end

%% determine desired distribution
new_r_dist = normrnd(des_radius,sigma_radius,[Nradii,1]);

% adjust to ensure the mean value is mu. SETH: is this necessary?
new_r_dist = new_r_dist + (des_radius-mean(new_r_dist)); 


if flag_do_debug
    new_r_mu = mean(new_r_dist);
    new_r_sigma = std(new_r_dist);
    fprintf(1,'Ideal distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',new_r_mu);
    fprintf(1,'\tStd dev: %.4f\n',new_r_sigma);
    
    figure(fig_for_debug+2);
    histogram(new_r_dist,20)
    title(sprintf('Histogram of target radii. Mean: %.4f, Std-dev: %.4f. Targets are: %.4f and %.4f',...
        new_r_mu,new_r_sigma,...
        des_radius,...
        sigma_radius));
   
end

% Check to see if truncation will occur. THis is checked by using the
% definition of the maximum radius and minimum radius, and checking how
% many points fall outside of this range. If there are truncations, warn
% the user.

% % OLD: uses a single maximum for all polytopes
% max_r_dist = max(old_max_radii); % largest possible radius
% NEW: uses a each polytope's maximum. SETH: the above code gives wrong
% answers if there's a few polytopes that have very large initial radii,
% when other polytopes have small initial radii.

max_r_dist = old_max_radii; % largest possible radius for each polytope
min_r_dist = min_rad; % smallest possible radius
Ntruncations = sum((new_r_dist>max_r_dist)+(new_r_dist<min_r_dist));

% Warn the user:
if flag_do_debug
    if  Ntruncations> 0
        st = dbstack;
        warning('In function: %s, in file: %s\n\tThe standard deviation will skewed due to truncated distribution.\n\tThe number of truncations per total population of polytopes is: %.0d out of %.0d',...
            st(1).name,st(1).file,Ntruncations,Nradii);
    end
end

% Force the tails of the distribution that stick out of the truncation
% limits, to the truncation limits. So we can re-check how much this
% truncation affects the mean. Specifically: Truncating the edges will
% cause mean to shift. So move the mean around slightly, but each motion
% may require truncation again. Essentially, this piles all the tails of
% the distribution back onto the edges. (not very normal, btw). NOTE: if a
% single-valued max value is used for all polytopes (not just the max for
% EACH polytope), then the code below only works if max radius on all
% polytopes is about the same. If there are a few polytopes with large
% radii, while others are small, then this calculation will still give
% wrong answers.

new_r_dist(new_r_dist>max_r_dist) = max_r_dist(new_r_dist>max_r_dist); % truncate any values that are too large
new_r_dist(new_r_dist<min_r_dist) = min_r_dist; % truncate any values that are too small
while abs(mean(new_r_dist) - des_radius) > 1e-10
    new_r_dist = new_r_dist + (des_radius-mean(new_r_dist));
    new_r_dist(new_r_dist>max_r_dist) = max_r_dist(new_r_dist>max_r_dist);
    new_r_dist(new_r_dist<min_r_dist) = min_r_dist;
end

if flag_do_debug
    mu_final = mean(new_r_dist);
    sigma_final = std(new_r_dist);
    fprintf(1,'Target distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',mu_final);
    fprintf(1,'\tStd dev: %.4f\n',sigma_final);

   figure(fig_for_debug+2);  
   hold on;
   histogram(new_r_dist,20)
   title(sprintf('Histogram of bounded target radii. Mean: %.4f, Std-dev: %.4f. Targets are: %.4f and %.4f',...
       mu_final,sigma_final,...
       des_radius,...
       sigma_radius));

end

%% shrink polytopes to achieve the distribution
% Sort the radii changes by size. Effectively, this randomizes the changes
% against the old polytopes, since the size of the new radii were
% determined randomly.
% The following lines effectively forces the "sorted" radii and indices to
% match the new_r_dist.
%
% Sorting was used to ensure the biggest old polytope
% becomes the biggest new polytope. This is an easy way to ensure
% there is usually a large enough polytope to shrink for each new polytope.
% Othewise, we might have to grow some obstacles to achieve the
% distribution and this could cause overlapping since the original
% polytopes are tiled together.

[new_radii_sorted,ob_index] = sort(new_r_dist);
% % Uncomment to skip sorting:
% new_radii_sorted = new_r_dist;
% ob_index = find(new_r_dist>=0);

% Check that the old polytopes are large enough to shrink and 
% achieve the new radius distribution. Want all the changes to be smaller
% than -2 times the minimum radius, to ensure we do not get singular
% polytopes.
change_in_radii = sort(old_max_radii)'-sort(new_r_dist);
Num_goal_polys_smaller_than_start = sum(change_in_radii>=-2*min_rad);
if  Num_goal_polys_smaller_than_start < Nradii
    error('distribution is unachievable with generated map')
end

% Initialize the shrunk polytopes structure array, and tolerance for
% distance between vertices, below which vertices are merged into one.
shrunk_polytopes = polytopes;
tolerance = 1e-5; % Units are (implied) kilometers

% Loop through each polytope, shrinking it to the reference size
for ith_radii = 1:length(new_radii_sorted)
    shrinker = polytopes(ob_index(ith_radii)); % obstacle to be shrunk
    des_rad = new_radii_sorted(ith_radii);
        
    % assign to shrunk_polytopes
    shrunk_polytopes(ob_index(ith_radii)) = ...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad,tolerance);
end

% Fill in mu and sigma values from final result
final_max_radii = [shrunk_polytopes.max_radius]';
mu_final = mean(final_max_radii);
sigma_final = std(final_max_radii);

if flag_do_debug
    fprintf(1,'Final distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',mu_final);
    fprintf(1,'\tStd dev: %.4f\n',sigma_final);
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

if flag_do_plot
    figure(fig_num);
    hold on
    
    % Plot the input polytopes in red
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'r',2);
    
    % plot the shrunk in blue
    fcn_MapGen_plotPolytopes(shrunk_polytopes,fig_num,'b',2);

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


