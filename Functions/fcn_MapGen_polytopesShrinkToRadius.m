function [shrunkPolytopes,muFinal,sigmaFinal] = ...
    fcn_MapGen_polytopesShrinkToRadius(...
    polytopes,...
    desiredRadius,...
    sigmaRadius,...
    minRadius,...
    varargin)
% fcn_MapGen_polytopesShrinkToRadius shrinks the polytopes to achieve the
% desired mean radius and specified variance
%
% FORMAT:
% 
% [shrunkPolytopes,muFinal,sigmaFinal] = ...
%     fcn_MapGen_polytopesShrinkToRadius(...
%     polytopes,...
%     desiredRadius,...
%     sigmaRadius,...
%     minRadius,...
%     (fig_num))
%
% INPUTS:
%
%     polytopes: original polytopes with same fields as shrunkPolytopes
%
%     desiredRadius: desired average max radius   
%
%     sigmaRadius: desired variance in the radii 
%
%     minRadius: minimum acceptable radius
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
%     shrunkPolytopes: a structure of shrunken polytopes, 
%
%     muFinal: final average maximum radius achieved
%
%     sigmaFinal: final variance achieved
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_polytopeShrinkToRadius
%     fcn_MapGen_plotPolytopes
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
% 2025_07_14 by Sean Brennan
% -- cleaned up variable naming for clarity and avoid underscores
% 2025_07_17 by Sean Brennan
% -- standardized Debugging and Input checks area, Inputs area
% -- made codes use MAX_NARGIN definition at top of code, narginchk
% -- made plotting flag_do_plots and code consistent across all functions

% TO DO
% -- Vectorize the for loop if possible
% -- check inputs are positive numbers where appropriate (e.g. make a
%    % "positive number" check

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 5; % The largest Number of argument inputs to the function
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
    debug_fig_num = 999978; %#ok<NASGU>
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
        narginchk(4,MAX_NARGIN);

        % Check the polytopes input
        fcn_DebugTools_checkInputsToFunctions(polytopes, 'polytopes');

        % Check the desiredRadius input
        fcn_DebugTools_checkInputsToFunctions(desiredRadius, 'positive_1column_of_numbers',1);

        % Check the sigmaRadius input
        fcn_DebugTools_checkInputsToFunctions(sigmaRadius, '1column_of_numbers',1);

        % Check the minRadius input
        fcn_DebugTools_checkInputsToFunctions(minRadius, 'positive_1column_of_numbers',1);

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


%% find current distribution
old_max_radii = [polytopes.max_radius]';
Nradii = length(old_max_radii);
old_r_mu = mean(old_max_radii);
old_r_sigma = std(old_max_radii);

if flag_do_debug
    fprintf(1,'Target distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',desiredRadius);
    fprintf(1,'\tStd dev: %.4f\n',sigmaRadius);
    
    fprintf(1,'Input distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',old_r_mu);
    fprintf(1,'\tStd dev: %.4f\n',old_r_sigma);
    
    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_for_debug,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_for_debug)); %#ok<NASGU>
    
    figure(fig_for_debug+1);
    histogram(old_max_radii,20)
    title(sprintf('Histogram of input radii. Mean: %.4f, Std-dev: %.4f. Targets are: %.4f and %.4f',...
        old_r_mu,old_r_sigma,...
        desiredRadius,...
        sigmaRadius));
end


if old_r_mu < desiredRadius
    error('cannot achieve the desired radius by shrinking because average radius is already smaller than desired radius')
end

%% determine desired distribution
new_r_dist = normrnd(desiredRadius,sigmaRadius,[Nradii,1]);

% adjust to ensure the mean value is mu. SETH: is this necessary?
new_r_dist = new_r_dist + (desiredRadius-mean(new_r_dist)); 


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
        desiredRadius,...
        sigmaRadius));
   
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
min_r_dist = minRadius; % smallest possible radius
Ntruncations = sum((new_r_dist>max_r_dist)+(new_r_dist<min_r_dist));

% Warn the user:
if flag_do_debug
    if  Ntruncations> 0
        st = dbstack;
        warning('on','backtrace');
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
while abs(mean(new_r_dist) - desiredRadius) > 1e-10
    new_r_dist = new_r_dist + (desiredRadius-mean(new_r_dist));
    new_r_dist(new_r_dist>max_r_dist) = max_r_dist(new_r_dist>max_r_dist);
    new_r_dist(new_r_dist<min_r_dist) = min_r_dist;
end

if flag_do_debug
    muFinal = mean(new_r_dist);
    sigmaFinal = std(new_r_dist);
    fprintf(1,'Target distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',muFinal);
    fprintf(1,'\tStd dev: %.4f\n',sigmaFinal);

   figure(fig_for_debug+2);  
   hold on;
   histogram(new_r_dist,20)
   title(sprintf('Histogram of bounded target radii. Mean: %.4f, Std-dev: %.4f. Targets are: %.4f and %.4f',...
       muFinal,sigmaFinal,...
       desiredRadius,...
       sigmaRadius));

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
Num_goal_polys_smaller_than_start = sum(change_in_radii>=-2*minRadius);
if  Num_goal_polys_smaller_than_start < Nradii
    error('distribution is unachievable with generated map')
end

% Initialize the shrunk polytopes structure array, and tolerance for
% distance between vertices, below which vertices are merged into one.
shrunkPolytopes = polytopes;

% Loop through each polytope, shrinking it to the reference size
for ith_radii = 1:length(new_radii_sorted)
    shrinker = polytopes(ob_index(ith_radii)); % obstacle to be shrunk
    des_rad = new_radii_sorted(ith_radii);
        
    % assign to shrunkPolytopes
    shrunkPolytopes(ob_index(ith_radii)) = ...
        fcn_MapGen_polytopeShrinkToRadius(...
        shrinker,des_rad, -1);
end

% Fill in mu and sigma values from final result
final_max_radii = [shrunkPolytopes.max_radius]';
muFinal = mean(final_max_radii);
sigmaFinal = std(final_max_radii);

if flag_do_debug
    fprintf(1,'Final distrubution statistics:\n');
    fprintf(1,'\tMean: %.4f\n',muFinal);
    fprintf(1,'\tStd dev: %.4f\n',sigmaFinal);
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
    figure(fig_num);
    hold on
    
    % Plot the input polytopes in red
    % fcn_MapGen_OLD_plotPolytopes(polytopes,fig_num,'r',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [1 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytopes, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>
    
    % plot the shrunk in blue
    % fcn_MapGen_OLD_plotPolytopes(shrunkPolytopes,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(shrunkPolytopes, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

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


