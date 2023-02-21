function [shrunk_polytopes,mu_final,sigma_final] = ...
    fcn_MapGen_polytopesShrinkFromEdges(...
    polytopes,...
    des_gap_size,...
    varargin)
% fcn_MapGen_polytopesShrinkFromEdges shrinks the polytopes to achieve the
% desired gap size between polytopes
%
% FORMAT:
%
% [shrunk_polytopes] = ...
%     fcn_MapGen_polytopesShrinkFromEdges(...
%     polytopes,...
%     des_gap_size,...
%     (fig_num))
%
% INPUTS:
%
%     POLYTOPES: original polytopes with same fields as shrunk_polytopes
%
%     DES_GAP_SIZE: desired normal gap size between polytopes, note this is
%     relative, not absolute e.g. if the input polytopes have no gap and a
%     des_gap_size of 0.005 km is input, the resulting gap size will be
%     0.005 km.  However if the input polytopes already have a gap size of
%     0.002 km the resulting gap size will be 0.007 km for a des_gap_size
%     of 0.005 km.
%
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a
%     figure to be drawn showing results.
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
% DEPENDENCIES:
%
%     fcn_MapGen_checkInputsToFunctions
%     fcn_MapGen_plotPolytopes
%     fcn_MapGen_polytopeShrinkFromEdges
%
% EXAMPLES:
%
% For additional examples, see: script_test_fcn_MapGen_polytopesShrinkFromEdges
%
% This function was written on 2022_01_17 by Steve Harentt using fcn_MapGen_polytopesShrinkToRadius
% as a template
% Questions or comments? sjh6473@psu.edu
%

% Revision History:

% TO DO
% copied from TODOs in fcn_MapGen_polytopesShrinkToRadius
% -- Vectorize the for loop if possible
% -- check inputs are positive numbers where appropriate (e.g. make a
% "positive number" check
% -- add non uniform shrinking and allow specification of a gap size variance,
%    similar to how a radius variance is allowed in fcn_MapGen_polytopesShrinkToRadius.m

%% Debugging and Input checks
% set an environment variable on your machine with the getenv function...
% in the Matlab console.  Char array of '1' will be true and '0' will be false.
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
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');

    % Check the des_radius input
    fcn_MapGen_checkInputsToFunctions(...
        des_gap_size, 'positive_1column_of_numbers',1);
end


% Does user want to show the plots?
if  3 == nargin
    fig_num = varargin{end};
    figure(fig_num);
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

initial_stats = fcn_MapGen_polytopesStatistics(polytopes);
initial_average_max_rad = initial_stats.average_max_radius;

if flag_do_debug
    fprintf(1,'Target statistics:\n');
    fprintf(1,'\tGap size: %.4f\n',des_gap_size);

    fprintf(1,'Input field statistics:\n');
    fprintf(1,'\tAvg Max Rad: %.4f\n',initial_average_max_rad);

    fcn_MapGen_plotPolytopes(polytopes,fig_for_debug,'b',2);
end

%% shrink polytopes to achieve the distribution
% Check that the old polytopes are large enough to shrink and
% achieve the new radius distribution. Want all the changes to be smaller
% than -2 times the minimum radius, to ensure we do not get singular
% polytopes.

% Initialize the shrunk polytopes structure array, and tolerance for
% distance between vertices, below which vertices are merged into one.
clear shrunk_polytopes
shrinker = polytopes(1); % obstacle to be shrunk
temp = ...
    fcn_MapGen_polytopeShrinkFromEdges(...
    shrinker,des_gap_size/2);
shrunk_polytopes(length(polytopes)) = temp;

% Loop through each polytope, shrinking it to the reference size
for ith_poly = 1:length(polytopes)
    shrinker = polytopes(ith_poly); % obstacle to be shrunk

    % assign to shrunk_polytopes
    % gap_size over 2 is the normal distance to pull edges in
    shrunk_polytopes(ith_poly) = ...
        fcn_MapGen_polytopeShrinkFromEdges(...
        shrinker,des_gap_size/2);
end


final_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes);
final_average_max_rad = final_stats.average_max_radius;
if flag_do_debug
    fprintf(1,'Final distrubution statistics:\n');
    fprintf(1,'\tAvg mag rad: %.4f\n',final_average_max_rad);
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
