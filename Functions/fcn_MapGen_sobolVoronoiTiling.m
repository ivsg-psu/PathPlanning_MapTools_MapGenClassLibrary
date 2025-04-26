function [polytopes] = ...
    fcn_MapGen_sobolVoronoiTiling(Sobol_range,varargin)
% fcn_MapGen_sobolVoronoiTiling generates a map with
% obstacles perfectly tiled together using the Voronoi cells generated from
% the Sobol sequence. See more about this at:
% https://www.mathworks.com/help/stats/generating-quasi-random-numbers.html
%
% FORMAT:
% 
% [polytopes] = ...
%    fcn_MapGen_sobolVoronoiTiling(Sobol_range,(stretch),(fig_num))
%
% INPUTS:
%
%    Sobol_range: 1 x 2 vector of integers specifying the [low high] range
%    of Sobol point indices to use to generate the tiling, where  
%    low: the lowest point index to use in the Sobol sequence
%    high: the highest point index to use in the Sobol sequence
%     
%    (OPTIONAL INPUTS)
%
%    stretch: [x,y] scaling factors to allow the values from the Sobol set
%    to be scaled to fit maps shapes other than square
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
% POLYTOPES: a 1-by-n seven field structure, where n <= number of polytopes
%   with fields:
%   vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%   xv: a 1-by-m vector of vertice x-coordinates
%   yv: a 1-by-m vector of vertice y-coordinates
%   distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%   mean: centroid xy coordinate of the polytope
%   area: area of the polytope
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_MapGen_polytopeCentroidAndArea
%      fcn_MapGen_plotPolytopes
%
% EXAMPLES:
%      
% See the script: script_test_fcn_MapGen_sobolVoronoiTiling
% for a full test suite.
%
% This function was written on 2019_06_13 by Seth Tau
% Questions or comments? sbrennan@psu.edu 

% REVISION HISTORY:
% 2021_06_06 
% -- edited by S. Brennan to put it into MapGen format
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions

% TO DO:
% -- check cross product around entire polytope
% -- add unit normal vectors for each edge
% -- buffer the Sobol set to be sure polytopes are well formed on edges
% -- force correct number of polytopes?

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(1,3);

        % Check the Sobol_range input
        fcn_DebugTools_checkInputsToFunctions(...
            Sobol_range, '2column_of_integers');

    end
end

% check variable argument
stretch = [1 1]; % default stretch value
if 2 == nargin
    stretch = varargin{1};
    
    % Check the stretch input
    fcn_DebugTools_checkInputsToFunctions(...
        stretch, '2column_of_numbers',1);
       
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (3 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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
%% pull Sobol set
Sobol_points = sobolset(2);
points_scrambled = scramble(Sobol_points,'MatousekAffineOwen'); % scramble values

%% pick values from Sobol set
low_pt = Sobol_range(1,1);
high_pt = Sobol_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
% V = V.*stretch;

%% fill polytopes from tiling
AABB = [0 0 1 1];
polytopes = fcn_MapGen_generatePolysFromTiling(seed_points,V,C,AABB, stretch);


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
 
    % plot the polytopes
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2,[0 1 0 1]);

    % plot the seed points
    plot(seed_points(:,1),seed_points(:,2),'r.','Markersize',10);

    
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function


