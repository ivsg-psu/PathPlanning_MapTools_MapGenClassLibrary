function [polytopes] = ...
    fcn_MapGen_mixedSetVoronoiTiling(mixedSet, varargin)

warning('on','backtrace');
warning('fcn_MapGen_mixedSetVoronoiTiling is being deprecated. Use fcn_MapGen_voronoiTiling instead.');

% fcn_MapGen_mixedSetVoronoiTiling generates a map with
% obstacles perfectly tiled together using the Voronoi cells generated from
% the Halton sequence
%
% FORMAT:
%
% [polytopes] = ...
%    fcn_MapGen_mixedSetVoronoiTiling(Halton_range,(stretch),(fig_num))
%
% INPUTS:
%
%    mixedSet: a structure with the following format
%
%             mixedSet(1).name = 'haltonset';
%             mixedSet(1).settings = set_range;
%             mixedSet(1).AABB = [0 0 1 1];
%
%         where:
%
%             name: the name of the type of set to use, which can be:
%                  'haltonset','sobolset','lhsdesign','rand','randn'
%
%             settings: the [low high] values for each of the methods
%                   above. The number of points is high-low.
%
%             AABB: the axis-aligned bounding box, defining the low and
%                   high range for the number generation. The box is
%                   defined as: [xlow ylow xhigh yhigh]
%
%
%    (OPTIONAL INPUTS)
%
%    stretch: [x,y] scaling factors to allow the values from the combined
%    set to be scaled to fit maps shapes other than square
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
% See the script: script_test_fcn_MapGen_mixedSetVoronoiTiling
% for a full test suite.
%
% This function was written on 2021_07_08 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 2021_07_08
% -- Created by S. Brennan
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions


% TO DO:
% -- (add here)


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
        if nargin < 1 || nargin > 3
            error('Incorrect number of input arguments')
        end

        % Check the mixedSet input
        fcn_DebugTools_checkInputsToFunctions(...
            mixedSet, 'mixedSet');

    end
end

% check variable argument
stretch = [1 1]; % default stretch value
if 2 <= nargin
    stretch = varargin{1};

    if (0==flag_max_speed)
        % Check the stretch input
        fcn_DebugTools_checkInputsToFunctions(...
            stretch, '2column_of_numbers',1);
    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  3 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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


% Fill in starting values
set_points{length(mixedSet)} = [];
all_points = [];
min_xy_offset = [ inf  inf];
max_xy_offset = [-inf -inf];

% Fill in set points
for ith_set = 1:length(mixedSet)
    set_points{ith_set} = INTERNAL_fcn_generateSetPoints(mixedSet(ith_set));
    all_points = [all_points; set_points{ith_set}];  %#ok<AGROW>
    min_xy_offset = min(min_xy_offset,mixedSet(ith_set).AABB(1:2));
    max_xy_offset = max(max_xy_offset,mixedSet(ith_set).AABB(3:4));
end

% Calculate the axis-aligned bounding box for the sets
AABB = [min_xy_offset max_xy_offset];

% Calculate the Voronoi polytopes
[V,C] = voronoin(all_points);

% fill polytopes from tiling
polytopes = fcn_MapGen_generatePolysFromTiling(all_points,V,C, AABB, stretch);


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
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2);

    % plot the seed points
    for i_seed = 1:length(set_points)
        moved_seed_points = set_points{i_seed}.*stretch;
        plot(moved_seed_points(:,1),moved_seed_points(:,2),'.','Markersize',10);
    end

    % plot the means in black
    temp = zeros(length(polytopes),2);
    for ith_poly = 1:length(polytopes)
        temp(ith_poly,:) = polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);

end

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


function seed_points = INTERNAL_fcn_generateSetPoints(set)
% 'haltonset','sobolset','lhsdesign','rand','randn'
set_range = set.settings;
lower_value = set.AABB(1:2);
scale =  set.AABB(3:4) -  set.AABB(1:2);

switch lower(set.name)
    case{'haltonset'}  % pull halton set
        halton_points = haltonset(2);
        points_scrambled = scramble(halton_points,'RR2'); % scramble values

        % pick values from halton set
        low_pt = set_range(1,1);
        high_pt = set_range(1,2);

        % Produce seed points
        seed_points = points_scrambled(low_pt:high_pt,:).*scale + lower_value;
    case{'sobolset'} % pull Sobel set
        Sobol_points = sobolset(2);
        points_scrambled = scramble(Sobol_points,'MatousekAffineOwen'); % scramble values

        % pick values from Sobol set
        low_pt = set_range(1,1);
        high_pt = set_range(1,2);
        seed_points = points_scrambled(low_pt:high_pt,:).*scale + lower_value;
    case{'lhsdesign'}
        % pull latin set
        Npoints = max(set_range) - min(set_range) + 1;
        latin_points = lhsdesign(Npoints,2);

        % pick values from latin set
        seed_points = latin_points.*scale + lower_value;
    case{'rand'}
        % pull rand set
        Npoints = max(set_range) - min(set_range) + 1;
        rand_points = rand(Npoints,2);

        % pick values from rand set
        seed_points = rand_points.*scale + lower_value;
    case{'randn'}  % pull rand set
        Npoints = max(set_range) - min(set_range) + 1;
        rand_points = 0.15*randn(Npoints*2,1)+0.5;

        % Prevent random points from being outside the range. If they are, we
        % resample just those points
        flag_points_are_good = 0;
        while flag_points_are_good==0
            bad_point_indices = find(rand_points>1 | rand_points<0);
            if ~isempty(bad_point_indices)
                rand_points(bad_point_indices) = 0.15*randn(length(bad_point_indices),1)+0.5;
            else
                flag_points_are_good = 1;
            end
        end
        rand_points = [rand_points(1:Npoints,1),rand_points(Npoints+1:end,1)];
        seed_points = rand_points.*scale + lower_value;

    otherwise
        error('Unknown method given');
end


end
