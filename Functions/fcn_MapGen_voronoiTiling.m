function [polytopes, eachGeneratorPolytopes] = fcn_MapGen_voronoiTiling(seedGeneratorNames, seedGeneratorRanges, varargin)
% fcn_MapGen_voronoiTiling generates a map fully covered by polytope
% obstacles that are tiled via Voronoi tilings. The tilings are specified
% by one or more sets of seed points. Each seed point set can be specified
% independently in type, settings, and AABB range for each set. A mapStretch
% can be added to the final map.
%
%     FORMAT:
%
%     [polytopes] = fcn_MapGen_voronoiTiling(...
%       seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
%       seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
%       (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
%       (mapStretch),...        % vector or cellArrayOf_vectors to specify how to stretch X and Y axis
%       (fig_num));
%
%     INPUTS:
%
%     seedGeneratorNames: a single string or cell array of strings
%     specifying the names of the seed generator type. Allowable names
%     include:
%       'haltonset'    - generates points using sampling from MATLAB's
%                        built-in haltonset
%       'sobolset'     - generates points using sampling from MATLAB's
%                        built-in sobelset
%       'lhsdesign'    - generates points using call to lhsdesign, latin
%                        hyper square sampling
%       'rand'         - generates points using call to rand
%       'randn'        - generates points using call to randn (NOTE:
%                        resamples points that are outside user-given AABB
%                        to force all points to be within AABB, so points
%                        are technically not normal)
%       'latin'        - same as lhsdesign
%      See more about each at:
%      https://www.mathworks.com/help/stats/generating-quasi-random-numbers.html
%
%     seedGeneratorRanges: a vector or cell array of vectors
%     specifying the point ranges or settings for the seed generator type.
%     Usually, the settings are the [low high] values for each of the
%     methods above. The number of points is high-low. See the examples.
%     Defaults are listed above.
%
%     (OPTIONAL INPUTS)
%
%     AABB: a vector or cell array of vectors specifying the axis-aligned
%     bounding box, defining the low and high range each seed pointr
%     generation. Each AABB is defined as: [xlow ylow xhigh yhigh].
%     Defaults are [0 0 1 1]
%
%     mapStretch: a vector or cell array of vectors specifying the [x,y]
%     scaling factor to allow the values from the combined set to be scaled
%     to fit maps shapes other than square. Default is a unit mapStretch:
%     [1 1]
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose. Default is not to plot.
%
%     OUTPUTS:
%
%     polytopes: a polytope structure. See fcn_MapGen_polytopeFillEmptyPoly
%
%     eachGeneratorPolytopes: a cell array of polytopes for each generator
%
%     DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     fcn_MapGen_polytopeCentroidAndArea
%     fcn_MapGen_plotPolytopes
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_voronoiTiling
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
% 2025_07_12 by S. Brennan
% -- renamed function from fcn_MapGen_mixedSetVoronoiTiling
% -- changed to cell array input types, rather than structures (easier to
%    use)

% TO DO:
% -- (add here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; %       Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %       Flag to plot the results for debugging
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
        if nargin < 2 || nargin > 5
            error('Incorrect number of input arguments')
        end

        % Check the seedGeneratorNames input
        assert(iscell(seedGeneratorNames) || isstring(seedGeneratorNames) || ischar(seedGeneratorNames));

        % Check the seedGeneratorNames input
        assert(iscell(seedGeneratorRanges) || isnumeric(seedGeneratorRanges));
        if iscell(seedGeneratorNames)
            assert(iscell(seedGeneratorRanges))
            assert(length(seedGeneratorRanges) == length(seedGeneratorNames))
        end
        if isstring(seedGeneratorNames) || ischar(seedGeneratorNames)
            assert(isnumeric(seedGeneratorRanges))
        end

    end
end

% check variable argument AABBs
AABBs = []; % use default AABBs
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        AABBs = temp;
        if flag_check_inputs
            % Check the AABBs input
            assert(iscell(AABBs) || isnumeric(AABBs));
            if iscell(seedGeneratorNames)
                assert(iscell(AABBs))
                assert(length(AABBs) == length(seedGeneratorNames))
            end
            if isstring(seedGeneratorNames) || ischar(seedGeneratorNames)
                assert(isnumeric(AABBs))
            end
        end
    end
end


% check variable argument mapStretch
mapStretch = []; % use default mapStretch
if 4 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        mapStretch = temp;
        if flag_check_inputs
            % Check the mapStretch input
            assert(iscell(mapStretch) || isnumeric(mapStretch));
            if iscell(seedGeneratorNames)
                assert(iscell(mapStretch))
                assert(length(mapStretch) == length(seedGeneratorNames))
            end
            if isstring(seedGeneratorNames) || ischar(seedGeneratorNames)
                assert(isnumeric(mapStretch))
            end
        end
    end
end

% STOPPED INPUT EDITS HERE

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

totalMapStretch = [1 1];

% Convert everything to cell arrays
if ~iscell(seedGeneratorNames)
    seedGeneratorNames = {seedGeneratorNames};
    seedGeneratorRanges = {seedGeneratorRanges};
    AABBs = {AABBs};
    mapStretch = {mapStretch};
end

Ngenerators = length(seedGeneratorNames);

% Fill in starting values
eachGeneratorSeedPoints = cell(Ngenerators,1);
eachGeneratorPolytopes = cell(Ngenerators,1);
all_points = [];
min_xy_offset = [ inf  inf];
max_xy_offset = [-inf -inf];

% Fill in set points
for ith_generator = 1:Ngenerators
    ith_seedPointsRaw = fcn_INTERNAL_generateSetPoints(seedGeneratorNames{ith_generator}, seedGeneratorRanges{ith_generator}, AABBs{ith_generator});
    eachGeneratorSeedPoints{ith_generator} = ith_seedPointsRaw.*mapStretch{ith_generator};
    all_points = [all_points; eachGeneratorSeedPoints{ith_generator}];  %#ok<AGROW>
    min_xy_offset = min(min_xy_offset,AABBs{ith_generator}(1:2).*[mapStretch{ith_generator}(1) mapStretch{ith_generator}(2)]);
    max_xy_offset = max(max_xy_offset,AABBs{ith_generator}(3:4).*[mapStretch{ith_generator}(1) mapStretch{ith_generator}(2)]);
end

% Calculate the axis-aligned bounding box for the sets
AABB_final = [min_xy_offset max_xy_offset];

% Calculate the Voronoi polytopes, % fill polytopes from tiling
[V,C] = voronoin(all_points);
polytopes = fcn_MapGen_generatePolysFromTiling(all_points,V,C, AABB_final, totalMapStretch,[], -1);

% Repeat for each generator set
if Ngenerators>1
    for ith_generator = 1:Ngenerators
        [V,C] = voronoin(eachGeneratorSeedPoints{ith_generator});
        eachGeneratorPolytopes{ith_generator,1} = fcn_MapGen_generatePolysFromTiling(...
            eachGeneratorSeedPoints{ith_generator},V,C, AABB_final, totalMapStretch, [], -1);
    end
else
    % the result will be the same
    eachGeneratorPolytopes{1,1} = polytopes;
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

    % plot the polytopes
    % FIX THIS FOR BETTER CLARITY
    fcn_MapGen_plotPolytopes(polytopes,fig_num,'b',2);

    % plot the seed points for each Generator
    colorEachGenerator = zeros(Ngenerators,3);
    for ith_generator = 1:Ngenerators
        % Plot the seed points
        moved_seed_points = eachGeneratorSeedPoints{ith_generator}.*totalMapStretch;
        h_plot = plot(moved_seed_points(:,1),moved_seed_points(:,2),'.','Markersize',10,'DisplayName',sprintf('Seed points %.0f',ith_generator));
        colorEachGenerator(ith_generator,:) = get(h_plot,'Color');

        % % Plot the polys for each genrator
        % (COMMENTED OUT BECAUSE IT LOOKS MESSY)
        % fcn_MapGen_plotPolytopes(eachGeneratorPolytopes{ith_generator,1},fig_num,'b',1,colorEachGenerator(ith_generator,:));

        % % plot the means
        % (COMMENTED OUT BECAUSE THE MEANS INCLUDE MESSY POLYTOPES THAT ARE
        % NOT INSIDE CORRECT AREAS)
        % NpolysThisGenerator = length(eachGeneratorPolytopes{ith_generator,1});
        % temp = zeros(NpolysThisGenerator,2);
        % for ith_poly = 1:NpolysThisGenerator
        %     temp(ith_poly,:) = eachGeneratorPolytopes{ith_generator,1}(ith_poly).mean;
        % end
        % plot(temp(:,1),temp(:,2),'o','Markersize',3,'Color',colorEachGenerator(ith_generator,:),'DisplayName',sprintf('Means %.0f',ith_generator))

    end

    legend('Interpreter','none');

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

%% fcn_INTERNAL_generateSetPoints
function seedPoints = fcn_INTERNAL_generateSetPoints(name, set_range, AABB)

lower_value = AABB(1:2);
scale =  AABB(3:4) -  AABB(1:2);

switch lower(name)
    case{'haltonset'}  % pull halton set
        halton_points = haltonset(2);
        points_scrambled = scramble(halton_points,'RR2'); % scramble values

        % pick values from halton set
        low_pt = set_range(1,1);
        high_pt = set_range(1,2);

        % Produce seed points
        seedPoints = points_scrambled(low_pt:high_pt,:).*scale + lower_value;

    case{'sobolset'} % pull Sobel set
        Sobol_points = sobolset(2);
        points_scrambled = scramble(Sobol_points,'MatousekAffineOwen'); % scramble values

        % pick values from Sobol set
        low_pt = set_range(1,1);
        high_pt = set_range(1,2);

        % Produce seed points
        seedPoints = points_scrambled(low_pt:high_pt,:).*scale + lower_value;

    case{'lhsdesign', 'latin'}
        % pull latin set
        Npoints = max(set_range) - min(set_range) + 1;
        latin_points = lhsdesign(Npoints,2);

        % pick values from latin set
        seedPoints = latin_points.*scale + lower_value;

    case{'rand'}
        % pull rand set
        Npoints = max(set_range) - min(set_range) + 1;
        rand_points = rand(Npoints,2);

        % pick values from rand set
        seedPoints = rand_points.*scale + lower_value;

    case{'randn'}  % pull rand set
        Npoints = max(set_range) - min(set_range) + 1;
        rand_points = 0.15*randn(Npoints*2,1)+0.5;

        % Prevent random normal points from being outside the range. If they are, we
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
        seedPoints = rand_points.*scale + lower_value;

    otherwise
        warning('on','backtrace');
        warning('fcn_MapGen_voronoiTiling was given an unrecognized geenrator type: %s. An error will be thrown.', lower(name));
        error('Unknown method given');
end

end % Ends fcn_INTERNAL_generateSetPoints
