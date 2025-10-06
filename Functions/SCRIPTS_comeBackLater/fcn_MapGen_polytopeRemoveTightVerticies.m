function [ ...
    cleanedPolytope ...
    ] = ...
    fcn_MapGen_polytopeRemoveTightVerticies( ...
    polytope, ...
    tolerance, ...
    varargin...
    )
% fcn_MapGen_polytopeRemoveTightVerticies
% removes verticies of polytopes that are too close to each other,
% measured by a tolerance
%
% Sometimes, when shrinking, the new verticies are particularly close to
% each other to where an edge has a trivial length. To prevent this, we
% get rid of one of any vertices that are too close to each other. This
% proximity is set by a user-defined tolerance.
%
% FORMAT:
%
%    [ ...
%    cleanedPolytope ...
%    ] = ...
%    fcn_MapGen_polytopeRemoveTightVerticies( ...
%    polytope, ...
%    tolerance, ...
%    (fig_num) ...
%    )
%
% INPUTS:
%
%     polytope: an individual structure or structure array of 'polytopes'
%     type that stores the polytopes to be evaluated
%
%     tolerance: a numeric value that defines how close points should be
%     to be removed
%
%     (optional inputs)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
%
% OUTPUTS:
%
%     cleanedPolytope: the resulting polytope after close edges are
%     removed.
%
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
%     fcn_MapGen_fillPolytopeFieldsFromVerticies
%
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_polytopeRemoveTightVerticies
% for a full test suite.
%
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of function
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed call to fcn_MapGen_fillPolytopeFieldsFromVertices
% 2025_10_05 by Sean Brennan
% -- removed call to fcn_MapGen_fillPolytopeFieldsFromVertices 
%    % replaced with fcn_MapGen_polytopesFillFieldsFromVertices


% TO DO
% -- none

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
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the polytopes input, make sure it is 'polytopes' type
        fcn_DebugTools_checkInputsToFunctions(...
            polytope, 'polytopes');

        % Check the tolerance input, make sure it is '1column_of_numbers' type,
        % e.g. 1x1 numerical data.
        fcn_DebugTools_checkInputsToFunctions(...
            tolerance, '1column_of_numbers',[1 1]);

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





% pull values
vertices = polytope.vertices;
centroid = polytope.mean;


% Work with tolerance squared, since it avoids square-root calculations
tol_squared = tolerance.^2;

% find all vertices that are larger away from next one, relative to
% tolerance
good_ind = ...
    sum((vertices(1:end-1,:)-vertices(2:end,:)).^2,2)>(tol_squared);

% Check how many are good. If only 2 points are left, then the
% result is just a line. If 1 or zero, then result is just a point.
if sum(good_ind)>2 % sufficient good points to make a shape
    new_vert = vertices(good_ind,:);
elseif sum(good_ind)==2 % line shape
    new_vert = vertices(good_ind,:);

    % The line may not go through the centroid, which is odd. We force
    % this by removing the point closest to the centroid
    distances_to_centroid = sum((new_vert-centroid).^2,2).^0.5;
    if distances_to_centroid(1,1)>distances_to_centroid(2,1)
        new_vert(2,:)=centroid;
    else
        new_vert(1,:)=centroid;
    end

    new_vert = [new_vert; flipud(new_vert)];
else % singular shape (i.e. point) or no shape
    new_vert = [centroid; centroid; centroid];
end

% adjust polytopes
cleanedPolytope.vertices = [new_vert; new_vert(1,:)];
cleanedPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(cleanedPolytope);

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

    % Plot the cetroid in black
    plot(centroid(:,1),centroid(:,2),'ko','Markersize',10);

    % Plot the input polytope in red
    % fcn_MapGen_OLD_plotPolytopes(polytope,fig_num,'r',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [1 0 0];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(polytope, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

    % plot the output polytope in blue
    % fcn_MapGen_OLD_plotPolytopes(cleanedPolytope,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(cleanedPolytope, (plotFormat), (fillFormat), (fig_num)); %#ok<NASGU>

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

