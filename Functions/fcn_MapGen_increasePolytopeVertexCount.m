function interpolatedPolytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes, resolution, varargin)
% fcn_MapGen_increasePolytopeVertexCount
% Given polytope field and a desired resolution distance, n, returns an equivalent
% polytope field with colinear vertices added to each polytope side such that
% there is a vertex every n units
% The utility of this is that if path planning is restricted to using polytope
% vertices as waypoints, this increases the number of options the planner has
% while keeping the obstacle field the same. This allows testing of whether
% the path plan is benefited by using waypoints that are within edges,
% rather than on just the vertices.
%
% FORMAT:
%     interpolatedPolytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes,resolution, (fig_num))
%
% INPUTS:
%
%     polytopes - the initial polytope field
%     resolution - the desired linear spacing between vertices along each polytope side
%
%     (optional inputs)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%
%     interpolatedPolytopes - a polytope field equivalent to the input but with vertices added
%     along the polytopes sides every RESOLUTION units such that each polytope now has more vertices
%
% DEPENDENCIES:
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_increasePolytopeVertexCount.m
% for a full test suite.
%
% Questions or comments? contact sjh6473@psu.edu or Sean Brennan,
% sbrennan@psu.edu

% REVISION HISTORY:
% 2021_10_13
% -- first written by Steve Harnett
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_07_15 by Sean Brennan
% -- improved docstrings to clarify fig_num input
% -- cleaned up variable names to avoid underscores.

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
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the polytopes input
        fcn_DebugTools_checkInputsToFunctions(polytopes, 'polytopes');

        % Check the resolution input, make sure it is [1 1]
        fcn_DebugTools_checkInputsToFunctions(resolution, '1column_of_numbers',[1 1]);

        
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

plotting_vertices = [];
plotting_updated_vertices = [];

% loop through all polytopes
for i = 1:length(polytopes)
    new_verts = [];

    if flag_do_plot
        plotting_vertices = [plotting_vertices; polytopes(i).vertices; nan nan]; %#ok<AGROW>
    end

    % loop through all vertices
    for j = 1:(size(polytopes(i).vertices,1)-1)
        % note the side range in x
        delta_x = polytopes(i).vertices(j,1) - polytopes(i).vertices(j+1,1);
        % note hte side range in y
        delta_y = polytopes(i).vertices(j,2) - polytopes(i).vertices(j+1,2);
        % find the side length
        side_length = sqrt(delta_x^2+delta_y^2);
        % determine number of vertices to place on this side per resolution
        num_verts_needed_this_side = side_length/resolution;
        % evenly space vertices in x since interp1 interpolates in x
        resolution_in_x = abs(delta_x)/num_verts_needed_this_side;
        % we expect the lower index vertex to have a lower x value than the higher
        % index vertex however...
        if polytopes(i).vertices(j,1) < polytopes(i).vertices(j+1,1)
            % sample x points for interp1
            xq = polytopes(i).vertices(j,1):resolution_in_x:polytopes(i).vertices(j+1,1);
            % see matlab documentation on interp1 for further examples of what this does
            vq = interp1([polytopes(i).vertices(j,1),polytopes(i).vertices(j+1,1)],...
                [polytopes(i).vertices(j,2),polytopes(i).vertices(j+1,2)],...
                xq);
            % ...if this side goes backwards (from larger x to smaller x) we
            % have to flip the interpolated points, additionally...
        elseif polytopes(i).vertices(j,1) > polytopes(i).vertices(j+1,1)
            xq = polytopes(i).vertices(j+1,1):resolution_in_x:polytopes(i).vertices(j,1);
            vq = interp1([polytopes(i).vertices(j+1,1),polytopes(i).vertices(j,1)],...
                [polytopes(i).vertices(j+1,2),polytopes(i).vertices(j,2)],...
                xq);
            % if this is a "backwards side" flip the interpolated vectors back
            vq = flip(vq);
            xq = flip(xq);
            % ...if this side is completely vertical (i.e. the vertices have the same x)...
        elseif polytopes(i).vertices(j,1) == polytopes(i).vertices(j+1,1)
            % ...then we can interpolate in y instead
            % we don't need to use interp1 because we know the slope of the side is inf.
            % so simple linearly space between the low y and high y
            vq = linspace(min(polytopes(i).vertices(j,2),polytopes(i).vertices(j+1,2)),...
                max(polytopes(i).vertices(j,2),polytopes(i).vertices(j+1,2)),...
                num_verts_needed_this_side);
            % and make an x vector of the appropriate size containing the only x value possible
            xq = ones(1,length(vq)).*polytopes(i).vertices(j,1);
        else
            error(sprintf('polytope in position %i could not be interpolated between vertex %i and %i. Does this vertex pair correctly represent a polytope side?',i,j,j+1)) %#ok<SPERR>
        end

        % log new vertices for this side
        new_verts_this_side = [xq; vq]';

        % append to array of new vertices for this polytope
        new_verts = [new_verts; new_verts_this_side]; %#ok<AGROW>

    end

    % update this polytope's fields based on the new vertices found
    polytopes(i).vertices = [new_verts;new_verts(1,:)];
    polytopes(i).xv = new_verts(:,1)';
    polytopes(i).yv = new_verts(:,2)';

    % plot the new vertices in a different color for comparison
    if flag_do_plot
        plotting_updated_vertices = [plotting_updated_vertices; polytopes(i).vertices; nan nan]; %#ok<AGROW>
    end
end
interpolatedPolytopes = polytopes;
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
    clf;
    hold on;

    % plot original polytope vertices
    plot(plotting_vertices(:,1),plotting_vertices(:,2),'k.-','DisplayName','original vertices','MarkerSize',10);
    %plot(plotting_vertices(:,1),plotting_vertices(:,2),'kx','DisplayName','original vertices');

    % Plot updated vertices
    plot(plotting_updated_vertices(:,1),plotting_updated_vertices(:,2),'r.','DisplayName','dense vertices','MarkerSize',5);

    legend;

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





