function [filled_polytopes] = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes, varargin)
% fcn_MapGen_fillPolytopeFieldsFromVertices
% Given a polytoope structure array where the vertices field is filled,
% calculates the values for all the other fields.
%
% FORMAT:
%
%    [filled_polytopes] = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes,
%    (is_nonconvex),...
%    (fig_num),...
%    )
%
% INPUTS:
%
%     polytopes: an individual structure or structure array of 'polytopes'
%     type that stores the polytopes to be filled
%
%     (optional inputs)
%
%     is_nonconvex - boolean flag indicating if there are or are not non-convex polytopes
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.

% OUTPUTS:
%
%     filled_polytopes: the polytopes array with all fields completed
%
%
% DEPENDENCIES:
%
%     fcn_MapGen_polytopeCentroidAndArea
%     fcn_DebugTools_checkInputsToFunctions
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
% 2023_03_13 by Sean Brennan
% -- added check and fix for ensuring verticies are counter-clockwise
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% -- fixed argument list to make figure number last to be consistent with
% other functions, and changed all functions/scripts that call this
% function to correct this
% 2025_07_15 by Sean Brennan
% -- fixed bug where parent_poly_id was not being filled

% TO DO
% 2025_07_09 - S. Brennan and K. Hayes
% --  need a tool to check if polytope is convex. This is causing some of
% the codes in Bounded_AStar to break

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


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 1 || nargin > 3
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it has vertices
    if ~isfield(polytopes,'vertices')
        error('Field of vertices was not found');
    end

    % Check the vertices input to have 4 or more rows, 2 columns
    %     fcn_DebugTools_checkInputsToFunctions(...
    %         polytopes.vertices, '2column_of_numbers',[4 5]);


end

% Does user specify if not convex?
is_nonconvex = 0;
if nargin >= 2
    temp = varargin{1};
    if ~isempty(temp) % Make sure the user is not giving empty input
        is_nonconvex = temp;
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

% Initialize variables
filled_polytopes = polytopes;
num_poly = length(polytopes);

% Loop over each polytope, filling in data for each
for ith_poly = 1:num_poly

    % check that verticies are counter-clockwise by calculating the angles.
    [angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(...
        filled_polytopes(ith_poly).vertices, -1);

    % Confirm that all angles are positive
    if ~all(angles>=0)
        if any(isnan(angles)) % This happens when there is a repeating point, which is a degenerate poly
            filled_polytopes(ith_poly).vertices = nan(length(filled_polytopes(ith_poly).vertices),2);
        elseif all(angles<=0)
            filled_polytopes(ith_poly).vertices = flipud(filled_polytopes(ith_poly).vertices);
        else
            if ~is_nonconvex
                fprintf(1,'Verticies:\n');
                for ith_vertex = 1:length(filled_polytopes(ith_poly).vertices)
                    fprintf(1,'%.2f %.2f\n',filled_polytopes(ith_poly).vertices(ith_vertex,1),filled_polytopes(ith_poly).vertices(ith_vertex,2))
                end
                fprintf(1,'\nAngles:\n');
                for ith_angle = 1:length(angles)
                    fprintf(1,'%.2f\n',angles(ith_angle));
                end
                warning('on','backtrace');
                warning('Non-convex polytope encountered. Code may break, forcing an error.');
                error('All vertices must be organized counter-clockwise, e.g. with positive cross-products. This error will get thrown for nonconvex obstacles.  Did you mean intentionally create nonconvex obstalces? If so you must pass the "is_nonconvex" flag to fcn_MapGen_fillPolytopeFieldsFromVertices to turn off this error.');
            end
        end
    end

    % adjust polytopes
    filled_polytopes(ith_poly).xv        = (polytopes(ith_poly).vertices(1:end-1,1)');
    filled_polytopes(ith_poly).yv        = (polytopes(ith_poly).vertices(1:end-1,2)');
    filled_polytopes(ith_poly).distances = ...
        sum((polytopes(ith_poly).vertices(1:end-1,:) - ...
        polytopes(ith_poly).vertices(2:end,:)).^2,2).^0.5;

    % Calculate the mean and area
    [filled_polytopes(ith_poly).mean,filled_polytopes(ith_poly).area] = ...
        fcn_MapGen_polytopeCentroidAndArea(polytopes(ith_poly).vertices, -1);

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

    % Fill in empty parent_poly_id
    filled_polytopes(ith_poly).parent_poly_id = [];
end


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
    % fcn_MapGen_OLD_plotPolytopes(filled_polytopes,fig_num,'b',2);
    plotFormat.LineWidth = 2;
    plotFormat.MarkerSize = 10;
    plotFormat.LineStyle = '-';
    plotFormat.Color = [0 0 1];
    fillFormat = [];
    h_plot = fcn_MapGen_plotPolytopes(filled_polytopes, (plotFormat),(fillFormat),(fig_num)); %#ok<NASGU>

    % plot the means in black
    temp = zeros(length(filled_polytopes),2);
    for ith_poly = 1:length(filled_polytopes)
        temp(ith_poly,:) = filled_polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3,'DisplayName','mean');

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



