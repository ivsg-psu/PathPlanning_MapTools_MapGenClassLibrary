function [cropped_vertices] = ...
    fcn_MapGen_polytopeRemoveColinearVertices(...,
    vertices,...
    varargin)

% Given a set of vertices that may have multiple points in a line on the
% same edge, this function returns an equivalent polytope with no repeats
% on any edge.
%
% FORMAT:
% 
% [cropped_vertices] = ...
%    fcn_MapGen_polytopeRemoveColinearVertices(...
%     vertices,...
%     (fig_num))
%
% INPUTS:
%
%     vertices: an Nx2 matrix of [x y] vertices
%
%    (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     cropped_vertices: an Mx2 matrix of [x y] points, where no vertices
%     are repeated
%   
% DEPENDENCIES:
% 
%     fcn_DebugTools_checkInputsToFunctions
% 
% EXAMPLES:
%      
% For additional examples, see:
% script_test_fcn_MapGen_polytopeRemoveColinearVertices
%
% This function was written on 2021_08_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu 
%

% Revision History:
% 2021_08_03 - S. Brennan
% -- first write of code
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions


% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2column_of_numbers');

    end
end
    

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  2 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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


% Remove repeats
[cropped_vertices,~,~] = unique(vertices,'rows','stable');

% Use the cross-product to eliminate co-linear points
Npoints = length(cropped_vertices(:,1));
good_indices = zeros(Npoints,1);

for ith_point = 1:Npoints
    
    if ith_point>1
        previous_vector = [cropped_vertices(ith_point,:)-cropped_vertices(ith_point-1,:), 0];
    else % must wrap around
        previous_vector = [cropped_vertices(end,:)-cropped_vertices(1,:), 0];
    end
    
    if ith_point<Npoints
        subsequent_vector = [cropped_vertices(ith_point+1,:)-cropped_vertices(ith_point,:), 0];
    else % must wrap around
        subsequent_vector = [cropped_vertices(1,:)-cropped_vertices(Npoints,:), 0];
    end
    
    cross_result = cross(previous_vector,subsequent_vector);
    
    if cross_result(1,3)~=0
        good_indices(ith_point,:) = 1;
    end
end

% Final polytope
cropped_vertices = cropped_vertices(good_indices>0,:);



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
    grid on
    grid minor
    hold on
    axis equal
    
    % Plot the input vertices
    plot(vertices(:,1),vertices(:,2),'r-');
    
    % Plot the cropped_vertices
    plot(cropped_vertices(:,1),cropped_vertices(:,2),'ko');    

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


