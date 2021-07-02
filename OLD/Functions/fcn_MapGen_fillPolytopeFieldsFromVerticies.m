function [ ...
filled_polytopes ...
] = ...
fcn_MapGen_fillPolytopeFieldsFromVerticies( ...
polytopes, ...
varargin...
)
% fcn_MapGen_fillPolytopeFieldsFromVerticies
% Given a polytoope structure array where the verticies field is filled, 
% calculates the values for all the other fields.
% 
% 
% 
% FORMAT:
% 
%    [ ...
%    filled_polytopes ...
%    ] = ...
%    fcn_MapGen_fillPolytopeFieldsFromVerticies( ...
%    polytopes, ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     polytopes: an individual structure or structure array of 'polytopes' 
%     type that stores the polytopes to be filled
% 
%     (optional inputs)
%
%     fig_num: any number that acts as a figure number output, causing a 
%     figure to be drawn showing results.
% 
% 
% OUTPUTS:
% 
%     filled_polytopes: the polytopes array with all fields completed
% 
% 
% DEPENDENCIES:
% 
%     fcn_MapGen_polytopeCentroidAndArea
% 
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_fillPolytopeFieldsFromVerticies
% for a full test suite.
% 
% This function was written on 2021_07_02 by Sean Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of function

% 
% TO DO:
% 
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments 
flag_do_plot = 0;      % Set equal to 1 for plotting 
flag_do_debug = 0;     % Set equal to 1 for debugging 

if flag_do_debug
    fig_for_debug = 605;
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


if 1 == flag_check_inputs

    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the polytopes input, make sure it is 'polytopes' type
    fcn_MapGen_checkInputsToFunctions(...
        polytopes, 'polytopes');
 
end

% Does user want to show the plots?
if  2== nargin
    fig_num = varargin{end};
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_for_debug = fig.Number;
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






% Apply the stretch
num_poly = length(polytopes);
for ith_poly = 1:num_poly % pull each polytope
    
    % adjust polytopes
    filled_polytopes(ith_poly).xv        = (polytopes(ith_poly).vertices(1:end-1,1)');
    filled_polytopes(ith_poly).yv        = (polytopes(ith_poly).vertices(1:end-1,2)');
    filled_polytopes(ith_poly).distances = ...
        sum((polytopes(ith_poly).vertices(1:end-1,:) - ...
        polytopes(ith_poly).vertices(2:end,:)).^2,2).^0.5;
    
    % Calculate the mean and area
    [filled_polytopes(ith_poly).mean,filled_polytopes(ith_poly).area] = ...
        fcn_MapGen_polytopeCentroidAndArea(polytopes(ith_poly).vertices);
    
    % Find max radius
    radii = sum(...
        (polytopes(ith_poly).vertices(1:end-1,:) - ...
        ones(length(polytopes(ith_poly).xv),1)*polytopes(ith_poly).mean).^2,2).^0.5;
    filled_polytopes(ith_poly).max_radius = ...
        max(radii);
end

% 

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
    % Nothing to plot here
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

end

