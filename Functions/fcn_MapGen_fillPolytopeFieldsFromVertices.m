function [ ...
filled_polytopes ...
] = ...
fcn_MapGen_fillPolytopeFieldsFromVertices( ...
polytopes, ...
varargin...
)
% fcn_MapGen_fillPolytopeFieldsFromVertices
% Given a polytoope structure array where the vertices field is filled, 
% calculates the values for all the other fields.
% 
% 
% 
% FORMAT:
% 
%    [ ...
%    filled_polytopes ...
%    ] = ...
%    fcn_MapGen_fillPolytopeFieldsFromVertices( ...
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
%     fcn_MapGen_checkInputsToFunctions
% 
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

% 
% TO DO:
% 
% -- fill in to-do items here.

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    fig_for_debug = 3;
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

    % Check the polytopes input, make sure it has vertices   
    if ~isfield(polytopes,'vertices')
        error('Field of vertices was not found');
    end
    
    % Check the vertices input to have 4 or more rows, 2 columns
    %     fcn_MapGen_checkInputsToFunctions(...
    %         polytopes.vertices, '2column_of_numbers',[4 5]);
    
 
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

% Initialize variables
filled_polytopes = polytopes;
num_poly = length(polytopes);

% Loop over each polytope, filling in data for each
for ith_poly = 1:num_poly 
    
    % check that verticies are counter-clockwise by calculating the angles. 
    [angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(...
        filled_polytopes(ith_poly).vertices);

    % Confirm that all angles are positive
    if ~all(angles>=0)
        if any(isnan(angles)) % This happens when there is a repeating point, which is a degenerate poly
            filled_polytopes(ith_poly).vertices = nan(length(filled_polytopes(ith_poly).vertices),2);
        elseif all(angles<=0)
            filled_polytopes(ith_poly).vertices = flipud(filled_polytopes(ith_poly).vertices);
        else            
            fprintf(1,'Verticies:\n');
            for ith_vertex = 1:length(filled_polytopes(ith_poly).vertices)
                fprintf(1,'%.2f %.2f\n',filled_polytopes(ith_poly).vertices(ith_vertex,1),filled_polytopes(ith_poly).vertices(ith_vertex,2))
            end
            fprintf(1,'\nAngles:\n');
            for ith_angle = 1:length(angles)                
                fprintf(1,'%.2f\n',angles(ith_angle));
            end
            error('All vertices must be organized counter-clockwise, e.g. with positive cross-products');
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
        fcn_MapGen_polytopeCentroidAndArea(polytopes(ith_poly).vertices);
    
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
    fcn_MapGen_plotPolytopes(filled_polytopes,fig_num,'b',2);
      
    % plot the means in black
    temp = zeros(length(filled_polytopes),2);
    for ith_poly = 1:length(filled_polytopes)
        temp(ith_poly,:) = filled_polytopes(ith_poly).mean;
    end
    plot(temp(:,1),temp(:,2),'ko','Markersize',3);
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



