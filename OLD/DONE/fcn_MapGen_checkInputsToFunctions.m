function fcn_MapGen_checkInputsToFunctions(...
    variable,variable_type_string,varargin)

% fcn_MapGen_checkInputsToFunctions
% Checks the variable types commonly used in the MapGen codes to
% ensure they are correctly formed.
%
% This function is typically called at the top of most functions. The input
% is a variable and a string defining the "type" of the variable. This
% function checks to see that they are compatible. For example, say there
% 'column_vector' type of variables used in the function that is always a N
% x 1 array; if someone had a variable called "test_example", they could
% check that this fit the 'column_vector' type by calling
% fcn_MapGen_checkInputsToFunctions(test_example,'column_vector').
% This function would then check that the array was N x 1, and if it was
% not, it would send out an error warning.
%
% FORMAT:
%
%      fcn_MapGen_checkInputsToFunctions(...
%      variable,variable_type_string,(optional_arguments))
%
% INPUTS:
%
%      variable: the variable to check
%
%      variable_type_string: a string representing the variable type to
%      check. The current strings include:
%
%            'column_of_numbers' - checks that the input type is N x 1 and
%            is a number. Optional input: an integer forcing the value
%            of N, giving an error if the input variable does not have
%            length N.
%
%            '2column_of_numbers' - checks that the input type is N x 2 and
%            is a number. Optional input: an integer forcing the value
%            of N, giving an error if the input variable does not have
%            length N. Another optional input is a rwo vector [A B] where,
%            if B is greater than A, then the vector must be A or longer.
%            If B is less than A, then the vector must be A or shorter. If
%            B = A, then the vector must be length A, and no shorter or
%            greater.
%
%            '4column_of_numbers' - checks that the input type is N x 4 and
%            is a number. Optional input: an integer forcing the value
%            of N, giving an error if the input variable does not have
%            length N. Another optional input is a rwo vector [A B] where,
%            if B is greater than A, then the vector must be A or longer.
%            If B is less than A, then the vector must be A or shorter. If
%            B = A, then the vector must be length A, and no shorter or
%            greater.
%
%            '2or3column_of_numbers'  - checks that the input type is N x 2
%            or N x 3 and is a number. Optional input: an integer forcing
%            the value of N, giving an error if the input variable does not
%            have length N. Another optional input is a rwo vector [A B]
%            where, if B is greater than A, then the vector must be A or
%            longer. If B is less than A, then the vector must be A or
%            shorter. If B = A, then the vector must be length A, and no
%            shorter or greater.
%
%            '2column_of_integers' - checks that the input type is N x 2
%            and is an integer. Optional input: an integer forcing the
%            value of N, giving an error if the input variable does not
%            have length N. Another optional input is a rwo vector [A B]
%            where, if B is greater than A, then the vector must be A or
%            longer. If B is less than A, then the vector must be A or
%            shorter. If B = A, then the vector must be length A, and no
%            shorter or greater.
%
%             'polytopes' - a 1-by-n seven field structure of polytopes
%             within the boundaries, where n <= number of polytopes with
%             fields: 
%             vertices: a m+1-by-2 matrix of xy points with row1 =
%             rowm+1, where m is the number of the individual polytope
%             vertices 
%             xv: a 1-by-m vector of vertice x-coordinates 
%             yv: a 1-by-m vector of vertice y-coordinates 
%             distances: a 1-by-m vector of perimeter distances from one
%             point to the next point, distances(i) = distance from
%             vertices(i) to vertices(i+1)
%             mean: centroid xy coordinate of the polytope
%             area: area of the polytope
%             max_radius: the largest distance from the centroid to any
%             vertex
%
%
%      Note that the variable_type_string is not case sensitive: for
%      example, 'station' and 'Station' or 'STAtion' all give the same
%      result.
%
% OUTPUTS:
%
%      No explicit outputs, but produces MATLAB error outputs if conditions
%      not met, with explanation within the error outputs of the problem.
%
% EXAMPLES:
%
% See the script: script_test_fcn_MapGen_checkInputsToFunctions
% for a full test suite.
%
% DEPENDENCIES:
%
%      Uses MATLABs dbstack feature to trace dependencies 
%
% This function was written on 2021_01_06 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2021_06_06:
%      -- first write of the code copying functionality from fcn_Path_checkInputsToFunctions

% TO-DO:
% -- finish polytope checks for completeness
% -- add polytopes and 4column_of_numbers to test script

flag_do_debug = 0; % Flag to debug the results
flag_do_plot = 0; % Flag to plot the results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
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

if flag_check_inputs == 1   
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the string input, make sure it is characters
    if ~ischar(variable_type_string)
        error('The variable_type_string input must be a string type, for example: ''Path'' ');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab the variable name
variable_name = inputname(1);


%% column_of_numbers
if strcmpi(variable_type_string,'column_of_numbers')
    % Check the station input
    if (length(variable(1,:))~=1) || ~isnumeric(variable)
        error('The %s input must be a column vector type, namely an N x 1 vector with N>=1',variable_name);
    end
    
    if any(isnan(variable),'all')
        error('The %s input must be a column vector type, namely an N x 1 vector that has no NaN values.',variable_name);
    end
    
    if nargin == 3
        required_length = varargin{1};
        if length(variable(:,1))~=required_length
            error('The %s input must be a column vector (N x 1) with N == %.0d',variable_name,required_length);
        end
    end
end

%% 2column_of_numbers
if strcmpi(variable_type_string,'2column_of_numbers')
    % Check the station input
    if (length(variable(1,:))~=2) || ~isnumeric(variable)
        error('The %s input must be a 2column_of_numbers type, namely an N x 2 vector with N>=1',variable_name);
    end
    
    if any(isnan(variable),'all')
        error('The %s input must be a 2column_of_numbers type, namely an N x 2 vector that has no NaN values.',variable_name);
    end
    
    if nargin == 3
        required_length = varargin{1};
        if length(required_length(1,:))==1
            if length(variable(:,1))~=required_length
                error('The %s input must be a 2column_of_numbers, namely (N x 2) with N == %.0d',variable_name,required_length);
            end
        else
            if required_length(1,2)>required_length(1,1)
                min_length = required_length(1,1);
                if length(variable(:,1))<min_length
                    error('The %s input must be a 2column_of_numbers, namely (N x 2) with N >= %.0d',variable_name,min_length);
                end
            elseif required_length(1,2)<required_length(1,1)
                max_length = required_length(1,1);
                if length(variable(:,1))>max_length
                    error('The %s input must be a 2column_of_numbers, namely (N x 2) with N <= %.0d',variable_name,max_length);
                end
            else % It has to be equal
                required_length = required_length(1,1);
                if length(variable(:,1))~=required_length
                    error('The %s input must be a 2column_of_numbers, namely (N x 2) with N = %.0d',variable_name,required_length);
                end
            end
        end
    end
end

%% 4column_of_numbers
if strcmpi(variable_type_string,'4column_of_numbers')
    % Check the length of input
    if (length(variable(1,:))~=4) || ~isnumeric(variable)
        error('The %s input must be a 4column_of_numbers type, namely an N x 4 vector with N>=1',variable_name);
    end
    
    if any(isnan(variable),'all')
        error('The %s input must be a 4column_of_numbers type, namely an N x 4 vector that has no NaN values.',variable_name);
    end
    
    if nargin == 3
        required_length = varargin{1};
        if length(required_length(1,:))==1
            if length(variable(:,1))~=required_length
                error('The %s input must be a 4column_of_numbers, namely (N x 4) with N == %.0d',variable_name,required_length);
            end
        else
            if required_length(1,2)>required_length(1,1)
                min_length = required_length(1,1);
                if length(variable(:,1))<min_length
                    error('The %s input must be a 4column_of_numbers, namely (N x 4) with N >= %.0d',variable_name,min_length);
                end
            elseif required_length(1,2)<required_length(1,1)
                max_length = required_length(1,1);
                if length(variable(:,1))>max_length
                    error('The %s input must be a 4column_of_numbers, namely (N x 4) with N <= %.0d',variable_name,max_length);
                end
            else % It has to be equal
                required_length = required_length(1,1);
                if length(variable(:,1))~=required_length
                    error('The %s input must be a 4column_of_numbers, namely (N x 4) with N = %.0d',variable_name,required_length);
                end
            end
        end
    end
end

%% 2or3column_of_numbers
if strcmpi(variable_type_string,'2or3column_of_numbers')
    % Check the station input
    if ((length(variable(1,:))<2) || (length(variable(1,:))>3)) || ~isnumeric(variable)
        error('The %s input must be a 2or3column_of_numbers type, namely an N x 2 or N x 3 vector with N>=1',variable_name);
    end
    
    if any(isnan(variable),'all')
        error('The %s input must be a 2or3column_of_numbers type, namely an N x 2 or N x 3 vector that has no NaN values.',variable_name);
    end
    
    % Is there a specified required length?
    if nargin == 3
        required_length = varargin{1};
        if length(required_length(1,:))==1
            if length(variable(:,1))~=required_length
                error('The %s input must be a 2or3column_of_numbers, namely (N x 2) or (N x 3) with N == %.0d',variable_name,required_length);
            end
        else
            if required_length(1,2)>required_length(1,1)
                min_length = required_length(1,1);
                if length(variable(:,1))<min_length
                    error('The %s input must be a 2or3column_of_numbers, namely (N x 2) or (N x 3) with N >= %.0d',variable_name,min_length);
                end
            elseif required_length(1,2)<required_length(1,1)
                max_length = required_length(1,1);
                if length(variable(:,1))>max_length
                    error('The %s input must be a 2or3column_of_numbers, namely (N x 2) or (N x 3) with N <= %.0d',variable_name,max_length);
                end
            else % It has to be equal
                required_length = required_length(1,1);
                if length(variable(:,1))~=required_length
                    error('The %s input must be a 2or3column_of_numbers, namely (N x 2) or (N x 3) with N = %.0d',variable_name,required_length);
                end
            end
        end
    end
end

%% 2column_of_integers
if strcmpi(variable_type_string,'2column_of_integers')
    % Check the station input
    if (length(variable(1,:))~=2) 
        error('The %s input must be a 2column_of_integers type, namely an N x 2 vector with N>=1. The number of columns is not 2.',variable_name);
    end
    
    if ~all(round(variable)==variable)
        error('The %s input must be a 2column_of_integers type, namely an N x 2 vector with N>=1. A non-integer type was found.',variable_name);
    end
    
    if any(isnan(variable),'all')
        error('The %s input must be a 2column_of_integers type, namely an N x 2 vector that has no NaN values.',variable_name);
    end
    
    if nargin == 3
        required_length = varargin{1};
        if length(required_length(1,:))==1
            if length(variable(:,1))~=required_length
                error('The %s input must be a 2column_of_integers, namely (N x 2) with N == %.0d',variable_name,required_length);
            end
        else
            if required_length(1,2)>required_length(1,1)
                min_length = required_length(1,1);
                if length(variable(:,1))<min_length
                    error('The %s input must be a 2column_of_integers, namely (N x 2) with N >= %.0d',variable_name,min_length);
                end
            elseif required_length(1,2)<required_length(1,1)
                max_length = required_length(1,1);
                if length(variable(:,1))>max_length
                    error('The %s input must be a 2column_of_integers, namely (N x 2) with N <= %.0d',variable_name,max_length);
                end
            else % It has to be equal
                required_length = required_length(1,1);
                if length(variable(:,1))~=required_length
                    error('The %s input must be a 2column_of_integers, namely (N x 2) with N = %.0d',variable_name,required_length);
                end
            end
        end
    end
end

%% polytopes
if strcmpi(variable_type_string,'polytopes')
    % Check the polytopes subfields exist
    %             vertices: a m+1-by-2 matrix of xy points with row1 =
    %             rowm+1, where m is the number of the individual polytope
    %             vertices
    %             xv: a 1-by-m vector of vertice x-coordinates
    %             yv: a 1-by-m vector of vertice y-coordinates
    %             distances: a 1-by-m vector of perimeter distances from one
    %             point to the next point, distances(i) = distance from
    %             vertices(i) to vertices(i+1)
    %             mean: centroid xy coordinate of the polytope
    %             area: area of the polytope
    try
        vertices        = variable(1).vertices;
        xv              = variable(1).xv;
        yv              = variable(1).yv;
        distances       = variable(1).distances;
        mean            = variable(1).mean;
        area            = variable(1).area;
        max_radius      = variable(1).max_radius;
    catch
        error('The %s input must be a polytopes type, namely being a structure with fields: vertices, xv, yv, distances, mean, area, max_radius. One or more of the fields were not found. ',variable_name);
    end
    
    % Check that all are numeric
    if  ...
            ~isnumeric(vertices) || ...
            ~isnumeric(xv) || ...
            ~isnumeric(yv) ||  ...
            ~isnumeric(distances) || ...
            ~isnumeric(mean) || ...
            ~isnumeric(area) || ...
            ~isnumeric(max_radius)
        error('The %s input must be a polytope type. All the fields must be numeric, but at least one data field is non-numeric.',variable_name);
    end
    
    %     % Check that all are 1-dimensional columns
    %     if (length(X_central(1,:))~=1) || (length(Y_central(1,:))~=1) || (length(Z_central(1,:))~=1) ||  (length(Station_central(1,:))~=1)
    %         error('The %s input must be a traversal type, namely a structure with fields X, Y, Z, and Station, each N x 1 numeric arrays. At least one data field has multiple columns.',variable_name);
    %     end
    %
    %     % Check that their lengths are all the same
    %     if (length(X_central(:,1))~=length(Y_central(:,1))) || ((length(X_central(:,1))~=length(Z_central(:,1))))  || ((length(X_central(:,1))~=length(Station_central(:,1))))
    %         error('The %s input must be a traversal type, namely a structure with fields X, Y, Z, and Station, each N x 1 numeric arrays. The lengths do not match.',variable_name);
    %     end
    %
    %     % Make sure the station field is sorted
    %     if ~issorted(Station_central,'strictascend')
    %         error('The Station field on the %s input must be strictly increasing',variable_name);
    %     end
end


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
    fprintf(1,'The variable: %s was checked that it meets type: %s, and no errors were detected.\n',variable_name,variable_type_string);
end
if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends the function

