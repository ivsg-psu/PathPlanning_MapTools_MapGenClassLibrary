function [ ...
varargout...
] = ...
fcn_MapGen_checkInputsToFunctions( ...
variable, ...
variable_type_string, ...
varargin...
)
% fcn_MapGen_checkInputsToFunctions
% Checks the variable types commonly used in the FuncE codes to ensure 
% they are correctly formed.
% 
% This is a template function which is built for each class of functions. 
% It is typically called at the top of most functions in a particular 
% class. The input is a variable and a string defining the "type" of the 
% variable. This function checks to see that they are compatible. For 
% example, say there is a 'column_vector' type of variables used in the 
% function that is always a N x 1 array; if someone had a variable called 
% "test_example", one could check that this fit the 'column_vector' type 
% by calling 
% fcn_MapGen_checkInputsToFunctions(test_example,'column_vector'). This 
% function would then check that the array was N x 1, and if it was not, 
% it would send out an error warning.
% 
% FORMAT:
% 
%    [ ...
%    (AllowableInputs) ...
%    ] = ...
%    fcn_MapGen_checkInputsToFunctions( ...
%    variable, ...
%    variable_type_string, ...
%    (required_length), ...
%    (fig_num) ...
%    )
% 
% INPUTS:
% 
%     variable: the variable to check
% 
%     variable_type_string: a string representing the variable type to 
%     check. Call the function with any figure number to see allowable 
%     options.
% 
%     (optional inputs)
%
%     required_length: an integer forcing the value of N, giving an error 
%     if the input variable does not have length N. Another optional input 
%     is a rwo vector [A B] where, if B is greater than A, then the vector 
%     must be A or longer. If B is less than A, then the vector must be A 
%     or shorter. If B = A, then the vector must be length A, and no 
%     shorter or greater.
% 
%     fig_num: any number that acts somewhat like a figure number output. 
%     If given, this forces the variable types to be displayed as output 
%     and as well makes the input check process verbose.
% 
% 
% OUTPUTS:
% 
%     (optional outputs)
%
%     AllowableInputs: This is a structure output that lists all the 
%     allowable types, and a description of each. As well, if the output 
%     argument is given, the same information is printed within the 
%     workspace.
% 
% 
% DEPENDENCIES:
% 
%    (none)
% 
% EXAMPLES:
% 
% See the script: script_test_fcn_MapGen_checkInputsToFunctions
% for a full test suite.
% 
% This function was written on 2021_06_20 by S. Brennan
% Questions or comments? contact sbrennan@psu.edu

% 
% REVISION HISTORY:
% 
% 2021_06_20 by S. Brennan
% -- first write of function
% 2021_07_07 by S. Brennan
% -- modified to allow general input types, e.g. 8column_of_integers
% 
% TO DO:
% 
% -- fill in to-do items here.

return;

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments 
flag_do_plot = 0;      % Set equal to 1 for plotting 
flag_do_debug = 0;     % Set equal to 1 for debugging 

if flag_do_debug
    fig_for_debug = 159;
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
    if nargin < 2 || nargin > 4
        error('Incorrect number of input arguments')
    end

    % Check the variable_type_string input, make sure it is characters
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§


% Grab the variable name
variable_name = inputname(1);

% Check to see if output argument given
if nargout>0
    varargout{1}= INTERNAL_fcn_showPossibleFields;
    return
end

% Set default flags all to "off" mode
flags = INTERNAL_fcn_setDefaultFlagsToOff;

% See if special inputs:
if nargin == 3
    flags.check_requiredRowLength = 1;     % Must check for required length
    flags.rowLengthRangeRequired = varargin{1}; % Set to [x y]. Variable must be x or greater if y>x, =x if y=x, x or less if y<x
end

% Grab flag settings for current input
flags = INTERNAL_fcn_setFlagsByType(flags,variable_type_string);

% Check that variable meets requirements
INTERNAL_confirmVariable(flags,variable,variable_name);


%

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


function flags = INTERNAL_fcn_setDefaultFlagsToOff
% Set default flags all to "off" mode
flags.check_if_isnumeric = 0; % Check to see if isnumeric
flags.check_if_strictly_positive = 0; % Check to see if number is greater than zero, and not zero
flags.check_required_columns  = 0; % Check the number of columns
flags.minNrequiredcolumns  = 0; % No check
flags.maxNrequiredcolumns  = 0; % No check
flags.check_if_noNaN = 0; % Check that there are no NaN 
flags.check_if_integer = 0; % Check that the variable is an integer
flags.check_requiredRowLength = 0;     % Don't check for required length
flags.rowLengthRangeRequired = [0 0]; % Set to [x y]. Variable must be x or greater if y>x, =x if y=x, x or less if y<x
flags.check_likeStructure = 0; % Check that result is like a particular structure
template_structure = ...
    struct(...
    'vertices',[],...
    'xv',[],...
    'yv',[],...
    'distances',[],...
    'mean',[],...
    'area',[],...
    'max_radius',[]);
flags.structureToBeLike = template_structure;
end

%%
function flags = INTERNAL_fcn_setFlagsByType(flags, variable_type_string)

flag_pattern_was_matched = 0;

% Check the "Ncolumn_of" pattern, where N is a digit
pattern = digitsPattern(1)+"column_of";
if contains(variable_type_string,pattern)
    match = extract(variable_type_string,pattern);
    string_result = match{1};
    ncols_max = str2double(string_result(1));

    % Check for NorMcolumn_of format
    pattern = digitsPattern(1)+"or"+digitsPattern(1)+"column_of";
    if contains(variable_type_string,pattern)
        match = extract(variable_type_string,pattern);
        string_result = match{1};
        ncols_min = str2double(string_result(1));
    else
        ncols_min = ncols_max;
    end
    
    flags.check_if_isnumeric = 1; % Must be a number
    flags.check_required_columns  = 1; % Check the number of columns
    flags.minNrequiredcolumns  = ncols_min; % Must be 1 columns
    flags.maxNrequiredcolumns  = ncols_max; % Must be 1 columns
    flags.check_if_noNaN = 1; % Check that there are no NaN
    
    flag_pattern_was_matched = 1;
end

% positive_XXX
pattern = 'positive_';
if contains(variable_type_string,pattern)
    flags.check_if_strictly_positive = 1; % Must be a number    
end

% XXX_of_integers
pattern = '_of_integers';
if contains(variable_type_string,pattern)
    flags.check_if_integer = 1; % Check that the variable is an integer
end

% XXX_of_mixed
pattern = '_of_mixed';
if contains(variable_type_string,pattern)
    flags.check_if_noNaN = 0; % Removes check that it be numeric
end

% polytopes
if strcmpi(variable_type_string,'polytopes')
    flags.check_likeStructure = 1; % Check that result is like a particular structure
    template_structure = ...
        struct(...
        'vertices',[],...
        'xv',[],...
        'yv',[],...
        'distances',[],...
        'mean',[],...
        'area',[],...
        'max_radius',[]);
    flags.structureToBeLike = template_structure;

    flag_pattern_was_matched = 1;
end

% mixedset
if strcmpi(variable_type_string,'mixedset')
    flags.check_likeStructure = 1; % Check that result is like a particular structure
    template_structure = ...
        struct(...
        'name',[],...
        'settings',[],...
        'AABB',[]);
    flags.structureToBeLike = template_structure;

    flag_pattern_was_matched = 1;
end

if 0==flag_pattern_was_matched
    error('The variable type: %s is not defined in context of error checking.',variable_type_string);
end

end % Ends INTERNAL_fcn_setFlagsByType

function INTERNAL_confirmVariable(flags,variable,variable_name)

% Numeric?
if flags.check_if_isnumeric   
    if ~isnumeric(variable)
        error('The %s input must be numeric.',variable_name);
    end
end

% Strictly positive?
if flags.check_if_strictly_positive   
    if any(variable<=0)
        error('The %s input must be strictly positive, e.g. greater than zero and not equal to zero.',variable_name);
    end
end

% NaN?
if flags.check_if_noNaN   
    if any(isnan(variable),'all')
        error('The %s input must have no NaN values.',variable_name);
    end
end

% Integer?
if flags.check_if_integer   
    if ~all(round(variable)==variable)
        error('The %s input must be an integer.',variable_name);
    end
end

% Column length?
if flags.check_required_columns    
    if flags.minNrequiredcolumns==0
        error('Need to set minimum number of columns for variable type: %s.',variable_name);
    end
    if flags.maxNrequiredcolumns==0
        error('Need to set maximum number of columns for variable type: %s.',variable_name);
    end
    
    % Exactly a number of columns?
    if flags.minNrequiredcolumns==flags.maxNrequiredcolumns
        if length(variable(1,:))~=flags.minNrequiredcolumns
            error('The %s input must have exactly %.0d columns.',variable_name,flags.minNrequiredcolumns);
        end
    end
    
    % A minimum number of columns
    if length(variable(1,:))<flags.minNrequiredcolumns
        error('The %s input must have at least %.0d columns.',variable_name,flags.minNrequiredcolumns);
    end

    % A maximum number of columns
    if length(variable(1,:))>flags.maxNrequiredcolumns
        error('The %s input must have no more than %.0d columns.',variable_name,flags.maxNrequiredcolumns);
    end
    
end

% Row length?
if flags.check_requiredRowLength    
    required_length = flags.rowLengthRangeRequired;
    if length(required_length(1,:))==1  % Exact, given number of rows
        if length(variable(:,1))~=required_length
            error('The %s input must have exactly %.0d rows',variable_name,required_length);
        end
    else
        if required_length(1,2)>required_length(1,1) % Must be at least given number of rows, or more
            min_length = required_length(1,1);
            if length(variable(:,1))<min_length
                error('The %s input must have %.0d rows or more',variable_name,min_length);
            end
        elseif required_length(1,2)<required_length(1,1) % Must be no more than given number of rows
            max_length = required_length(1,1);
            if length(variable(:,1))>max_length
                error('The %s input must have no more than %.0d rows',variable_name,max_length);
            end
        else % It has to be equal
            required_length = required_length(1,1);
            if length(variable(:,1))~=required_length
                error('The %s input must have %.0d rows or more',variable_name,required_length);
            end
        end
    end
end

% Structure?
if flags.check_likeStructure 
    template_fields = fieldnames(orderfields(flags.structureToBeLike));
    reference_fields = fieldnames(orderfields(variable));
    if ~isequal(template_fields,reference_fields)
        fprintf(1,'The template has fields of:\n');
        for ith_field = 1:length(template_fields)
            fprintf(1,'\t %s\n',string(template_fields(ith_field)));
        end
        fprintf(1,'The %s input has fields of:\n',variable_name);
        for ith_field = 1:length(reference_fields)
            fprintf(1,'\t %s\n',string(reference_fields(ith_field)));
        end
        error('The %s input must be a structure type. All the fields must be match the reference structure.',variable_name);
    end
    
end 
    

end % Ends INTERNAL_confirmVariable

function allowable_inputs = INTERNAL_fcn_showPossibleFields
num_inputs = 0;

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = '1column_of_numbers';
allowable_inputs(num_inputs).description = 'checks that the input type is N x 1 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N.';

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = 'positive_1column_of_numbers';
allowable_inputs(num_inputs).description = 'checks that the input type is N x 1 and is a strictly positive number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N.';

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = '2column_of_numbers';
allowable_inputs(num_inputs).description = 'checks that the input type is N x 2 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.';

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = '4column_of_numbers';
allowable_inputs(num_inputs).description = 'checks that the input type is N x 4 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.';

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = '2or3column_of_numbers';
allowable_inputs(num_inputs).description = 'checks that the input type is N x 2 or N x 3 and is a number. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.';

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = '2column_of_integers';
allowable_inputs(num_inputs).description = 'checks that the input type is N x 2 and is an integer. Optional input: an integer forcing the value of N, giving an error if the input variable does not have length N. Another optional input is a rwo vector [A B] where, if B is greater than A, then the vector must be A or longer. If B is less than A, then the vector must be A or shorter. If B = A, then the vector must be length A, and no shorter or greater.';

num_inputs = num_inputs+1;
allowable_inputs(num_inputs).name = 'polytopes';
allowable_inputs(num_inputs).description = 'a 1-by-n seven field structure of polytopes within the boundaries, where n <= number of polytopes with fields:  vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is the number of the individual polytope vertices  xv: a 1-by-m vector of vertice x-coordinates  yv: a 1-by-m vector of vertice y-coordinates  distances: a 1-by-m vector of perimeter distances from one point to the next point, distances(i) = distance from vertices(i) to vertices(i+1) mean: centroid xy coordinate of the polytope area: area of the polytope max_radius: the largest distance from the centroid to any vertex';
end



