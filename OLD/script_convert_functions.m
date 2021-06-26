%% script_convert_functions
% Converts functions into standard form

% Create a demo file
FuncEGen(1).fileNameSuffix = 'polytopesExpandEvenly';
FuncEGen(1).class = 'MapGen';
FuncEGen(1).short_description = 'Expands an obstacle out by exp_dist on all sides.';
FuncEGen(1).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
FuncEGen(1).filename_main = 'fcn_polytope_editing_expand_polytopes_evenly.m';
FuncEGen(1).filename_script = 'script_polytope_editing_expand_polytopes_evenly.m';


% Input Arguments
FuncEGen(1).Inputs(1).Name = 'polytopes';
FuncEGen(1).Inputs(1).Type = 'polytopes'; % of polytopes type
FuncEGen(1).Inputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(1).Inputs(1).Description = 'the structure of ''polytopes'' type that stores the polytopes to be expanded'; 

FuncEGen(1).Inputs(2).Name = 'delta';
FuncEGen(1).Inputs(2).Type = 'column_of_numbers'; % Forces it to have only one column. See fcn_MapGen_checkInputsToFunctions
FuncEGen(1).Inputs(2).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(1).Inputs(2).Required = 1; % 1 is required, 0 is not
FuncEGen(1).Inputs(2).Description = 'a small number relative to vehicle size to determine the inside of an obstacle'; 

FuncEGen(1).Inputs(3).Name = 'exp_dist';
FuncEGen(1).Inputs(3).Type = 'column_of_numbers'; % Blank means it will not be checked
FuncEGen(1).Inputs(3).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(1).Inputs(3).Required = 1; % 1 is required, 0 is not
FuncEGen(1).Inputs(3).Description = 'distance to expand the obstacle'; 

FuncEGen(1).Inputs(4).Name = 'fig_num';
FuncEGen(1).Inputs(4).Type = 'numeric'; % Must be a number
FuncEGen(1).Inputs(4).Required = 0; % 1 is required, 0 is not
FuncEGen(1).Inputs(4).Description = 'any number that acts somewhat like a figure number output. If given, this forces the variable types to be displayed as output and as well makes the input check process verbose.'; 

FuncEGen(1).N_RequiredInputs = 3;
FuncEGen(1).N_OptionalInputs = 1;

% Dependencies
FuncEGen(1).Dependencies(1).Name = 'fcn_MapGen_polytopeCentroidAndArea';

% Output Arguments
FuncEGen(1).Outputs(1).Name = 'exp_polytopes';
FuncEGen(1).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
FuncEGen(1).Outputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(1).Outputs(1).Description = 'structure of expanded polytopes'; 

% Authorship
FuncEGen(1).Author = 'Seth Tau';
FuncEGen(1).Date = '2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan';
FuncEGen(1).Contact = 'sbrennan@psu.edu and sat5340@psu.edu';

% Common Flags
FuncEGen(1).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
FuncEGen(1).Flags.DoPlot = 1; % Adds the plotting
FuncEGen(1).Flags.CheckInputs = 1; % Adds input checking, turns this on by default


%% 





% Call the function to build the files
fcn_FuncEGen_buildFunctionAndScriptFromDescription(FuncEGen);





function C = ImportLongText(FileName)
S = fileread(FileName);
S(S == char(13)) = [];    % DOS linebreaks to unix linebreaks
S(S == newline) = ' ';   % Replace linebreaks by space - NOTE: this is char(10)
S = strrep(S, '  ', ' '); % Replace double spaces by single space
C = strsplit(S, '§');
end