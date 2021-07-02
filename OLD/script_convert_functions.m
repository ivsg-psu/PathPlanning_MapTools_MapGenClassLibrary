%% script_convert_functions
% Converts functions into standard form


%% Create a polytopesExpandEvenly file
filenum = 1;

FuncEGen(filenum).fileNameSuffix = 'polytopesExpandEvenly';
FuncEGen(filenum).class = 'MapGen';
FuncEGen(filenum).short_description = 'Expands an obstacle out by exp_dist on all sides.';
FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
FuncEGen(filenum).filename_main = 'codeCore_fcn_polytope_editing_expand_polytopes_evenly.m';
FuncEGen(filenum).filename_script = 'codeCore_script_polytope_editing_expand_polytopes_evenly.m';


% Input Arguments
FuncEGen(filenum).Inputs(1).Name = 'polytopes';
FuncEGen(filenum).Inputs(1).Type = 'polytopes'; % of polytopes type
FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(1).Description = 'the structure of ''polytopes'' type that stores the polytopes to be expanded'; 

FuncEGen(filenum).Inputs(2).Name = 'delta';
FuncEGen(filenum).Inputs(2).Type = 'column_of_numbers'; % Forces it to have only one column. See fcn_MapGen_checkInputsToFunctions
FuncEGen(filenum).Inputs(2).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(2).Description = 'a small number relative to vehicle size to determine the inside of an obstacle'; 

FuncEGen(filenum).Inputs(3).Name = 'exp_dist';
FuncEGen(filenum).Inputs(3).Type = 'column_of_numbers'; % Blank means it will not be checked
FuncEGen(filenum).Inputs(3).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(filenum).Inputs(3).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(3).Description = 'distance to expand the obstacle'; 

FuncEGen(filenum).Inputs(4).Name = 'fig_num';
FuncEGen(filenum).Inputs(4).Type = 'numeric'; % Must be a number
FuncEGen(filenum).Inputs(4).Required = 0; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(4).Description = 'any number that acts somewhat like a figure number output. If given, this forces the variable types to be displayed as output and as well makes the input check process verbose.'; 

FuncEGen(filenum).N_RequiredInputs = 3;
FuncEGen(filenum).N_OptionalInputs = 1;

% Dependencies
FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_polytopeCentroidAndArea';

% Output Arguments
FuncEGen(filenum).Outputs(1).Name = 'exp_polytopes';
FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Outputs(1).Description = 'structure of expanded polytopes'; 

% Authorship
FuncEGen(filenum).Author = 'Seth Tau';
FuncEGen(filenum).Date = '2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan';
FuncEGen(filenum).Contact = 'sbrennan@psu.edu and sat5340@psu.edu';

% Common Flags
FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default


%% Create generateOneRandomPolytope file
filenum = 2;


% Create a polytopesExpandEvenly file
FuncEGen(filenum).fileNameSuffix = 'generateOneRandomPolytope';
FuncEGen(filenum).class = 'MapGen';
FuncEGen(filenum).short_description = 'Generates a single random polytope by selecting one randomly from the Halton Set Voronoi method (HSV) tiling.';
FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
FuncEGen(filenum).filename_main = 'codeCore_generateOneRandomPolytope.m';
FuncEGen(filenum).filename_script = ''; % An empty starter script will be created


% Input Arguments
FuncEGen(filenum).Inputs(1).Name = 'fig_num';
FuncEGen(filenum).Inputs(1).Type = 'numeric'; % Must be a number
FuncEGen(filenum).Inputs(1).Required = 0; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(1).Description = 'any number that acts somewhat like a figure number output. If given, this forces the variable types to be displayed as output and as well makes the input check process verbose.'; 

FuncEGen(filenum).N_RequiredInputs = 0;
FuncEGen(filenum).N_OptionalInputs = 1;

% Dependencies
FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_haltonVoronoiTiling';
FuncEGen(filenum).Dependencies(2).Name = 'fcn_MapGen_polytopeCropEdges';
FuncEGen(filenum).Dependencies(3).Name = 'fcn_MapGen_plotPolytopes';

% Output Arguments
FuncEGen(filenum).Outputs(1).Name = 'one_polytope';
FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Outputs(1).Description = 'one randomly generated polytope'; 

% Authorship
FuncEGen(filenum).Author = 'Sean Brennan';
FuncEGen(filenum).Date = '2021_06_27';
FuncEGen(filenum).Contact = 'sbrennan@psu.edu';

% Common Flags
FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default



%% Create snapToAABB file
filenum = 3;


% Create a snapToAABB file
FuncEGen(filenum).fileNameSuffix = 'snapToAABB';
FuncEGen(filenum).class = 'MapGen';
FuncEGen(filenum).short_description = 'Given an axis-aligned bounding box (AABB), and a test point, returns a snap point representing the contact point on the closest wall to the test point.';
FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
FuncEGen(filenum).filename_main = 'codeCore_snapToAABB.m';
FuncEGen(filenum).filename_script = 'codeCore_snapToAABB_script.m'; 

% Input Arguments
FuncEGen(filenum).Inputs(1).Name = 'axis_aligned_bounding_box';
FuncEGen(filenum).Inputs(1).Type = '4column_of_numbers'; % 4 columns of numbers
FuncEGen(filenum).Inputs(1).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(1).Description = 'the axis-aligned bounding box, in format [xmin ymin xmax ymax]'; 

FuncEGen(filenum).Inputs(2).Name = 'test_point';
FuncEGen(filenum).Inputs(2).Type = '2column_of_numbers'; % 2 columns of numbers
FuncEGen(filenum).Inputs(2).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(2).Description = 'the test point, in format [x y]'; 

FuncEGen(filenum).Inputs(3).Name = 'fig_num';
FuncEGen(filenum).Inputs(3).Type = 'numeric'; % Must be a number
FuncEGen(filenum).Inputs(3).Required = 0; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(3).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 

FuncEGen(filenum).N_RequiredInputs = 2;
FuncEGen(filenum).N_OptionalInputs = 1;

% Dependencies
FuncEGen(filenum).Dependencies(1).Name = '(none)';

% Output Arguments
FuncEGen(filenum).Outputs(1).Name = 'snap_point';
FuncEGen(filenum).Outputs(1).Type = '2column_of_numbers'; % See fcn_MapGen_checkInputsToFunctions
FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Outputs(1).Description = 'the resulting snap point, in format [x y]'; 

% Authorship
FuncEGen(filenum).Author = 'Sean Brennan';
FuncEGen(filenum).Date = '2021_07_02';
FuncEGen(filenum).Contact = 'sbrennan@psu.edu';

% Common Flags
FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default



% Call the function to build the files
fcn_FuncEGen_buildFunctionAndScriptFromDescription(FuncEGen);





function C = ImportLongText(FileName)
S = fileread(FileName);
S(S == char(13)) = [];    % DOS linebreaks to unix linebreaks
S(S == newline) = ' ';   % Replace linebreaks by space - NOTE: this is char(10)
S = strrep(S, '  ', ' '); % Replace double spaces by single space
C = strsplit(S, '§');
end