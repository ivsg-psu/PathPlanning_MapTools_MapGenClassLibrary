%% script_convert_functions
% Converts functions into standard form

clear FuncEGen

%% Create a polytopesExpandEvenly file
filenum = 1;

FuncEGen(filenum).fileNameSuffix = 'isWithinABBB';
FuncEGen(filenum).class = 'MapGen';
FuncEGen(filenum).short_description = 'Checks if the points are within the AABB, returning a vector of 1'' or 0''s the same length as the nubmer of rows of points.';
FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
FuncEGen(filenum).filename_main = 'codeCore_isWithinABBB.m';
FuncEGen(filenum).filename_script = 'codeCore_isWithinABBB_script.m';


% Input Arguments
FuncEGen(filenum).Inputs(1).Name = 'AABB';
FuncEGen(filenum).Inputs(1).Type = '4column_of_numbers'; % of polytopes type
FuncEGen(filenum).Inputs(1).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(1).Description = 'the Axis-Aligned Bounding Box, defined in form of [xmin ymin xmax ymax]'; 

FuncEGen(filenum).Inputs(2).Name = 'test_points';
FuncEGen(filenum).Inputs(2).Type = '2column_of_numbers'; % Forces it to have only one column. See fcn_MapGen_checkInputsToFunctions
FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(2).Description = 'the test points to check, in form of [x y] where x and y are scalar or column vectors'; 

FuncEGen(filenum).Inputs(3).Name = 'fig_num';
FuncEGen(filenum).Inputs(3).Type = 'numeric'; % Must be a number
FuncEGen(filenum).Inputs(3).Required = 0; % 1 is required, 0 is not
FuncEGen(filenum).Inputs(3).Description = 'any number that acts somewhat like a figure number output. If given, this forces the variable types to be displayed as output and as well makes the input check process verbose.'; 

FuncEGen(filenum).N_RequiredInputs = 2;
FuncEGen(filenum).N_OptionalInputs = 1;

% Dependencies
FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_checkInputsToFunctions';

% Output Arguments
FuncEGen(filenum).Outputs(1).Name = 'isInside';
FuncEGen(filenum).Outputs(1).Type = '1column_of_numbers'; % See fcn_MapGen_checkInputsToFunctions
FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
FuncEGen(filenum).Outputs(1).Description = 'a column of 1''s or 0''s, one for each test point, with 1 meaning that the test point is within the AABB'; 

% Authorship
FuncEGen(filenum).Author = 'Sean Brennan';
FuncEGen(filenum).Date = '2021_07_11';
FuncEGen(filenum).Contact = 'sbrennan@psu.edu';

% Common Flags
FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default



% 
% FuncEGen(filenum).fileNameSuffix = 'polytopesExpandEvenly';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'Expands an obstacle out by exp_dist on all sides.';
% FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
% FuncEGen(filenum).filename_main = 'codeCore_fcn_polytope_editing_expand_polytopes_evenly.m';
% FuncEGen(filenum).filename_script = 'codeCore_script_polytope_editing_expand_polytopes_evenly.m';
% 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'polytopes';
% FuncEGen(filenum).Inputs(1).Type = 'polytopes'; % of polytopes type
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'the structure of ''polytopes'' type that stores the polytopes to be expanded'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'delta';
% FuncEGen(filenum).Inputs(2).Type = 'column_of_numbers'; % Forces it to have only one column. See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Inputs(2).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
% FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'a small number relative to vehicle size to determine the inside of an obstacle'; 
% 
% FuncEGen(filenum).Inputs(3).Name = 'exp_dist';
% FuncEGen(filenum).Inputs(3).Type = 'column_of_numbers'; % Blank means it will not be checked
% FuncEGen(filenum).Inputs(3).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
% FuncEGen(filenum).Inputs(3).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(3).Description = 'distance to expand the obstacle'; 
% 
% FuncEGen(filenum).Inputs(4).Name = 'fig_num';
% FuncEGen(filenum).Inputs(4).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(4).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(4).Description = 'any number that acts somewhat like a figure number output. If given, this forces the variable types to be displayed as output and as well makes the input check process verbose.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 3;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_polytopeCentroidAndArea';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'exp_polytopes';
% FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'structure of expanded polytopes'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Seth Tau';
% FuncEGen(filenum).Date = '2018_11_17, Adjusted example code on 2021_04_28 by Seth Tau, Rebased on 2021_06_26 by S. Brennan';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu and sat5340@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% 
% %% Create generateOneRandomPolytope file
% filenum = 2;
% 
% 
% % Create a polytopesExpandEvenly file
% FuncEGen(filenum).fileNameSuffix = 'generateOneRandomPolytope';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'Generates a single random polytope by selecting one randomly from the Halton Set Voronoi method (HSV) tiling.';
% FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
% FuncEGen(filenum).filename_main = 'codeCore_generateOneRandomPolytope.m';
% FuncEGen(filenum).filename_script = ''; % An empty starter script will be created
% 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'fig_num';
% FuncEGen(filenum).Inputs(1).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(1).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'any number that acts somewhat like a figure number output. If given, this forces the variable types to be displayed as output and as well makes the input check process verbose.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 0;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_haltonVoronoiTiling';
% FuncEGen(filenum).Dependencies(2).Name = 'fcn_MapGen_polytopeCropEdges';
% FuncEGen(filenum).Dependencies(3).Name = 'fcn_MapGen_plotPolytopes';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'one_polytope';
% FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'one randomly generated polytope'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_06_27';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% 
% 
% %% Create snapToAABB file
% filenum = 3;
% 
% 
% % Create a snapToAABB file
% FuncEGen(filenum).fileNameSuffix = 'snapToAABB';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'Given an axis-aligned bounding box (AABB), and a test point, returns a snap point representing the contact point on the closest wall to the test point.';
% FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
% FuncEGen(filenum).filename_main = 'codeCore_snapToAABB.m';
% FuncEGen(filenum).filename_script = 'codeCore_snapToAABB_script.m'; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'axis_aligned_bounding_box';
% FuncEGen(filenum).Inputs(1).Type = '4column_of_numbers'; % 4 columns of numbers
% FuncEGen(filenum).Inputs(1).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'the axis-aligned bounding box, in format [xmin ymin xmax ymax]'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'test_point';
% FuncEGen(filenum).Inputs(2).Type = '2column_of_numbers'; % 2 columns of numbers
% FuncEGen(filenum).Inputs(2).TypeOptions = 1; % Forces it to have only one row. See fcn_MapGen_checkInputsToFunctions  
% FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'the test point, in format [x y]'; 
% 
% FuncEGen(filenum).Inputs(3).Name = 'fig_num';
% FuncEGen(filenum).Inputs(3).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(3).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(3).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 2;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = '(none)';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'snap_point';
% FuncEGen(filenum).Outputs(1).Type = '2column_of_numbers'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'the resulting snap point, in format [x y]'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_02';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% 
% %% Create fillPolytopeFieldsFromVerticies file
% 
% filenum = 4;
% 
% 
% % Create a snapToAABB file
% FuncEGen(filenum).fileNameSuffix = 'fillPolytopeFieldsFromVerticies';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'Given a polytoope structure array where the verticies field is filled, calculates the values for all the other fields.';
% FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
% FuncEGen(filenum).filename_main = 'codeCore_fillPolytopeFieldsFromVerticies.m';
% FuncEGen(filenum).filename_script = 'codeCore_fillPolytopeFieldsFromVerticies_script.m'; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'polytopes';
% FuncEGen(filenum).Inputs(1).Type = 'polytopes'; % of polytopes type
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'an individual structure or structure array of ''polytopes'' type that stores the polytopes to be filled'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'fig_num';
% FuncEGen(filenum).Inputs(2).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(2).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 1;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_polytopeCentroidAndArea';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'filled_polytopes';
% FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'the polytopes array with all fields completed'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_02';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% %% Create fcn_MapGen_polytopeCentroidAndArea file
% 
% filenum = 5;
% 
% 
% % Create a snapToAABB file
% FuncEGen(filenum).fileNameSuffix = 'polytopeCentroidAndArea';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'calculates the centroid and area of a closed polytope.';
% FuncEGen(filenum).long_description = ''; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
% FuncEGen(filenum).filename_main = 'codeCore_polytopeCentroidAndArea.m';
% FuncEGen(filenum).filename_script = 'codeCore_polytopeCentroidAndArea_script.m'; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'vertices';
% FuncEGen(filenum).Inputs(1).Type = '2column_of_numbers'; % of polytopes type
% FuncEGen(filenum).Inputs(2).TypeOptions = [4 5]; % Forces it to have 4 or more rows
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'the list of verticies used to perform calculation, in format [x y] where x and y are column vectors. X: x coordinates of the polytope (with the same first and last point)  Y: y coordinates of the polytope (with the same first and last point)'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'fig_num';
% FuncEGen(filenum).Inputs(2).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(2).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 1;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_checkInputsToFunctions';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'Centroid';
% FuncEGen(filenum).Outputs(1).Type = '2column_of_numbers'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'the calculated centroid of the polytope, given as [x-coordinate y_coordinate]'; 
% 
% FuncEGen(filenum).Outputs(2).Name = 'Area';
% FuncEGen(filenum).Outputs(2).Type = 'numeric'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(2).Description = 'the unsigned area enclosed by the polytope'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_02';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% 
% %% Create fcn_MapGen_polytopeRemoveTightVerticies file
% 
% filenum = 6;
% 
% 
% % Create a snapToAABB file
% FuncEGen(filenum).fileNameSuffix = 'polytopeRemoveTightVerticies';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'removes verticies of polytopes that are too close to each other, measured by a tolerance';
% FuncEGen(filenum).long_description = 'Sometimes, when shrinking, the new verticies are particularly close to each other to where an edge has a trivial length. To prevent this, we get rid of one of any vertices that are too close to each other. This proximity is set by a user-defined tolerance.'; % ImportLongText('fcn_polytope_editing_expand_polytopes_evenly.txt');
% FuncEGen(filenum).filename_main = 'codeCore_fcn_polytopeRemoveTightVerticies.m';
% FuncEGen(filenum).filename_script = 'codeCore_fcn_polytopeRemoveTightVerticies_script.m'; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'polytopes';
% FuncEGen(filenum).Inputs(1).Type = 'polytopes'; % of polytopes type
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'an individual structure or structure array of ''polytopes'' type that stores the polytopes to be evaluated'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'tolerance';
% FuncEGen(filenum).Inputs(2).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'a numeric value that defines how close points should be to be removed'; 
% 
% FuncEGen(filenum).Inputs(3).Name = 'fig_num';
% FuncEGen(filenum).Inputs(3).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(3).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(3).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 2;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_checkInputsToFunctions';
% FuncEGen(filenum).Dependencies(2).Name = 'fcn_MapGen_fillPolytopeFieldsFromVerticies';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'cleaned_polytope';
% FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'the resulting polytope after close edges are removed.'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_02';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% 
% %% Create fcn_MapGen_generatePolysFromTiling file
% 
% filenum = 7;
% 
% 
% % Create a fcn_MapGen_generatePolysFromTiling file
% FuncEGen(filenum).fileNameSuffix = 'generatePolysFromTiling';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'creates polytopes given seed points, V and C matrices from Voronoi tiling, and stretch matrix';
% FuncEGen(filenum).long_description = '';
% FuncEGen(filenum).filename_main = 'codeCore_generatePolysFromTiling.m';
% FuncEGen(filenum).filename_script = ''; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'seed_points';
% FuncEGen(filenum).Inputs(1).Type = '2column_of_numbers'; % of polytopes type
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'the list of seed points in [x y] format, where x and y are columns'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'V';
% FuncEGen(filenum).Inputs(2).Type = ''; % No check
% FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'the V matrix resulting from Voronoi calculations'; 
% 
% FuncEGen(filenum).Inputs(3).Name = 'C';
% FuncEGen(filenum).Inputs(3).Type = ''; % No check
% FuncEGen(filenum).Inputs(3).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(3).Description = 'the C matrix resulting from Voronoi calculations'; 
% 
% FuncEGen(filenum).Inputs(4).Name = 'stretch';
% FuncEGen(filenum).Inputs(4).Type = '2column_of_numbers'; % of polytopes type
% FuncEGen(filenum).Inputs(4).TypeOptions = [1 1]; % Forces it to have only 1 row
% FuncEGen(filenum).Inputs(4).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(4).Description = 'the list of seed points in [x y] format, where x and y are columns'; 
% 
% FuncEGen(filenum).Inputs(5).Name = 'fig_num';
% FuncEGen(filenum).Inputs(5).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(5).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(5).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 4;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_checkInputsToFunctions';
% FuncEGen(filenum).Dependencies(2).Name = 'fcn_MapGen_cropPolytopeToRange';
% FuncEGen(filenum).Dependencies(3).Name = 'fcn_MapGen_fillPolytopeFieldsFromVerticies';
% 
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'polytopes';
% FuncEGen(filenum).Outputs(1).Type = 'polytopes'; % See fcn_MapGen_checkInputsToFunctions
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'the resulting polytopes after converting to polytope form.'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_02';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% %% Create fcn_MapGen_generatePolysFromTiling file
% 
% filenum = 8;
% 
% 
% % Create a fcn_MapGen_generatePolysFromTiling file
% FuncEGen(filenum).fileNameSuffix = 'ugvSensorError';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'calculates error in sensed locations from a UGV perspective';
% FuncEGen(filenum).long_description = ImportLongText('codeCore_err_ugv_v3_header.m');
% FuncEGen(filenum).filename_main = 'codeCore_err_ugv_v3.m';
% FuncEGen(filenum).filename_script = 'codeCore_err_ugv_v3_script.m'; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'Scanning_Results';
% FuncEGen(filenum).Inputs(1).Type = 'positive_3column_of_numbers'; % of polytopes type
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'a Nx3 matrix, N being the number of received lidar points or sensor scans'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'Position_Uncertainty';
% FuncEGen(filenum).Inputs(2).Type = 'positive_3column_of_numbers'; % No check
% FuncEGen(filenum).Inputs(2).TypeOptions = [1 1]; % Forces it to have only 1 row
% FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'a 1x3 vector of constants defining uncertainty in the x, y, and z directions'; 
% 
% FuncEGen(filenum).Inputs(3).Name = 'Angular_Uncertainty';
% FuncEGen(filenum).Inputs(3).Type = 'positive_3column_of_numbers'; % No check
% FuncEGen(filenum).Inputs(3).TypeOptions = [1 1]; % Forces it to have only 1 row
% FuncEGen(filenum).Inputs(3).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(3).Description = 'a 1x3 vector of constants defining uncertainty in the x, y, and z pointing angle (in degrees)'; 
% 
% FuncEGen(filenum).Inputs(4).Name = 'Laser_Uncertainty';
% FuncEGen(filenum).Inputs(4).Type = 'positive_2column_of_numbers'; % No check
% FuncEGen(filenum).Inputs(4).TypeOptions = [1 1]; % Forces it to have only 1 row
% FuncEGen(filenum).Inputs(4).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(4).Description = 'a 1x2 vector of constants defining uncertainty in the LIDAR (see long description)'; 
% 
% FuncEGen(filenum).Inputs(5).Name = 'fig_num';
% FuncEGen(filenum).Inputs(5).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(5).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(5).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 4;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_checkInputsToFunctions';
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'DX_err';
% FuncEGen(filenum).Outputs(1).Type = ''; 
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'Error in the x position'; 
% 
% FuncEGen(filenum).Outputs(2).Name = 'DY_err';
% FuncEGen(filenum).Outputs(2).Type = ''; 
% FuncEGen(filenum).Outputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(2).Description = 'Error in the y position'; 
% 
% FuncEGen(filenum).Outputs(3).Name = 'DZ_err';
% FuncEGen(filenum).Outputs(3).Type = ''; 
% FuncEGen(filenum).Outputs(3).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(3).Description = 'Error in the z position'; 
% 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_07';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default
% 
% %% Create fcn_MapGen_generatePolysFromTiling file
% 
% filenum = 9;
% 
% 
% % Create a fcn_MapGen_ugvSensorErrorBubble file
% FuncEGen(filenum).fileNameSuffix = 'ugvSensorErrorBubble';
% FuncEGen(filenum).class = 'MapGen';
% FuncEGen(filenum).short_description = 'calculates error in sensed locations from a UGV perspective';
% FuncEGen(filenum).long_description = '';
% FuncEGen(filenum).filename_main = 'codeCore_err_ugv_bubble_v3.m';
% FuncEGen(filenum).filename_script = 'codeCore_err_ugv_v3_script.m'; 
% 
% % Input Arguments
% FuncEGen(filenum).Inputs(1).Name = 'Polytopes';
% FuncEGen(filenum).Inputs(1).Type = 'polytopes'; % of polytopes type
% FuncEGen(filenum).Inputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(1).Description = 'an individual structure or structure array of ''polytopes'' type that stores the polytopes to be evaluated'; 
% 
% FuncEGen(filenum).Inputs(2).Name = 'Heading_Angle'; 
% FuncEGen(filenum).Inputs(2).Type = '1column_of_numbers'; 
% FuncEGen(filenum).Inputs(2).TypeOptions = [1 1]; % Forces it to have only 1 row
% FuncEGen(filenum).Inputs(2).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(2).Description = 'is a singular value in degrees of the robot''s heading angle'; 
% 
% FuncEGen(filenum).Inputs(3).Name = 'Bubble_Resolution';
% FuncEGen(filenum).Inputs(3).Type = '1column_of_numbers'; % No check
% FuncEGen(filenum).Inputs(3).TypeOptions = [1 1]; % Forces it to have only 1 row
% FuncEGen(filenum).Inputs(3).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(3).Description = 'is a singular value (unitless)'; 
% 
% FuncEGen(filenum).Inputs(4).Name = 'fig_num';
% FuncEGen(filenum).Inputs(4).Type = 'numeric'; % Must be a number
% FuncEGen(filenum).Inputs(4).Required = 0; % 1 is required, 0 is not
% FuncEGen(filenum).Inputs(4).Description = 'any number that acts as a figure number output, causing a figure to be drawn showing results.'; 
% 
% FuncEGen(filenum).N_RequiredInputs = 3;
% FuncEGen(filenum).N_OptionalInputs = 1;
% 
% % Dependencies
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_checkInputsToFunctions';
% FuncEGen(filenum).Dependencies(1).Name = 'fcn_MapGen_ugvSensorError';
% 
% 
% % Output Arguments
% FuncEGen(filenum).Outputs(1).Name = 'err';
% FuncEGen(filenum).Outputs(1).Type = ''; 
% FuncEGen(filenum).Outputs(1).Required = 1; % 1 is required, 0 is not
% FuncEGen(filenum).Outputs(1).Description = 'Error result'; 
% 
% % Authorship
% FuncEGen(filenum).Author = 'Sean Brennan';
% FuncEGen(filenum).Date = '2021_07_08';
% FuncEGen(filenum).Contact = 'sbrennan@psu.edu';
% 
% % Common Flags
% FuncEGen(filenum).Flags.DoDebug = 1; % Adds the debug flag, header/tailer options
% FuncEGen(filenum).Flags.DoPlot = 1; % Adds the plotting
% FuncEGen(filenum).Flags.CheckInputs = 1; % Adds input checking, turns this on by default


%% Call the function to build the files
fcn_FuncEGen_buildFunctionAndScriptFromDescription(FuncEGen);





function C = ImportLongText(FileName)
S = fileread(FileName);
S(S == char(13)) = [];    % DOS linebreaks to unix linebreaks
S(S == newline) = ' ';   % Replace linebreaks by space - NOTE: this is char(10)
S = strrep(S, '  ', ' '); % Replace double spaces by single space
C = strsplit(S, '§');
end