% script_test_fcn_MapGen_polytopeShrinkEvenly
% Tests: fcn_MapGen_polytopeShrinkEvenly

% 
% REVISION HISTORY:
% 
% 2025_07_29 - S. Brennan, sbrennan@psu.edu
% -- first write of script 
%    % * using fcn_MapGen_polytopeGenerateOneRandomPoly as starter

%% Set up the workspace
close all

%% Code demos start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                              ____   __    _____          _
%  |  __ \                            / __ \ / _|  / ____|        | |
%  | |  | | ___ _ __ ___   ___  ___  | |  | | |_  | |     ___   __| | ___
%  | |  | |/ _ \ '_ ` _ \ / _ \/ __| | |  | |  _| | |    / _ \ / _` |/ _ \
%  | |__| |  __/ | | | | | (_) \__ \ | |__| | |   | |___| (_) | (_| |  __/
%  |_____/ \___|_| |_| |_|\___/|___/  \____/|_|    \_____\___/ \__,_|\___|
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Demos%20Of%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: basic demo of fcn_MapGen_polytopeShrinkEvenly
fig_num = 10001;
titleString = sprintf('DEMO case: basic demo of fcn_MapGen_polytopeShrinkEvenly');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 0.02;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: repeated cutting
fig_num = 10002;
titleString = sprintf('DEMO case: repeated cutting');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [0 0; 0.4 0.1; 1 1; 0 1; 0 0]*5;
unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

step =  0.1;
for edgeCut = step:step:1

    % Call the function
    [shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
        fcn_MapGen_polytopeShrinkEvenly(...
        unshrunkPolytope,...
        edgeCut,...
        (fig_num));
    legend('off');
end

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% DEMO case: repeated cutting of random poly
fig_num = 10003;
titleString = sprintf('DEMO case: repeated cutting of random poly');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set up polytopes
seedGeneratorNames = 'haltonset';
seedGeneratorRanges = [1 100];
AABBs = [0 0 1 1];
mapStretchs = [1 1];
[polytopes] = fcn_MapGen_generatePolysFromSeedGeneratorNames(...
    seedGeneratorNames,...  % string or cellArrayOf_strings with the name of the seed generator to use
    seedGeneratorRanges,... % vector or cellArrayOf_vectors with the range of points from generator to use
    (AABBs),...             % vector or cellArrayOf_vectors with the axis-aligned bounding box for each generator to use
    (mapStretchs),...       % vector or cellArrayOf_vectors to specify how to stretch X and Y axis for each set
    (-1));


bounding_box = [0,0, 1,1];
trim_polytopes = fcn_MapGen_polytopesDeleteByAABB(polytopes,bounding_box);

% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
unshrunkPolytope = trim_polytopes(rand_poly);

edge_cut_step = 0.002;
for edgeCut = edge_cut_step:edge_cut_step:(unshrunkPolytope.max_radius/1.5+edge_cut_step)

    % Call the function
    [shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
        fcn_MapGen_polytopeShrinkEvenly(...
        unshrunkPolytope,...
        edgeCut,...
        (fig_num));
    legend('off')

end

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(all(isnan(shrunkPolytope.vertices),'all')); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Test cases start here. These are very simple, usually trivial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _______ ______  _____ _______ _____
% |__   __|  ____|/ ____|__   __/ ____|
%    | |  | |__  | (___    | | | (___
%    | |  |  __|  \___ \   | |  \___ \
%    | |  | |____ ____) |  | |  ____) |
%    |_|  |______|_____/   |_| |_____/
%
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases\n');

%% TEST case: vertical wall polytope
fig_num = 20001;
titleString = sprintf('TEST case: vertical wall polytope');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [
    0 0; 2 0; 1 2; 0 1; 0 0]*5;
unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 0.1;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 
% extract new vertices
new_vertices = shrunkPolytope.vertices;
% assert that new vertices are within 5% error of having the same x
% position
error_tolerance = 1E-6;
vertical_error = abs(new_vertices(1,1)-new_vertices(4,1))/new_vertices(4,1);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
    ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
    new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: make sure not to create concave polytope
fig_num = 20002;
titleString = sprintf('TEST case: make sure not to create concave polytope');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;

unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 3.6;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% assert that the polytope is convex to start
[angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(unshrunkPolytope.vertices);
interior_angles = 180-angles*180/pi;
assert(~any(interior_angles>180));

% assert that the polytope is convex after shrinking
new_vertices = shrunkPolytope.vertices;
[new_angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(new_vertices);
new_interior_angles = 180-new_angles*180/pi;
assert(~any(new_interior_angles>180),['All interior angles must be < 180 ',...
    'polytope to be convex']);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: NaN produced if overcut
fig_num = 20003;
titleString = sprintf('TEST case: NaN produced if overcut');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;

unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 4;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(all(isnan(shrunkPolytope.vertices),'all')); 

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: a square
fig_num = 20004;
titleString = sprintf('TEST case: a square');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*5;

unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 1;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% assert that the polytope is convex to start
[angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(unshrunkPolytope.vertices);
interior_angles = 180-angles*180/pi;
assert(~any(interior_angles>180));

% assert that the polytope is convex after shrinking
new_vertices = shrunkPolytope.vertices;
[new_angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(new_vertices);
new_interior_angles = 180-new_angles*180/pi;
assert(~any(new_interior_angles>180),['All interior angles must be < 180 ',...
    'polytope to be convex']);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% TEST case: a triangle
fig_num = 20005;
titleString = sprintf('TEST case: a triangle');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [0 0; 1 1; 0 1; 0 0]*5;

unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 1;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% assert that the polytope is convex to start
[angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(unshrunkPolytope.vertices);
interior_angles = 180-angles*180/pi;
assert(~any(interior_angles>180));

% assert that the polytope is convex after shrinking
new_vertices = shrunkPolytope.vertices;
[new_angles, ~, ~] = fcn_MapGen_polytopeFindVertexAngles(new_vertices);
new_interior_angles = 180-new_angles*180/pi;
assert(~any(new_interior_angles>180),['All interior angles must be < 180 ',...
    'polytope to be convex']);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 0.02;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    ([]));

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 0.02;

% Call the function
[shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
    fcn_MapGen_polytopeShrinkEvenly(...
    unshrunkPolytope,...
    edgeCut,...
    (-1));

% Check variable types
assert(isstruct(shrunkPolytope));
assert(isfield(shrunkPolytope,'vertices'));
assert(isfield(shrunkPolytope,'xv'));
assert(isfield(shrunkPolytope,'yv'));
assert(isfield(shrunkPolytope,'distances'));
assert(isfield(shrunkPolytope,'mean'));
assert(isfield(shrunkPolytope,'area'));
assert(isfield(shrunkPolytope,'max_radius'));
assert(isfield(shrunkPolytope,'min_radius'));
assert(isfield(shrunkPolytope,'mean_radius'));
assert(isfield(shrunkPolytope,'radii'));
assert(isfield(shrunkPolytope,'cost'));
assert(isfield(shrunkPolytope,'parent_poly_id'));
assert(iscell(newVertices));
assert(iscell(newProjectionVectors));
assert(isnumeric(cutDistance));

% Check variable sizes
assert(isequal(length(unshrunkPolytope),length(shrunkPolytope))); 
Ncuts = length(cutDistance(:,1));
assert(Ncuts==length(newVertices));
assert(Ncuts==length(newProjectionVectors));

% Check variable values
assert(shrunkPolytope.area<unshrunkPolytope.area); 
assert(all(shrunkPolytope.max_radius<unshrunkPolytope.max_radius)); 
assert(all(shrunkPolytope.mean_radius<unshrunkPolytope.mean_radius)); 

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Set parameters
unshrunkPolytopeVerticesOnly = fcn_MapGen_polytopeFillEmptyPoly(-1);
unshrunkPolytopeVerticesOnly.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.8921    0.5627
    1.0000    0.5217
];
unshrunkPolytope = fcn_MapGen_polytopesFillFieldsFromVertices(unshrunkPolytopeVerticesOnly,([]),(-1));

edgeCut = 0.02;

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
        fcn_MapGen_polytopeShrinkEvenly(...
        unshrunkPolytope,...
        edgeCut,...
        ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
        fcn_MapGen_polytopeShrinkEvenly(...
        unshrunkPolytope,...
        edgeCut,...
        (-1));
end
fast_method = toc;

% Do calculation with pre-calculation, FAST_MODE on, with precalculation
precalcVertices = newVertices;
precalcProjectionVectors = newProjectionVectors;
precalcCutDistance = cutDistance;
tic;
for ith_test = 1:Niterations
    % Call the function
    [shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
        fcn_MapGen_polytopeShrinkEvenly(...
        unshrunkPolytope,...
        edgeCut,...
        precalcVertices, precalcProjectionVectors, precalcCutDistance,...
        (-1));
end
fast_fast_method = toc;

% Do calculation with pre-calculation, FAST_MODE on, with precalculation,
% fast function
precalcVertices = newVertices;
precalcProjectionVectors = newProjectionVectors;
precalcCutDistance = cutDistance;
tic;
for ith_test = 1:Niterations
    % Call the function
    [shrunkPolytope, newVertices, newProjectionVectors, cutDistance]= ...
        fcn_MapGen_polytopeShrinkEvenly_FAST(...
        unshrunkPolytope,...
        edgeCut,...
        precalcVertices, precalcProjectionVectors, precalcCutDistance,...
        (-1));
end
fast_fast_fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

categoryStrings = {'Normal mode','Fast mode','Fast mode with precalculation','Fast mode, precalculated, fast fcn version'};
X = categorical(categoryStrings);
X = reordercats(X,categoryStrings); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method fast_fast_method fast_fast_fast_method]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG

%% Fail conditions
if 1==0
    %

end


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


% %% fcn_INTERNAL_loadExampleData
% function [seed_points, V, C] = fcn_INTERNAL_loadExampleData
%
%
% % pull halton set
% halton_points = haltonset(2);
% points_scrambled = scramble(halton_points,'RR2'); % scramble values
%
% % pick values from halton set
% Halton_range = [1801 1901];
% low_pt = Halton_range(1,1);
% high_pt = Halton_range(1,2);
% seed_points = points_scrambled(low_pt:high_pt,:);
% [V,C] = voronoin(seed_points);
% % V = V.*stretch;
% end % Ends fcn_INTERNAL_loadExampleData