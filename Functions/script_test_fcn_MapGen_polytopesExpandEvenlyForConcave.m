% script_test_fcn_MapGen_polytopesExpandEvenly
% Tests: fcn_MapGen_polytopesExpandEvenly

%
% REVISION HISTORY:
%§
% 2018_11_17, Seth Tau
% -- first write of script
% 2021_04_28, Seth Tau
% -- Adjusted example code ,
% 2021_06_26 S. Brennan
% -- Rebased code

% Prep the workspace
close all;
clear polytopes;
polytopes = fcn_MapGen_generateOneRandomPolytope;

% xv = [-2 -1 1 2 2 1 -1 -2];
% yv = [-1 -2 -2 -1 1 2 2 1];
% polytopes.vertices = [[xv xv(1)]' [yv yv(1)]'];
% polytopes.xv = xv;
% polytopes.yv = yv;
%
% polytopes.distances = sum((polytopes(1).vertices(1:end-1,:)-polytopes(1).vertices(2:end,:)).^2,2).^0.5;
% [Cx,Cy,polytopes.area] = fcn_MapGen_polytopeCentroidAndArea([xv xv(1)],[yv yv(1)]);
% polytopes.mean = [Cx, Cy];
% polytopes.max_radius = max(sum((polytopes.vertices(1:end-1,:)-ones(length(xv),1)*polytopes.mean).^2,2).^0.5);


%% concave polytope
polytopes.vertices = [
    1.0000    0.5217
    1.0000    0.5242
    0.9300    0.6329
    0.8472    0.6479
    0.85    0.52
    0.9    0.58
    0.9    0.62
    1.0000    0.5217
];

polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);

% Set parameters
% delta = 0.01; % Set the delta value (what is this used for?)
exp_dist = 0.04; % Set the expansion distance
fig_num = 221; % Set the figure number

% notice how the concavity is not well represented when enlarging the polytope
exp_polytopes=fcn_MapGen_polytopesExpandEvenly(polytopes,exp_dist,fig_num);
hold on; box on;
scale_factor = 1.2;
for p = 1:length(polytopes)
    this_polytope = polytopes(p);
    this_polyshape = polyshape(this_polytope.vertices);
    scaled_polyshape = scale(this_polyshape,scale_factor,this_polytope.mean);
    plot(scaled_polyshape);
end

hold on; box on;
scale_factor = 1.2;
clear new_polytopes;
for p = 1:length(polytopes)
    this_polytope = polytopes(p);
    this_polyshape = polyshape(this_polytope.vertices);
    scaled_polyshape = polybuffer(this_polyshape,exp_dist,'JointType','miter','MiterLimit',2);
    plot(scaled_polyshape);
    new_vertices = scaled_polyshape.Vertices;
    new_vertices = [new_vertices; new_vertices(1,:)];
    new_polytopes(p).vertices = new_vertices;
end
new_polytopes= fcn_MapGen_fillPolytopeFieldsFromVertices(new_polytopes);
fig_num = fig_num + 1;
line_width = 2;
fcn_MapGen_plotPolytopes(new_polytopes,fig_num,'r-',line_width);

