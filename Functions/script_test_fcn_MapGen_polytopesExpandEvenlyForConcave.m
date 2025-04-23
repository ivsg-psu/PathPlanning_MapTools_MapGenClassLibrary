% script_test_fcn_MapGen_polytopesExpandEvenlyForConcave
% Tests: fcn_MapGen_polytopesExpandEvenlyForConcave

%
% REVISION HISTORY:
%§
% 2024_02_05, Steve Harnett
% -- first write of script

% Prep the workspace
close all;
clear polytopes;
polytopes = fcn_MapGen_generateOneRandomPolytope;


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

polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes,0,1);
fcn_MapGen_plotPolytopes(polytopes,1,'r',2);
% Set parameters
exp_dist = 0.04; % Set the expansion distance
fig_num = 221; % Set the figure number

%% use legacy polytope expansion
% notice how the concavity is not well represented when enlarging the polytope with *ExpandEvenly
% this code is taken directly from fcn_MapGen_polytopesExpandEvenly
% code is reproduced here because this function requires slight modification to work with nonconvex
% obstacles and this modification should not be made to the actual function

exp_polytopes = polytopes; % both structures will be the same size

for ith_poly = 1:size(polytopes,2) % check each obstacle

    % pull values
    vertices = polytopes(ith_poly).vertices;
    centroid = polytopes(ith_poly).mean;
    rad = polytopes(ith_poly).max_radius;

    % Calculate scale
    my_scale = (rad+exp_dist)/rad;

    % Calculate new vertices
    exp_polytopes(ith_poly).vertices = centroid + my_scale*(vertices-centroid);

    % fill in other fields from the vertices field
    exp_polytopes(ith_poly) = fcn_MapGen_fillPolytopeFieldsFromVertices(exp_polytopes(ith_poly),1,1);

end
hold on; box on;

%% use scale method
% notice how the concavity is not well represented when enlarging the polytope with the scale polyshape method
scale_factor = 1.2;
for p = 1:length(polytopes)
    this_polytope = polytopes(p);
    this_polyshape = polyshape(this_polytope.vertices);
    scaled_polyshape = scale(this_polyshape,scale_factor,this_polytope.mean);
    plot(scaled_polyshape);
end

%% use concave expansion function
exp_polytopes=fcn_MapGen_polytopesExpandEvenlyForConcave(polytopes,exp_dist);
line_width =2;
fig_num = 1;
fcn_MapGen_plotPolytopes(exp_polytopes,fig_num,'g-',line_width);
legend('original','*ExpandEvenly','','scale method','*ExpandEvenlyForConvex');
