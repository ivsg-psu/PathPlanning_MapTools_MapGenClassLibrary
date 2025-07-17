% script_test_fcn_MapGen_polytopesExpandEvenlyForConcave
% Tests: fcn_MapGen_polytopesExpandEvenlyForConcave

%
% REVISION HISTORY:
%§
% 2024_02_05, Steve Harnett
% -- first write of script

% Prep the workspace
close all;

%% Show how different polytope expansion methods do NOT work

%%%%%
% Fill in one polytope and plot it
fig_num = 1;

if 1==1
    rng(1);
    polytopes = fcn_MapGen_polytopeGenerateOneRandomPoly(fig_num);

else
    % concave polytope
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

    polytopes = fcn_MapGen_polytopesFillFieldsFromVertices(polytopes,fig_num);
end
% fcn_MapGen_plotPolytopes(polytopes,1,'r',2);
% Set parameters
exp_dist = 0.2*polytopes.max_radius; % Set the expansion distance to get a 20 percent increase

%%%%%
% What if we enlarge by simply moving vertices away from the centroid?
%
% The following code shows how the concavity is not well represented when
% enlarging the polytope with *ExpandEvenly this code is taken directly
% from fcn_MapGen_polytopesExpandEvenly code is reproduced here because
% this function requires slight modification to work with nonconvex
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
    exp_polytopes(ith_poly) = fcn_MapGen_polytopesFillFieldsFromVertices(exp_polytopes(ith_poly),1,fig_num);

end
hold on; box on;
% The plot should illustrate that some sides get expanded FAR more than
% other sides. In other words, the scaling factor pushes "outward"
% differently on different parts of the polytope. So this is not a very
% good solution

%%%%
% Would a scaling method work better?
% The following example shows how the concavity is not well represented
% when enlarging the polytope with the scale polyshape method
scale_factor = 1.2;
for p = 1:length(polytopes)
    this_polytope = polytopes(p);
    this_polyshape = polyshape(this_polytope.vertices);
    scaled_polyshape = scale(this_polyshape,scale_factor,this_polytope.mean);
    plot(scaled_polyshape);
end

% Notice that this gives (usually) EXACTLY the same answer

%%%%
% What if we use a concave expansion function?
exp_polytopes=fcn_MapGen_polytopesExpandEvenlyForConcave(polytopes,exp_dist, fig_num);

% Notice now how all the walls are pushed out evenly

legend('original','*ExpandEvenly','','scale method','*ExpandEvenlyForConvex');
