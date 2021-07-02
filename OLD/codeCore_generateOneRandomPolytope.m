
%% main code §

% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box);

% Pick a random polytope
Npolys = length(trim_polytopes);
rand_poly = 1+floor(rand*Npolys);
one_polytope = trim_polytopes(rand_poly);

% Plot results
fig_num = 2;
fcn_MapGen_plotPolytopes(one_polytope,fig_num,'-',2,[0.5 0 0])

% §
% Debug
%
% Functions §