% script_test_fcn_MapGen_removeInfiniteVertices
% Tests: fcn_MapGen_removeInfiniteVertices

% 
% REVISION HISTORY:
% 
% 2021_07_02 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§


%% pull halton set
halton_points = haltonset(2);
points_scrambled = scramble(halton_points,'RR2'); % scramble values

%% pick values from halton set
Halton_range = [3901        4001];
low_pt = Halton_range(1,1);
high_pt = Halton_range(1,2);
seed_points = points_scrambled(low_pt:high_pt,:);
[V,C] = voronoin(seed_points);
% V = V.*stretch;

%% fill polytopes from tiling
fig_num = 1;
AABB = [0 0 1 1];
stretch = [1 1];


num_poly = size(seed_points,1);
polytopes(num_poly) = ...
    struct(...
    'vertices',[],...
    'xv',[],...
    'yv',[],...
    'distances',[],...
    'mean',[],...
    'area',[],...
    'max_radius',[]);

Npolys = length(polytopes);
Nvertices_per_poly = 20; % Maximum estimate
Nvertices_per_map = Npolys*Nvertices_per_poly;
all_vertices = nan(Nvertices_per_map,3);
all_neighbors = nan(Nvertices_per_map,1);

% Loop through the polytopes, filling all_vertices matrix
for ith_poly = 1:Npolys
    vertices_open = V(C{ith_poly},:); 
    vertices = [vertices_open; vertices_open(1,:)]; % Close off the vertices
    Nvertices = length(vertices(:,1));
    if Nvertices>Nvertices_per_poly
        error('Need to resize the number of allowable vertices');
    else
        row_offset = (ith_poly-1)*Nvertices_per_poly;
        all_vertices(row_offset+1:row_offset+Nvertices,1) = ith_poly;
        all_vertices(row_offset+1:row_offset+Nvertices,2:3) = vertices;
    end
    
       
end


%% Remove infinite vertices
fig_num = 222;
[bounded_vertices] = ...
    fcn_MapGen_removeInfiniteVertices(...
    all_vertices,seed_points,AABB,Nvertices_per_poly,fig_num);

