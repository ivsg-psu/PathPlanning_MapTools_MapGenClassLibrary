% script_test_fcn_MapGen_flattenPolytopeMap
% Tests: fcn_MapGen_polytopesStatistics

%
% REVISION HISTORY:
%
% 2022_10_19 by Steve Harnett
% -- first write of script
% 2022_10_21 by Steve Harnett
% -- test fixtures and assertions added
% 2025_04_25 by Sean Brennan
% -- added test fixtures to global path during main code, removed relative
% dependencies here which was breaking the code on different OS versions
%%%%%%%%%%%%%%§

close all; 

%% tests two sets of overlapping polytopes
load('overlapped_real_world_polys'); % located in testFixtures subdirectory
flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes);
% this will produce 29 polys, 6 of which are very small so we assert that
% these are removed
assert(isequal(length(flattened_polytopes),23));
% assert that no polytope has an area smaller than 1e-10
assert(isequal(sum(extractfield(flattened_polytopes,'area')<1e-10),0));

close all; 

%% tests an large polytope with several enclave polytopes
% this is how free space would be broken into polytopes
load('bound_entire_map'); % located in testFixtures subdirectory
% flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes);
flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes(1:2));
assert(isequal(length(flattened_polytopes),12));
% assert that no polytope has an area smaller than 1e-10
assert(isequal(sum(extractfield(flattened_polytopes,'area')<1e-10),0));
