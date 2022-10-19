% script_test_fcn_MapGen_flattenPolytopeMap
% Tests: fcn_MapGen_polytopesStatistics

%
% REVISION HISTORY:
%
% 2022_10_19 by Steve Harnett
% -- first write of script
%%%%%%%%%%%%%%§

close all; clear all; clc;
% load('Figures/overlapped_real_world_polys')
% flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes_flattened)
% close all; clear all; clc;
load('Figures/bound_entire_map')
flattened_polytopes = fcn_MapGen_flattenPolytopeMap(polytopes_flattened)