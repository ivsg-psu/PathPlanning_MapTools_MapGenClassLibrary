% script_test_fcn_MapGen_increasePolytopeVertexCount
% Tests: fcn_MapGen_increasePolytopeVertexCount

%
% REVISION HISTORY:
%
% 2022_10_21 by Steve Harnett
% -- first write of script
%%%%%%%%%%%%%%§

close all; clear all; clc;
%% tests two sets of overlapping polytopes
axis_box = [0 1 0 1];
halton_range = [1 100]; % Sets the indices to use from halton set
[polytopes] = fcn_MapGen_haltonVoronoiTiling(halton_range);
fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],axis_box,'square',[1 0 0 0 0.5]);
resolution = 0.05;
interpolated_polytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes,resolution/2,100)
fcn_MapGen_plotPolytopes(interpolated_polytopes,1000,'-',2,[0 0 0],axis_box,'square',[1 0 0 0 0.5]);
% assert that there are more vertices now than there were before
assert(size(extractfield(polytopes,'vertices'),2)<size(extractfield(interpolated_polytopes,'vertices'),2));
