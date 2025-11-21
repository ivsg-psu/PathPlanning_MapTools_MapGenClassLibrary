% script_demo_removalOfSelfIntersection
% Tests function: fcn_MapGen_polytopeFindSelfIntersections

% REVISION HISTORY:
% 
% 2025_07_16 by Sean Brennan, sbrennan@psu.edu
% - first written by S. Brennan by copying out of
%  script_test_fcn_MapGen_polytopeRemoveColinearVertices
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - Updated rev history to be in Markdown format
% - Replaced fig_+num with figNum

% TO-DO:
% 
% 2025_11_20 by Sean Brennan, sbrennan@psu.edu
% - fill in to-do items here.


close all;


%% Basic example of self-intersection being "cleaned"
figNum = 1;
vertices = [0 0; 1 0; 0.5 1.5; 1 1; 0 1; 0 0];
vertices_with_self_intersects = fcn_MapGen_polytopeFindSelfIntersections(...
    vertices,figNum);

figNum = 11;
interior_point = [0.5 0.5];
[projected_points] = ...
    fcn_MapGen_polytopeProjectVerticesOntoWalls(...,
    interior_point,...
    vertices_with_self_intersects,...
    vertices_with_self_intersects(1:end-1,:),...
    vertices_with_self_intersects(2:end,:),...
    11);

assert(isequal(round(projected_points,4),[0,0; 0,0; 1,0; 0.75,0.75; 0.6667,1; 0.6667,1; 0.5,1; 0,1]));

    figNum = 12;
[cropped_vertices] = ...
    fcn_MapGen_polytopeRemoveColinearVertices(...,
    projected_points,...
    figNum);

assert(isequal(round(cropped_vertices,4),[0,0; 1,0; 0.75,0.75; 0.6667,1; 0,1]));