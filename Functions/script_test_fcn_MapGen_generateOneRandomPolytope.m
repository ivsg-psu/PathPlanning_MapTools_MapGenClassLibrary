% script_test_fcn_MapGen_generateOneRandomPolytope
% Tests: fcn_MapGen_generateOneRandomPolytope

% 
% REVISION HISTORY:
% 
% 2021_06_27 by Sean Brennan
% -- first write of script
%%%%%%%%%%%%%%§

fig_num = 1;
figure(fig_num);
clf;

one_polytope = fcn_MapGen_generateOneRandomPolytope(fig_num);
assert(isstruct(one_polytope))

