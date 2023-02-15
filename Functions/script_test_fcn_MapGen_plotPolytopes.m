% script_test_fcn_MapGen_plotPolytopes
% Tests function: fcn_MapGen_plotPolytopes

% REVISION HISTORY:
% 2021_06_07
% -- first written by S. Brennan.

%% BASIC example
clear polytopes;

polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
fig_num = 1;
line_spec = '-';
line_width = 3;
line_color = [0.5 0 0];

fcn_MapGen_plotPolytopes(polytopes,fig_num,line_spec,line_width,line_color);

%% Advanced example
axis_box = [0 1 0 1];
halton_range = [1 100]; % Sets the indices to use from halton set
[polytopes] = fcn_MapGen_haltonVoronoiTiling(halton_range);
fig1 = fcn_MapGen_plotPolytopes(polytopes,[],'-',2,[0.5 0 0]);
fig2 = fcn_MapGen_plotPolytopes(polytopes,998,'-',2,[0 0 0.5],axis_box);
fig3 = fcn_MapGen_plotPolytopes(polytopes,999,'-',2,[0 0.5 0],axis_box,'square');
fig4 = fcn_MapGen_plotPolytopes(polytopes,1000,'-',2,[0 0 0],axis_box,'square',[1 0 0 0 0.5]);

assert(true);

% To show that overplotting works, we redo plotting but with a different
% set of polytopes, on the same figures as before
halton_range = [101 200]; % Sets the indices to use from halton set
[polytopes2] = fcn_MapGen_haltonVoronoiTiling(halton_range);
fcn_MapGen_plotPolytopes(polytopes2,fig1,'r-',2)
fcn_MapGen_plotPolytopes(polytopes2,fig2,'b--',2,axis_box)
fcn_MapGen_plotPolytopes(polytopes2,fig3,'g-',3,axis_box,'square')
fcn_MapGen_plotPolytopes(polytopes2,fig4,'k-',3,axis_box,'square',[1 0 0 0 0.5])
