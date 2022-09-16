% script_ui_manuallyDefineMapLayers
%
% allow the user to draw polytopes on an opened figure
% store these polytopes in a 2D polytope matrix where each row is a map layer
% (i.e. a set of polytopes defining one type of feature like impassable obstacles
% water features, passable obstacles, rough patches, etc.) and each column is
% a single polytope
%
% REVISION HISTORY:
%
% 2022_15_Sept by S. Harnett
% -- first write of script
% 2022_16_Sept by S. Harnett
% -- allowing user to set polytope costs and changing from matrices of polytopes to cell array
% TODO(@sjharnett) enforce convexity or break concave polys into convex polys
%%%%%%%%%%%%%%§

clear all; close all; clc;
% prompt user to open a fig file
sprintf("Select a satellite photo or other map figure with appropriate coordinates\n")
[FileName,PathName,FilterIndex] = uigetfile('.fig')
% open the file of interest, hold on so we can draw on it
openfig(strcat(PathName,FileName))
hold on;

% want to create a 2D polytope array, of different "maps layers"
prompt_map_layers = "How many map layers would you like to make?\n";
map_layers = input(prompt_map_layers)
% each map layer should have its own color, if possible
colors = ['r','g','c','m','b','y'];
if map_layers > length(colors)
    warning("There are more map layers than colors to plot them in, colors will be recycled\n")
end

for map_layer = 1:map_layers
    prompt_polys = sprintf("How many polytopes would you like to make on layer %i of %i?\n",map_layer,map_layers);
    polys = input(prompt_polys)
    polytopes_this_layer = [];
    for poly = 1:polys
        sprintf("Click all the vertices you would like to make in polytope %i of %i,\n then press return. You do not need to close the polytope (this is done automatically \n between the first and last points.\n",poly,polys)
        pos = ginput

        % close the polytope by copying first vertex to end
        pos = [pos; pos(1,:)]

        % add vertices to polytope array
        polytopes_this_layer(poly).vertices = pos;

        % pull out x and y for plotting
        x = pos(:,1);
        y = pos(:,2);

        % want to plot all polys from each map layer in the same color, as they're drawn
        % check if we're out of colors, if so, use modulo of map layer to loop back through color array
        if map_layer <= length(colors)
            plot(x,y,'LineWidth',2,'Color',colors(map_layer))
        else
            plot(x,y,'LineWidth',2,'Color',colors(mod(map_layer,length(colors))))
        end
    end
    polytopes{map_layer} = polytopes_this_layer;
end
% fill out polytope struct from vertices
for i = 1:length(polytopes)
    prompt_polytope_cost = sprintf("What traversal cost would you like for polytopes in layer %i?\n",i);
    des_polytope_cost = input(prompt_polytope_cost)
    polytopes{i} = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes{i});
    polytopes{i} = fcn_polytope_editing_set_all_costs(polytopes{i},des_polytope_cost);
end

