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
% TODO(@sjharnett) enforce convexity or break concave polys into convex polys
% TODO(@sjharnett) should multiple map layers be stored in a nxm matrix or an array of arrays?
%%%%%%%%%%%%%%§


% prompt user to open a fig file
sprintf("Select a satellite photo or other map figure with appropriate coordinates")
[FileName,PathName,FilterIndex] = uigetfile('.fig')
% open the file of interest, hold on so we can draw on it
openfig(strcat(PathName,FileName))
hold on;

% want to create a 2D polytope array, of different "maps layers"
prompt_map_layers = "How many map layers would you like to make?";
map_layers = input(prompt_map_layers)
% each map layer should have its own color, if possible
colors = ['r','g','c','m','b','y'];
if map_layers > length(colors)
    warning("There are more map layers than colors to plot them in, colors will be recycled")
end

for map_layer = 1:map_layers
    prompt_polys = sprintf("How many polytopes would you like to make on layer %i of %i?",map_layer,map_layers);
    polys = input(prompt_polys)
    for poly = 1:polys
        sprintf("Click all the vertices you would like to make in polytope %i of %i,\n then press return. You do not need to close the polytope (this is done automatically \n between the first and last points.",poly,polys)
        pos = ginput

        % close the polytope by copying first vertex to end
        pos = [pos; pos(1,:)]

        % add vertices to polytope array
        polytopes(map_layer,poly).vertices = pos;

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
end
% fill out polytope struct from vertices
polytopes = fcn_MapGen_fillPolytopeFieldsFromVertices(polytopes);
