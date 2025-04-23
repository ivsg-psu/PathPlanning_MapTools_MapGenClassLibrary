function interpolated_polytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes,resolution)
    % fcn_MapGen_increasePolytopeVertexCount
    % Given polytope field and a desired resolution distance, n, returns an equivalent
    % polytope field with colinear vertices added to each polytope side such that
    % there is a vertex every n units
    % The utility of this is that if path planning is restricted to using polytope
    % vertices as waypoints, this increases the number of options the planner has
    % while keeping the obstacle field the same
    %
    %
    %
    % FORMAT:
    % interpolated_polytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes,resolution)
    %
    % INPUTS:
    %     polytopes - the initial polytope field
    %     resolution - the desired linear spacing between vertices along each polytope side
    %
    % OUTPUTS:
    %
    %
    %     interpolated_polytopes - a polytope field equivalent to the input but with vertices added
    %       along the polytopes sides every RESOLUTION units such that each polytope now has more vertices
    %
    % DEPENDENCIES:
    %
    % EXAMPLES:
    %
    % See the script: script_fcn_MapGen_increasePolytopeVertexCount.m
    % for a full test suite.
    %
    % Questions or comments? contact sjh6473@psu.edu

    % REVISION HISTORY:
    % 2021_10_13
    % -- first written by Steve Harnett


    %% Debugging and Input checks
    flag_do_plot = 0;      % Set equal to 1 for plotting
    flag_do_debug = 0;     % Set equal to 1 for debugging
    if flag_do_debug
        fig_for_debug = 747;
        fig_num = 2100;
        st = dbstack; %#ok<*UNRCH>
        fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    end

    %% Start of main code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   __  __       _
    %  |  \/  |     (_)
    %  | \  / | __ _ _ _ __
    %  | |\/| |/ _` | | '_ \
    %  | |  | | (_| | | | | |
    %  |_|  |_|\__,_|_|_| |_|
    %
    %See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
    % loop through all polytopes
    for i = 1:length(polytopes)
        new_verts = [];
        if flag_do_plot
            figure(102938)
            hold on
            % plot original polytope vertices
            plot(polytopes(i).vertices(:,1),polytopes(i).vertices(:,2),'k')
            plot(polytopes(i).vertices(:,1),polytopes(i).vertices(:,2),'kx')
        end
        % loop through all vertices
        for j = 1:(size(polytopes(i).vertices,1)-1)
            % note the side range in x
            delta_x = polytopes(i).vertices(j,1) - polytopes(i).vertices(j+1,1);
            % note hte side range in y
            delta_y = polytopes(i).vertices(j,2) - polytopes(i).vertices(j+1,2);
            % find the side length
            side_length = sqrt(delta_x^2+delta_y^2);
            % determine number of vertices to place on this side per resolution
            num_verts_needed_this_side = side_length/resolution;
            % evenly space vertices in x since interp1 interpolates in x
            resolution_in_x = abs(delta_x)/num_verts_needed_this_side;
            % we expect the lower index vertex to have a lower x value than the higher
            % index vertex however...
            if polytopes(i).vertices(j,1) < polytopes(i).vertices(j+1,1)
                % sample x points for interp1
                xq = polytopes(i).vertices(j,1):resolution_in_x:polytopes(i).vertices(j+1,1);
                % see matlab documentation on interp1 for further examples of what this does
                vq = interp1([polytopes(i).vertices(j,1),polytopes(i).vertices(j+1,1)],...
                                [polytopes(i).vertices(j,2),polytopes(i).vertices(j+1,2)],...
                                xq);
            % ...if this side goes backwards (from larger x to smaller x) we
            % have to flip the interpolated points, additionally...
            elseif polytopes(i).vertices(j,1) > polytopes(i).vertices(j+1,1)
                xq = polytopes(i).vertices(j+1,1):resolution_in_x:polytopes(i).vertices(j,1);
                vq = interp1([polytopes(i).vertices(j+1,1),polytopes(i).vertices(j,1)],...
                                [polytopes(i).vertices(j+1,2),polytopes(i).vertices(j,2)],...
                                xq);
                % if this is a "backwards side" flip the interpolated vectors back
                vq = flip(vq);
                xq = flip(xq);
            % ...if this side is completely vertical (i.e. the vertices have the same x)...
            elseif polytopes(i).vertices(j,1) == polytopes(i).vertices(j+1,1)
                % ...then we can interpolate in y instead
                % we don't need to use interp1 because we know the slope of the side is inf.
                % so simple linearly space between the low y and high y
                vq = linspace(min(polytopes(i).vertices(j,2),polytopes(i).vertices(j+1,2)),...
                        max(polytopes(i).vertices(j,2),polytopes(i).vertices(j+1,2)),...
                        num_verts_needed_this_side);
                % and make an x vector of the appropriate size containing the only x value possible
                xq = ones(1,length(vq)).*polytopes(i).vertices(j,1);
            end
            % log new vertecies for this side
            new_verts_this_side = [xq; vq]';
            % append to array of new vertices for this polytope
            new_verts = [new_verts; new_verts_this_side];
        end
        % update this polytope's fields based on the new vertices found
        polytopes(i).vertices = [new_verts;new_verts(1,:)];
        polytopes(i).xv = new_verts(:,1)';
        polytopes(i).yv = new_verts(:,2)';
        % plot the new vertices in a different color for comparison
        if flag_do_plot
            figure(102938)
            hold on
            plot(polytopes(i).vertices(:,1),polytopes(i).vertices(:,2),'rx')
            legend('','orig. vertices','dense vertices','','')
        end
    end
    interpolated_polytopes = polytopes;
end
