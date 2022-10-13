function interpolated_polytopes = fcn_MapGen_increasePolytopeVertexCount(polytopes,resolution)
    resolution = 0.005; % TODO something better than a hard coded resolution could be used here
    flag_do_plot = 1;
    for i = 1:length(shrunk_polytopes)
        new_verts = [];
        if flag_do_plot
            figure(102938)
            hold on
            plot(shrunk_polytopes(i).vertices(:,1),shrunk_polytopes(i).vertices(:,2),'k')
            plot(shrunk_polytopes(i).vertices(:,1),shrunk_polytopes(i).vertices(:,2),'kx')
        end
        for j = 1:(size(shrunk_polytopes(i).vertices,1)-1)
            % if this side goes backwards (from larger x to smaller x) we
            % have to flip the interpolated points
            delta_x = shrunk_polytopes(i).vertices(j,1) - shrunk_polytopes(i).vertices(j+1,1);
            delta_y = shrunk_polytopes(i).vertices(j,2) - shrunk_polytopes(i).vertices(j+1,2);
            side_length = sqrt(delta_x^2+delta_y^2);
            num_verts_needed_this_side = side_length/resolution;
            resolution_in_x = abs(delta_x)/num_verts_needed_this_side;
            if shrunk_polytopes(i).vertices(j,1) < shrunk_polytopes(i).vertices(j+1,1)
                xq = shrunk_polytopes(i).vertices(j,1):resolution_in_x:shrunk_polytopes(i).vertices(j+1,1);
                vq = interp1([shrunk_polytopes(i).vertices(j,1),shrunk_polytopes(i).vertices(j+1,1)],...
                                [shrunk_polytopes(i).vertices(j,2),shrunk_polytopes(i).vertices(j+1,2)],...
                                xq);
            else
                xq = shrunk_polytopes(i).vertices(j+1,1):resolution_in_x:shrunk_polytopes(i).vertices(j,1);
                vq = interp1([shrunk_polytopes(i).vertices(j+1,1),shrunk_polytopes(i).vertices(j,1)],...
                                [shrunk_polytopes(i).vertices(j+1,2),shrunk_polytopes(i).vertices(j,2)],...
                                xq);
                vq = flip(vq);
                xq = flip(xq);
            end
            new_verts_this_side = [xq; vq]';
            new_verts = [new_verts; new_verts_this_side];
        end
        shrunk_polytopes(i).vertices = [new_verts;new_verts(1,:)];
        shrunk_polytopes(i).xv = new_verts(:,1)';
        shrunk_polytopes(i).yv = new_verts(:,2)';
        if flag_do_plot
            figure(102938)
            hold on
            plot(shrunk_polytopes(i).vertices(:,1),shrunk_polytopes(i).vertices(:,2),'rx')
            legend('','orig. vertices','dense vertices')
        end
    end
end
