function [max_min_A, max_min_L, max_gap_A, max_gap_L]=eigen_plots(plot_dict, shape_query, colors_query, v, A, pos, counter, eig_A_sorted, max_min_A, eigen_gap_A, max_gap_A, ...
    x_zoom_adjacency, eig_L_sorted, max_min_L, eigen_gap_L, max_gap_L, x_zoom_laplacian, colour, normalised_vec_L, final_plot, names)
% Plotting function given dictionary
% Plots adjacency eigenspectrum, Laplacian eigenspectrum, Eigenvector
% centrality
% 
% plot_dict             boolean dictionary of plots to output
%                           > network plot
%                           > adjacency & laplacian spectra
%                           > adjacency & laplacian EIGENGAP spectra
%                           > eigenvector centralities
%                           > heatmap eigenvectors
    % create graph
    G = graph(A);
    G.Edges.LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    indeg = sum(A);
    n = length(A);
    
    if plot_dict(2)==1
        % ADJACENCY
        figure((counter-1)*5+1)
        hold on
        plot(1:n, eig_A_sorted ,shape_query, 'Color',colors_query,'LineWidth', 2)  
        max_min_A(1) = max(max(eig_A_sorted((1:x_zoom_adjacency))), max_min_A(1));
        max_min_A(2) = min(min(eig_A_sorted((1:x_zoom_adjacency))), max_min_A(2));
        if final_plot
            xlim([-inf x_zoom_adjacency])
            ylim([floor(max_min_A(2)*0.8) ceil(max_min_A(1)*1.2)])
            ylabel("Adjacency Eigenvalue, \lambda^A")
            xlabel("Rank index, i")
            grid on
            box on
            if ~sum(colors_query)==0
                legend(names,'Location', 'northeast');
            end
            fontsize(18,"points")
        end

        % LAPLACIAN
        figure((counter-1)*5+2)
        hold on  
        plot(1:n, eig_L_sorted,shape_query, 'Color',colors_query,'LineWidth', 2)
        max_min_L(1) = max(max(eig_L_sorted((1:x_zoom_laplacian))), max_min_L(1));
        max_min_L(2) = min(min(eig_L_sorted((1:x_zoom_laplacian))), max_min_L(2));
        if final_plot
            xlim([-inf x_zoom_laplacian])
            ylim([floor(max_min_L(2)*0.8) ceil(max_min_L(1)*1.2)])
            ylabel("Laplacian Eigenvalue, \lambda^L")
            xlabel("Rank index, i")
            grid on
            box on
            if ~sum(colors_query)==0
                legend(names,'Location', 'northeast');
            end
            fontsize(18,"points")
        end
    end
    
    if plot_dict(3)==1
        % ADJACENCY EIGENGAP
        figure((counter-1)*5+3)
        hold on
        plot(1:n-1, eigen_gap_A,shape_query, 'Color',colors_query,'LineWidth', 2)
        xlim([-inf x_zoom_adjacency])
        max_gap_A = min(min(eigen_gap_A((1:x_zoom_adjacency))), max_gap_A);
        if final_plot
            xlim([-inf x_zoom_adjacency])
            ylim([max_gap_A*1.2 inf])
            ylabel("Adjacency Eigengap, g^A")
            xlabel("Rank index, i")
            grid on
            box on
            if ~sum(colors_query)==0
                legend(names,'Location', 'southeast');
            end
            fontsize(18,"points")
        end
        
    
        % LAPLACIAN EIGENGAP
        figure((counter-1)*5+4)
        hold on 
        plot(1:n-1, eigen_gap_L,shape_query, 'Color',colors_query,'LineWidth', 2)
        xlim([-inf x_zoom_laplacian])
        max_gap_L = min(min(eigen_gap_L((1:x_zoom_laplacian))), max_gap_L);
        if final_plot
            xlim([-inf x_zoom_laplacian])
            ylim([max_gap_L*1.2 inf])
            ylabel("Laplacian Eigengap, g^L")
            xlabel("Rank index, i")
            grid on
            box on 
            if ~sum(colors_query)==0
                legend(names,'Location', 'southeast');
            end
            fontsize(18,"points")
        end
        
    end
    if plot_dict(1)==1
        %plot network coloured with indegrees
        figure()
        if isempty(pos)
            p=plot(G,'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=indeg, NodeLabel={});
            colormap("parula")
            cbh=colorbar;
            cbh.Label.String = 'Connectivity';
        elseif size(pos, 1)==2
            p=plot(G,'XData', pos(1,:),'YData', pos(2,:), 'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=indeg, NodeLabel={});
            xlabel("X")
            ylabel("Y")
            colormap("parula")
            cbh=colorbar;
            cbh.Label.String = 'Connectivity';
        elseif size(pos, 1)==3
            p=plot(G,'XData', pos(1,:),'YData', pos(2,:), 'ZData', pos(3,:),'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=indeg, NodeLabel={});
            xlabel("X position")
            ylabel("Y position")
            zlabel("Z position")
            colormap("parula")
            cbh=colorbar;
            cbh.Label.String = 'Connectivity';
            set(cbh,'YTick',[min(round(indeg)):1:max(round(indeg))])
        end        
        p.EdgeColor  = [200, 200, 200]./256;
        grid on 
        box on        
        fontsize(25,"points")
        if ~isempty(v)
            view(v)
        end
    end

    if plot_dict(4) ==1
        % CENTRALITY
        figure()
        if isempty(pos)
            p=plot(G,'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=colour, NodeLabel={});
        elseif size(pos, 1)==2
            p=plot(G,'XData', pos(1,:),'YData', pos(2,:), 'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=colour, NodeLabel={});
            xlabel("X")
            ylabel("Y")
        elseif size(pos, 1)==3
            p=plot(G,'XData', pos(1,:),'YData', pos(2,:), 'ZData', pos(3,:),'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=colour, NodeLabel={});
            xlabel("X")
            ylabel("Y")
            zlabel("Z")
        end        
        p.EdgeColor  = [200, 200, 200]./256;
        fontsize(25,"points")
        if ~isempty(v)
            view(v)
        end
        grid on 
        box on
        c = colorbar;
        c.Label.String = 'Centrality, C_i';
    end

    if plot_dict(5)==1
        % EIGENVECTOR HEATMAP
        figure()
        h = heatmap(normalised_vec_L(:, 1:x_zoom_laplacian));
        Ax = gca;
        Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
        Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
        colormap(flipud(gray))
        grid off 
        h.XLabel='Node IDs';
        h.YLabel = 'Sorted Eigenvectors, /tau';
        colorbar;
        hs = struct(h);
        ylabel(hs.Colorbar, "Absolute Magnitude of Eigenvectors, |\upsilon_i|");
        fontsize(18,"points")
        
    end

end
