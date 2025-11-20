clear all 
close all 


% tested ring connectivity
tested_conn = 4;

% selection of plots to output
plot_dict = [0;     % network plot
            1;      % adjacency & laplacian spectra
            1;      % adjacency & laplacian EIGENGAP spectra
            0;      % eigenvector centralities
            0];     % eigenvectors heatmap

% analysis parameters
% number of nodes per hub
nb_nodes_per_hub = [10, 10, 10, 10;
     10, 11, 12, 13;
     10,16,17,18;
     20, 10,11,12];
% hubs combination analysed (1 hub, then 2 hubs, ...)
%nb_hubs_analysed= [1,2,3,4];
nb_hubs_analysed= [4];

% ring coordinate dimensions 
r = 11000; % radius
n = 200; % #nodes

if length(nb_nodes_per_hub)<length(nb_hubs_analysed)
    error("Number of nodes for each hub does not correspond to number of hubs")
end

% plotting parameters 
separate_plots = true; % separate the ajdacency & laplacian plots
x_zoom_adjacency = 6*ones(1,size(nb_nodes_per_hub,1));
x_zoom_laplacian = 60*ones(1,size(nb_nodes_per_hub,1));
names = strcat(repmat("c = ",[length(nb_hubs_analysed) 1]), num2str(nb_hubs_analysed.'));
extra_zoom = false;
x_extra_zoom_laplacian = ((max(nb_nodes_per_hub)-1).*nb_hubs_analysed)-2;
k=1;
m=1;
angle = 0:(360/(n)):360;
x_data = 1*cosd(angle(1:end-1));
y_data = 1*sind(angle(1:end-1));
pos = [x_data;y_data];
v=[];
max_gap_A = 0;
max_gap_L = 0;
max_min_A = [0,0];
max_min_L = [0,0];

% path creation to use Communities of Dynamical Influence scripts
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% plot colour and marker shapes
colors_query=colormap(parula(size(nb_nodes_per_hub,2)+1));
close()
if isscalar(nb_hubs_analysed)
   shape_query= {'-o'};
   colors_query = [0, 0, 0];
else
    shape_query = {'-o', '-*', '-+', '-s', '->', '-^', '-.', '-<', '-|', '-diamond', '-v'};
    if length(tested_conn)>length(shape_query)
        ratio = ceil(length(tested_conn)/length(shape_query));
        shape_query = repmat(shape_query, [1, ratio]);
    end
    shape_query = shape_query(1,1:length(nb_hubs_analysed));
end

% run through hub combinations
for i=1:size(nb_nodes_per_hub,1)
    for j=1:length(nb_hubs_analysed)
        nb_hubs =nb_hubs_analysed(j);
        connectivity  = tested_conn;       
        
        % create ring network
        A = ring_graph(r, n, connectivity);
        % embed the hubs
        A(1:nb_nodes_per_hub(i, 1), 1:nb_nodes_per_hub(i, 1)) = ones(nb_nodes_per_hub(i, 1),nb_nodes_per_hub(i, 1));
        A = A - diag(diag(A));
        hub_pt = round(n/nb_hubs);
        for ii=1:nb_hubs
            A((ii-1)*hub_pt+1:(ii-1)*hub_pt+nb_nodes_per_hub(i,ii), (ii-1)*hub_pt+1:(ii-1)*hub_pt+nb_nodes_per_hub(i,ii)) = ones(nb_nodes_per_hub(i,ii),nb_nodes_per_hub(i,ii));
            A = A - diag(diag(A));               
        end
        % create laplacian
        indeg = sum(A,1); 
        D = sum(A,2).*toeplitz([1, zeros(1, n-1)]);
        L = D - A;

        % extract eigenspectra
        [eig_A_sorted, vec_A_sorted, eigen_gap_A] =  eigengap(A, "A");
        [eig_L_sorted, vec_L_sorted, eigen_gap_L] =  eigengap(L, "L");
        normalised_vec_L =  abs(vec_L_sorted);
        
        % assess relevant communities
        nb_sig_gaps = 2;
        % most significant eigengaps 
        g = mink(eigen_gap_L(1:round(n/2),1), nb_sig_gaps);
        g_rank(:,i) = find(ismember(eigen_gap_L,g));

        % eigenvector centrality
        colour = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);

        % if final plot --> add legends
        final_plot = false;
        if j==length(nb_hubs_analysed)
           final_plot=true;
        end    
        if i>1 && or(or(plot_dict(1)==1,plot_dict(4)==1), plot_dict(5)==1)
            disp("Cannot plot the ajdacency/laplacian AND centrality/heatmap")
            break
        end   
        
        colour_A = vecnorm(vec_A_sorted(:, 1:j), 2,2);
        G = graph(A);
        G.Edges.LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
        figure(200+i)
        p=plot(G,'XData', pos(1,:),'YData', pos(2,:), 'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=colour_A, NodeLabel={});
        xlabel("X")
        ylabel("Y")
        p.EdgeColor  = [200, 200, 200]./256;
        fontsize(25,"points")
        grid on 
        box on
        c = colorbar;
        c.Label.String = 'Centrality, C_i';
        
        % plot outputs
        [max_min_A, max_min_L, max_gap_A, max_gap_L]=eigen_plots(plot_dict, shape_query{j}, colors_query(j,:), v, A, pos, i, eig_A_sorted, max_min_A, eigen_gap_A, max_gap_A, ...
        x_zoom_adjacency(i), eig_L_sorted, max_min_L, eigen_gap_L, max_gap_L, x_zoom_laplacian(i), colour, normalised_vec_L, final_plot, names);
    
    end
end



