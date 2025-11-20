clear all
close all 


tested_conn = [2, 4, 6,8];
%tested_conn = 4;
type_comm = "3 Rings";
time_step = 250;
window_length_basic = 10;

% selection of plots to output
plot_dict = [0;     % network plot
            1;      % adjacency & laplacian spectra
            1;      % adjacency & laplacian EIGENGAP spectra
            1;      % eigengap comparison
            1;      % eigenvector centralities
            1];     % eigenvectors heatmap

% plotting parameters 
GM = 0.39860*10^(6);
RE = 6378;
v= [-149.4564,32.2465];
flip_plotting=false;
extra_zoom = false;
colors_query=colormap(parula(size(tested_conn,2)+1));
close()
if isscalar(tested_conn)
   shape_query= {'-o'};
   colors_query = [0, 0, 0];
else
    shape_query = {'-o', '-*', '-+', '-s', '->', '-^', '-.', '-<', '-|', '-diamond', '-v'};
    if length(tested_conn)>length(shape_query)
        ratio = ceil(length(tested_conn)/length(shape_query));
        shape_query = repmat(shape_query, [1, ratio]);
    end
    shape_query = shape_query(1,1:length(tested_conn));
end
max_gap_A = 0;
max_gap_L = 0;
max_min_A = [nan,nan];
max_min_L = [nan,nan];
if type_comm == "2 Rings"
    data_spec = "2_rings_%G_conn_ANALYSIS.mat";
    nb_rings = 2;
    alt_data = [2, 4, 6, 8;
        30000,25000,12000,8000]; 
    tested_conn = alt_data(1,ismember(alt_data(1,:),tested_conn));
    alt = alt_data(2,ismember(alt_data(1,:),tested_conn));
    if isscalar(alt)
        window_length = window_length_basic;
    else        
        orbital_speed = sqrt(GM./(RE+alt));
        ratio_speeds = max(orbital_speed)./orbital_speed;        
        window_length = window_length_basic.*ratio_speeds;
    end
    % zoom in on spectra from 0 to :
    x_zoom_laplacian = 50;
    x_zoom_adjacency = 10;
    names = strcat(repmat("2 Rings, \delta =  ",[length(tested_conn) 1]), num2str(tested_conn.'));
elseif type_comm == "3 Rings"
    data_spec = "3_rings_%G_conn_ANALYSIS.mat";
    nb_rings = 3;
    alt_data = [2, 4, 6;
        18000,10000,5000]; 
    tested_conn = alt_data(1,ismember(alt_data(1,:),tested_conn));
    alt = alt_data(2,ismember(alt_data(1,:),tested_conn));
    if isscalar(alt)
        window_length = window_length_basic;
    else        
        orbital_speed = sqrt(GM./(RE+alt));
        ratio_speeds = max(orbital_speed)./orbital_speed;        
        window_length = window_length_basic.*ratio_speeds;
    end
    % zoom in on spectra from 0 to :
    x_zoom_laplacian = 100;
    x_zoom_adjacency = 10;   
    names = strcat(repmat("3 Rings, \delta =  ",[length(tested_conn) 1]), num2str(tested_conn.'));
end

% path creation to use Communities of Dynamical Influence scripts
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


for i=1:length(tested_conn)
    % define average ring connectivity delta_R
    connectivity  = tested_conn(i);

    % Load network based on test parameters
    data = sprintf(data_spec, connectivity);
    sat_data=load(data);
    A_all = sat_data.A;
    A_type = class(A_all);
    n = sqrt(length(A_all(: ,1)));
    tof = zeros(n);
    pos_plot = sat_data.pos_plot;
        
    % extract the adjacency matrix 
    % average over multiple time steps
    A_window = zeros(n,n,round(window_length(i)));
    indeg_window = zeros(1,n,round(window_length(i)));
    pos_window = zeros(3,n,round(window_length(i)));
    for k=1:round(window_length(i))
        A_window(:,:,k) = time_step_adajacency(A_all, A_type, time_step+k, n, tof);
        pos_window(:,:,k) = [pos_plot(1, :, time_step);pos_plot(2, :, time_step);pos_plot(3, :, time_step)];
        A_window(:,:,k) = double(A_window(:,:,k));
        indeg_window(:,:,k) = sum(A_window(:,:,k));
    end
    A = mean(A_window,3);
    A(A<0.1) = 0;
    indeg = sum(A);
    pos = mean(pos_window,3);

    % create laplacian
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
    colour_comm_all = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);
    colour_comm_outer = vecnorm(vec_L_sorted(:, g_rank(1,i)+1:g_rank(2,i)), 2,2);
    colour_comm_inner = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);

    colour = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);
    
    % if final plot --> add legends
    final_plot = false;
    if i==length(tested_conn)
        final_plot=true;
    end    
    
    colour_A = vecnorm(vec_A_sorted(:, 1:find(eigen_gap_A==min(eigen_gap_A))), 2,2);
    G = graph(A);
    G.Edges.LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
    figure(200+i)
    p=plot(G,'XData', pos(1,:),'YData', pos(2,:),'ZData', pos(3,:), 'MarkerSize',15, 'LineWidth', G.Edges.LWidths, NodeCData=colour_A, NodeLabel={});
    xlabel("X")
    ylabel("Y")
    p.EdgeColor  = [200, 200, 200]./256;
    fontsize(25,"points")
    grid on 
    box on
    c = colorbar;
    c.Label.String = 'Centrality, C_i';

    [max_min_A, max_min_L, max_gap_A, max_gap_L]=eigen_plots(plot_dict, shape_query{i}, colors_query(i,:), v, A, pos, 1, eig_A_sorted, max_min_A, eigen_gap_A, max_gap_A, ...
    x_zoom_adjacency, eig_L_sorted, max_min_L, eigen_gap_L, max_gap_L, x_zoom_laplacian, colour, normalised_vec_L, final_plot, names);
    

    % test clustering
    %[CDI_L, ~, ~, leaders_export] = CDI_nopath(L,"L",1:30);

    % figure()
    % node_colours = CDI_L;
    % p1= plot(G,'XData', pos(1,:),'YData', pos(2,:), 'ZData', pos(3,:),'MarkerSize',15, 'LineWidth', G.Edges.LWidths,  'NodeCData',node_colours);    
    % p1.EdgeColor  = [150, 150, 150]./256;
    % p1.NodeLabel = {};
    % title("Laplacian")
    % 
   
end
