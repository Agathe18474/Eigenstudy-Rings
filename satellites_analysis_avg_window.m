close all 
clear all

% simulation parameters
connectivity = 4;
type_comm = "3 Rings";
time_step = 250;
step_avg_window = 1;
avg_window = 1:step_avg_window:20;

% selection of plots to output
plot_dict = [0;     % network plot
            0;      % adjacency & laplacian spectra
            0;      % adjacency & laplacian EIGENGAP spectra
            0;      % eigenvector centralities
            0];     % heatmap of eigenvectors

% plotting parameters 
GM = 0.39860*10^(6);
RE = 6378;
v= [-149.4564,32.2465];
flip_plotting=false;
extra_zoom = false;
colors_query=colormap(parula(length(avg_window)+1));
shape_query = {'-o', '-*', '-+', '-s', '->', '-^', '-.', '-<', '-|', '-diamond', '-v'};
max_gap_L=0;
max_gap_A=0;
if length(avg_window)>length(shape_query)
    ratio = ceil(length(avg_window)/length(shape_query));
    shape_query = repmat(shape_query, [1, ratio]);
end
shape_query = shape_query(1,1:length(avg_window));
if flip_plotting==true
    connectivity = flip(connectivity);
    colors_query = flip(colors_query(1:end-1,:));
    shape_query = flip(shape_query);
end

if type_comm == "2 Rings"
    data_spec = "2_rings_%G_conn_ANALYSIS.mat";
    nb_rings = 2;
    alt = 25000;
    x_zoom_laplacian = 60;
    x_zoom_adjacency = 10;
    gaps = nan(length(avg_window), 30);
elseif type_comm == "3 Rings"
    data_spec = "3_rings_%G_conn_ANALYSIS.mat";
    nb_rings = 3;
    alt = 10000;
    x_zoom_laplacian = 100;
    x_zoom_adjacency = 10;    
    gaps = nan(length(avg_window), 70);
end


% path creation to use Communities of Dynamical Influence scripts
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


for i=1:length(avg_window)
    % Load network based on test parameters
    data = sprintf(data_spec, connectivity);
    sat_data=load(data);
    A_all = sat_data.A;
    A_type = class(A_all);
    n = sqrt(length(A_all(: ,1)));
    tof = zeros(n);
    pos_plot = sat_data.pos_plot;
    names(1,i) = string(avg_window(i));   
        
    % extract the adjacency matrix 
    % average over multiple time steps
    A_window = zeros(n,n,round(avg_window(i)));
    indeg_window = zeros(1,n,round(avg_window(i)));
    pos_window = zeros(3,n,round(avg_window(i)));
    for k=1:round(avg_window(i))
        A_window(:,:,k) = time_step_adajacency(A_all, A_type, time_step+k, n, tof);
        pos_window(:,:,k) = [pos_plot(1, :, time_step);pos_plot(2, :, time_step);pos_plot(3, :, time_step)];
        A_window(:,:,k) = double(A_window(:,:,k));
        indeg_window(:,:,k) = sum(A_window(:,:,k));
    end
    A = mean(A_window,3);
    A(A<0.1) = 0;    
    pos = mean(pos_window,3);

    % create laplacian
    D = sum(A,2).*toeplitz([1, zeros(1, n-1)]);
    L = D - A;

    % extract eigenspectra
    [eig_A_sorted, vec_A_sorted, eigen_gap_A] =  eigengap(A, "A");
    [eig_L_sorted, vec_L_sorted, eigen_gap_L] =  eigengap(L, "L");
    normalised_vec_L =  abs(vec_L_sorted);

    gaps(i, :) = eigen_gap_L(1:size(gaps,2));
    gaps(i, gaps(i,:)>=-0.1)=nan;

    % assess relevant communities
    nb_sig_gaps = 2;
    % most significant eigengaps 
    g = mink(eigen_gap_L(1:round(n/2),1), nb_sig_gaps);
    g_rank(:,i) = find(ismember(eigen_gap_L,g));

    % eigenvector centrality    
    comm_all = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);
    comm_outer = vecnorm(vec_L_sorted(:, g_rank(1,i)+1:g_rank(2,i)), 2,2);
    comm_inner = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);

    colour = vecnorm(vec_L_sorted(:, 1:g_rank(2,i)), 2,2);

    final_plot = false;
    if i==length(avg_window)
        final_plot=true;
    end    
    
    [max_gap_A, max_gap_L]=eigen_plots(plot_dict, shape_query{i}, colors_query(i,:), v, A, pos, eig_A_sorted, eigen_gap_A, max_gap_A, ...
    x_zoom_adjacency, eig_L_sorted, eigen_gap_L, max_gap_L, x_zoom_laplacian, colour, normalised_vec_L, final_plot, names);

   
    figure(1)
    hold on
    scatter(avg_window(i)*ones(1, size(gaps, 2)), 1:size(gaps,2), 40, reshape(abs(gaps(i,:)), 1, []), 'filled')
    colormap(flip(gray))
    if final_plot
        xlim([0, avg_window(end)+1])
        xlabel("Averaging Window Duration")
        ylabel("Laplacian Eigengaps Ranks")
        c=colorbar;
        c.Label.String  = "Eigengap Magnitude";
        fontsize(18,"points")
        box on        
    end
   
   

end


