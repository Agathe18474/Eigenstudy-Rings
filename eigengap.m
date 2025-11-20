function [eig_M_sorted, vec_M_sorted, eigen_gap] =  eigengap(M, type)
% Calculate and sort the eigenvalues, eigengaps, and eigenvector of a
% matrix
% Inputs: 
%   M       matrix M,either the adjacency or Laplacian
%   type    indicate what matrix M is
%               > A 
%               > L
    [vec_M, eig_M] = eig(M);
    
    if type == "norm_L"
        [eig_M_sorted, Ind] = sort(diag(real(eig_M)), 'ascend');
        vec_M_sorted = real(vec_M(:,Ind));
    elseif type == "correlation" || type == "A" || type == "L" ||type == "B"
        [eig_M_sorted, Ind] = sort(diag(real(eig_M)), 'descend');
        vec_M_sorted = real(vec_M(:,Ind));
    end
        
    eigen_gap = diff(real(eig_M_sorted));

end