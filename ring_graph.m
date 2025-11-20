function [A_ring,tof, distances]  = ring_graph(r, n, k)
    % Create ring network with n nodes and k-neighbours
    % INPUTS:
    % r             radius of ring
    % n             number of nodes     
    % k             number of neighbours on each side
    % OUTPUTS: 
    % A_ring        adjacency matrix of ring graph
    % tof           time of flight
    
    % speed of light
    c = 299792.458;

    if k>n
        error("k value is too large. k<n")
    elseif rem(k,2)~=0
        error("k must be even.")
    end
    
    
    k=k/2;

    A_ring=zeros(n);
    tof = zeros(n);
    distances = zeros(n);
    angles = (2*pi/n).*(1:k);
    chord = 2*(r/1000)*sin(angles./2);
    
    for i=1:n
        range = min(n-i, k);
        A_ring(i,i+1:i+range)=1;
        A_ring(i+1:i+range, i)=1; 
        distances(i,i+1:i+range) = chord(1:range);
        distances(i+1:i+range, i) = chord(1:range);
        tof(i,i+1:i+range) = chord(1:range)./c;
        tof(i+1:i+range, i) = (chord(1:range)./c).';
    end

    % close ring
    for i=1:k+1
        A_ring(i, end-k+i:end)=1;
        A_ring(end-k+i:end, i)=1;
        distances(i, end-k+i:end) = flip(chord(i:end));
        distances(end-i+1, 1:k-i+1) = chord(i:end);
        tof(i, end-k+i:end) = flip(chord(i:end)./c);
        %tof(end-k+i:end, i) = chord(i:end)./c;
        tof(end-i+1, 1:k-i+1) = chord(i:end)./c;
    end
        
    

end