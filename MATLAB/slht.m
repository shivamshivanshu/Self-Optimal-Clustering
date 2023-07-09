function [S, GS] = slht(s, idx, n, m, nk)

    S_const(1:nk) = 0;

    for r = 1:n
        S_const(idx(r)) = S_const(idx(r)) +  s(r);     % to add up individual Silhouette value of data points in each of clusters to obtain Silhouette value for each cluster
    end

    % To calculate cluster Silhouette 
    for j = 1:nk
        S(j) = (1/m(j)) * S_const(j);
    end

    % To calculate a Global Silhouette value
    GS = (1/nk) * sum(S);

end

