function [result] = soc(x,nk,factor)

    % To implement Improved Mountain Clustering Technique(IMC) 
    % for images containing RGB data points or for any (n x k) matrix
    % WITH ALL THE DATA POINTS INCLUDED
    
    g = (1:nk);
    [n,k] = size(x); 

    % STEP 1
    % to Normalize each dimension of hyperspace
    x_min = min(x);
    x_max = max(x);
    z = x_max - x_min;

    for j = 1:n
        y = x(j,:) - x_min;
        if (z ~= 0)
            w = double(y)./double(z);
        else
            w(1,1:k) = 0;
        end
        for i = 1:k
            u(j,i) = w(1,i);                         % range of values is (0,1); to find normalized value of each unit of x
        end
    end
    U = u;
    u = double(u);

    m(g) = 0;                        % no. of data points in each cluster ('m' here corresponds to elements in 'SL')
    t(1,1:nk+1) = 0;                 % no. of data points in the beginning before each clustering ('t' here corresponds to elements in 'sl')
    t(1) = n;                        % initially before 1st clustering all x data points are counted in t(1) 
    sl(:,1) = (1:n);                 % initially before 1st clustering all x data points locations are counted in sl(:,1)

    for v = 1:nk
        if (t(v) ~= 0)
            d(v) = 0;
            
            % STEP 2
            % to Determine the parameter d1 for each window
            for j = 1:t(v)
                if (sum(x(sl(j,v),:)) ~= 0)
                    d(v) = d(v) + double(min(x(sl(j,v),:)))./double(sum(x(sl(j,v),:)));
                end
            end

            d1(v) = ((1/(2.*t(v))).*d(v)).*factor(v);

            % STEP 3
            % to Calculate the potential value of each point using a mountain function
            for r = 1:t(v)
                P(r,v) = 0;
                for j = 1:t(v)
                    P(r,v) = P(r,v) + exp(-(((u(r,:) - u(j,:))*(u(r,:) - u(j,:))') ./ (d1(v).^2)));
                end
            end

            % STEP 4
            % to Select the first cluster center according to the highest value of P1  
            [ymax(v),zmax(v)] = max(P(:,v));                          % z stores the location of P at which it is maximum every time
            cc_norm(v,:) = u(zmax(v),:);                              % cc_norm stores cluster center value
            cc_norm=double(cc_norm);
            
            % STEP 5
            % to Assign concerned data points to the first cluster 
            for r = 1:t(v)
                if (((u(r,:) - cc_norm(v,:))*(u(r,:) - cc_norm(v,:))') <= d1(v)) 
                    m(v) = m(v) + 1;                                % m(v) stores no. of data points in vth cluster
                    SL(m(v),v) = sl(r,v);                           % SL(m(v),v) stores positions/locations of data points in vth cluster
                    for i = 1:k
                        clst(m(v),i,v) = x(sl(r,v),i);              % only those x values(data points) being stored in clst which satisfie the condition
                        c_disp(sl(r,v),i,v) = x(sl(r,v),i);         % those x values stored in c_disp which satisfy the condition, else 0 is stored as value at those unsatisfied locations
                        ex(sl(r,v),i) =255;
                    end    
                else
                     t(v+1) = t(v+1) + 1;                                  % t(v+1) stores the no. of data points(dp) left after clustering v, or say dp before (v+1)th clustering 
                    sl(t(v+1),v+1) = sl(r,v);
                    for i = 1:k
                         u(t(v+1),i) = u(r,i);                             % u is updated here as having only left out data points after each clustering
                        ex(sl(r,v),i) = x(sl(r,v),i);
                        c_disp(sl(r,v),i,v) = 255;
                    end
                end
            end
            c_disp(:,:,v+1) = ex;
        end                                     % this 'end' is for 'if (t(v) ~= 0)'
    end                                         % this 'end' is for 'for v = 1:nk' loop

    % To distribute left out data points among already segmented clusters
    if (t(v+1) ~= 0)                                        % to avoid popping up of errors in case no data points are left
        for r = 1:t(nk+1)
            for v = 1:nk
                D(v,r) = ((u(r,:) - cc_norm(v,:))*(u(r,:) - cc_norm(v,:))');   
            end
        end
        [ymin,zmin] = min(D);                               % gives a row matrix having all min. distance values for each data point from existing cluster centers

        for r = 1:t(nk+1)
            v = zmin(r);                                    % zmin tells the cluster corresponding to each left out data points 
            m(v) = m(v) + 1;                    
            SL(m(v),v) = sl(r,nk+1);            
            for i = 1:k
                clst(m(v),:,v) = x(sl(r,nk+1),:);           % to append left out data points into respective clusters
                c_disp(sl(r,nk+1),:,v) = x(sl(r,nk+1),:);   % to append left our data points into respective displace matrix
            end
        end 
    end                                                     % 'end' is for 'if' loop (t(v+1) ~= 0)

    % -------------------------------------------------------------------------

    % to calculate idx
    for r = 1:n
        for v = 1:nk
            if (m(v) ~= 0)
                dd(r,v) =((U(r,:) - cc_norm(v,:))*(U(r,:) - cc_norm(v,:))');
                for j = 1:m(v)
                    if (x(r,:) == clst(j,:,v))
                        idx(r) = v;
                    end
                end
            end
        end
    end


    % To calculate partition matrix

    for v = 1:nk
        if (m(v) ~= 0)
          part = zeros(n,v);
        end
    end


    for j = 1:n
        [M,label] = min(dd(j,:));
        part(j,label) = 1;
    end
    
    % -------------------------------------------------------------------------

    result.dd=dd;
    result.part=part;
    result.cc_norm=cc_norm;
    result.idx=idx;
    result.c_disp=c_disp;
    result.clst=clst;
    result.m=m;
    result.SL=SL;
    result.n=n;
    result.d1=d1;
    
end                         
