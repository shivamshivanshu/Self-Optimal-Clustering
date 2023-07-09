function [fac]=factorcal(x,nk,iter)

    flag=0;
    factor=ones(10,nk);          
    while(iter<=10)
        iter
        [result] = soc(x,nk,factor(iter,:));
        % disp('d1 value');     result.d1 -> threshold distance

        [s] = silhouette(double(x), result.idx);
        % disp('S value');      S -> Global Silhouette
        [S, GS(iter)] = slht(s, result.idx, result.n, result.m, nk);
        
        % disp('PI and SI values');
        % [PI, SI] = valid(result.dd,result.cc_norm,result.part.^2,nk);

        if(min(result.m)==0)
            flag==1;
            break;
        end
        for g=1:nk-1
            for gg=g+1:nk
                if(result.d1(g)==result.d1(gg))
                    flag=1;
                end
            end
        end
        if(flag==1)
            break;
        end

        polym=lagrangepoly(result.d1,S);
        polym(1,nk)=polym(1,nk)-1;
        r=roots(polym);
        sumn=zeros(size(r));
        for i=1:nk-1
            for j=1:nk
            sumn(i)=sumn(i)+polym(nk-j+1)*r(i)^(j-1);
            end
        end
        [mm,label]=min(abs(sumn));
        % disp('dmax value');
        dmax=abs(r(label));
        factor(iter+1,:)=(dmax./result.d1);
        iter=iter+1;
    end
    [mm,label]=max(GS);
    fac=factor(label,:);
    
end
