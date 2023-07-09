% To implement Self Optimal Clustering(SOC) technique 
% for images containing RGB data points or for any (n x k) matrix
% WITH ALL THE DATA POINTS INCLUDED

clc
clear
disp('input the image name in the format --> imread(image name, image format); ');
a = input('Enter image data : ');
[f, h, k] = size(a);
nk = input('no. of clusters required : ');                      % nk is no. of clusters required

for j = 1:f
    for i = 1:h
        R(j,i) = a(j,i,1);                                  % extracting R value from original matrix a
        G(j,i) = a(j,i,2);                                  % extracting G value from original matrix a
        B(j,i) = a(j,i,3);                                  % extracting B value from original matrix a
    end
end
n = f*h;                                                    % n is total no. of data points
for j = 1:n
    x(j,:) = [R(j),G(j),B(j)];
end                                                         % So x is a (n x 3) matrix for RGB 

[fac]=factorcal(x,nk,1);
[result] = soc(x,nk,fac);

[s] = silhouette(double(x), result.idx);
[S, GSS] = slht(s, result.idx, result.n, result.m, nk);
GSS

% [PI, SI] = valid(result.dd,result.cc_norm,result.part.^2,nk)
% [ADI] = adu(result.clst,nk,result.m)
