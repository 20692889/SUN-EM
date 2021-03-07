function [lin_vec, Idx] = GeoSurface(X)











lin_vec = [];
temp_vec = X(1,:);
lin_vec(1,:) = X(1,:);
Idx(1) = 1;
d =0;
%i=2;
for i = 2: size(X,1) 
    
    dif_vec = temp_vec - X(i,:);
    d = find(abs(dif_vec) > 0.1);
    if d > 0
        Idx(end+1) = i-1;
        temp_vec = X(i,:);
        lin_vec(end+1,:) = temp_vec;
    end
end



    