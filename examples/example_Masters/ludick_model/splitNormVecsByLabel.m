function [Idx] = splitNormVecsByLabel(Solver_setup)

char = Solver_setup.metallic_triangles_labels;
% norms = Solver_setup.triangle_normal_vector;
Idx(1) = 1; 
temp_char = char(1);

for i = 2: size(char,1)
   bool = strcmp(char(i),temp_char);
    if bool == 0
        Idx(end+1)=i;
        temp_char = char(i);
        
    end
end

% normvecs = norms(Idx(1):Idx(2), :);
% k=1;
% for k = 1:size(Idx,1)
%     if 
% 
% for i = 1: size(norms,1)
%     normvecs1 = norms(Idx(i):Idx(i+1)-1,:);
%     normvecs2 = 
%     