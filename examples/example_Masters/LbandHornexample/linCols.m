%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%-------------------Linearly independent column compressor function------------------------
%---------------------------------------------------------------------------------------------
function [Xsub,idx]=linCols(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted rows of X
% idx:  The indices (into X) of the extracted columns
     if ~nnz(X) %X has no non-zeros and hence no independent columns
         Xsub=[]; idx=[];
         return
     end
     if nargin<2, tol=1e-30; end %later specify tol in Solver_setup.
       [Q, R, E] = qr(X,0);
       
       %Q = Orthogonal factor; R = Upper-triangular factor; E = Permutation information
       %X*E = Q*R
       
       if ~isvector(R)
        diagr = abs(diag(R)); 
       else
        diagr = R(1);   
       end
       %Rank estimation
       r = find(diagr >= tol*diagr(1),1, 'last'); %rank estimation
       idx=sort(E(1:r));
       Xsub=X(:,idx);
       
