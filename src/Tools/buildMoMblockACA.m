 function [ZBlock,u,v] = buildMoMblockACA (Const,num,Solver_setup, zMatrices, freq, m, n, ObservRWGs, SourceRWGs)

% ZBlock         = U x V;
% size(Z_approx) = length(ObservRWGs) x length(SourceRWGs)
%
% 16-Apr-2013 Rob Maaskant, Göteborg

% 24.03.2014: Danie made some changes to get the script compatible with
%             MATLAB V 6.5.0.180913a (R13):
%             - Cannot use [~,]
%             - Function find only has a single parameters
%             - Use a different implementation of buildMoMblock
%
% 01.04.2014: Danie added now the use of a fixed number of iterations.
% 15.04.2014: Danie added efficient implementation (see Const.ACAalg = 2).
%             Returns ZBlock = U x V;
% 15.04.2014: Danie added efficient implementation (see Const.ACAalg = 2)
%             Returns U, V;
% 25.06.2017: Danie added also 

% References:
% -----------
% [1] K. Zao, M.N. Vouvakis and J-F. Lee, "The Adaptive Cross Approximation Algorithm for 
%     Accelerated Method of Moments Computations of EMC Problems", 
%     IEEE Trans. on Electrom. Compat., Vol. 47, No. 4, Nov. 2005.

%TO-DOs: The ACA can be accelerated in some cases by selecting a different starting and ending
%        vector. (See notes made for CEM'17 course in Barcelona, as shown by Alex from Politechnica
%        de Barcelona)

Tol = Const.useACAtol;

row(1)          = 1; 
i       = buildMoMblock(Const,num,Solver_setup, zMatrices, freq, ObservRWGs(1), SourceRWGs);
v(1,:)         = i.values ;% Build one row of the MoM block
        
Rapprox        = sparse(length(ObservRWGs),length(SourceRWGs));
Rapprox(1,:)    = v(1,:);
[tmp,col]       = max(abs(Rapprox(1,:)));
v(1,:)          = Rapprox(1,:)/Rapprox(1,col);
i       =  buildMoMblock(Const,num,Solver_setup, zMatrices, freq, ObservRWGs, SourceRWGs(col));
u(:,1)          = i.values(:,col);% Build one column of the MoM block

Rapprox(:,col)  = u(:,1);


Mrwg = length(ObservRWGs);
Nrwg = length(SourceRWGs);

% If a fixed maximum number of iterations are used, then we can use an
% optimised version where no convergence criteria are chosen

% -----------------------------------------
% Use convergence criteria to terminate ACA
if (Const.useACAfixedIterations < 0) 

    if (Const.ACAalg == 1)
    % --- Original algorithm (non-eficient)    
        k = 1;    
        ZBlock = u*v; % Slow
        NormZ  = sum(sum(conj(ZBlock).*ZBlock)); % Slow (|~Z|^2)
        uv     = (u(:,1)'*u(:,1))*(v(1,:)*v(1,:)'); % (|u|^2|v|^2)
        
		while uv > Tol.^2*NormZ
            
            % Find a row that was not taken yet to build Zblock
            [tmp,Index]   = sort(abs(Rapprox(:,col(k))));    
            count        = 0;            
            while ~isempty(find(Index(end-count)==row)); count=count+1; end
            
            k           = k + 1;
            row(k)      = Index(end-count);
            i   =  buildMoMblock (Const,num,Solver_setup, zMatrices,freq, ObservRWGs(row(k)), SourceRWGs);
            v(k,:)      =  i.values(row(k),:);
  
            Rapprox(row(k),:)   = v(k,:)-ZBlock(row(k),:);
            [tmp,Index]  = sort(abs(Rapprox(row(k),:)));
            
            count=0;        
            while ~isempty(find(Index(end-count)==col)); count=count+1; end
            
            col(k)              = Index(end-count);
            v(k,:)              = Rapprox(row(k),:)/Rapprox(row(k),col(k));
            i             = buildMoMblock(Const,num,Solver_setup, zMatrices,freq, ObservRWGs, SourceRWGs(col(k)));
            u(:,k)              = i.values(:,col(k));
           
            Rapprox(:,col(k))   = u(:,k)-ZBlock(:,col(k));
            u(:,k)              = Rapprox(:,col(k));
            uv                  = (u(:,k)'*u(:,k))*(v(k,:)*v(k,:)');
            ZBlock              = u*v; % Slow
            NormZ               = sum(sum(conj(ZBlock).*ZBlock)); % Slow
		end
        
    elseif (Const.ACAalg >= 2)
    % --- Improved algorithm (eficient)
    %     k is rank
            k = 1;
            ZBlockNorm = norm(u(:,1),'fro')*norm(v(1,:),'fro');
            ZBlockNorm = ZBlockNorm*ZBlockNorm;
            uv = ZBlockNorm;
            
            while uv > Tol.^2*ZBlockNorm
                
                % Find a row that was not taken yet to build Zblock
                [tmp,Index]   = sort(abs(Rapprox(:,col(k))));    
                count        = 0;                
                while ~isempty(find(Index(end-count)==row)); count=count+1; end                
                k           = k + 1;
                row(k)      = Index(end-count); 
                
                % Update the (Ik)th row of the approximate error matrix (R=Zmom - Zapprox^(k-1))
                i        = buildMoMblock (Const,num,Solver_setup, zMatrices, freq, ObservRWGs(row(k)), SourceRWGs); 
                Zmom     = i.values(row(k),:);
     
                Zapprox  = u(row(k),:)*v;
                Rapprox(row(k),:) = Zmom - Zapprox;
                                 
                % Find a column that was not taken yet to build Zblock
                [tmp,Index] = sort(abs(Rapprox(row(k),:)));
                count=0;        
                while ~isempty(find(Index(end-count)==col)); count=count+1; end                
                col(k)      = Index(end-count); 
                                
                v(k,:)      = Rapprox(row(k),:)/Rapprox(row(k),col(k)); 
                vk          = v(k,:);                

                % Update the (Jk)th column of the approximate error matrix
                i  = buildMoMblock (Const,num,Solver_setup, zMatrices, freq, ObservRWGs, SourceRWGs(col(k)));
                Zmom     = i.values(:,col(k));
              
                Zapprox  = u*v(1:k-1,col(k));
                Rapprox(:,col(k)) = Zmom - Zapprox;
                
                u(:,k)      = Rapprox(:,col(k));                
                uk          = u(:,k);
                
                % Calculate the matrix norm
               	for i = 1:(k-1) % (kth Iteration)
				    ZBlockNorm=ZBlockNorm+2*norm(uk'*u(:,i),'fro') * norm(v(i,:)*vk','fro');
                end%for
                
			    uv=norm(u(:,k),'fro')*norm(v(k,:),'fro'); uv=uv*uv;
			    ZBlockNorm=ZBlockNorm+uv;
            end
 end

% -----------------------------------------
% Use fixed iterations for ACA
else 
    k = 1;
    ZBlock = u*v; % Slow
    while k < Const.useACAfixedIterations
        
        % [1] L(8): Find the next row index I_(k+1) such that
        %           |R(I_(k+1),Jk)|=max_i(|R(i,Jk))|
        [tmp,Index]  = sort(abs(Rapprox(:,col(k))));
        count        = 0;
        while ~isempty(find(Index(end-count)==row)); count=count+1; end        
        k           = k + 1;
        row(k)      = Index(end-count);
        i     = buildMoMblock (Const,num,Solver_setup, zMatrices, freq, ObservRWGs(row(k)), SourceRWGs);
        v(k,:)      = i.values(row(k),:);% [1] Related to L(3) (kth Iteration)

        Rapprox(row(k),:)   = v(k,:)-ZBlock(row(k),:); % [1] L(1) (kth Iteration)
        
        % [1] L(2) (kth Iteration)
        [tmp,Index]  = sort(abs(Rapprox(row(k),:)));        
        count=0;        
        while ~isempty(find(Index(end-count)==col)); count=count+1; end        
        col(k)              = Index(end-count);
        
        v(k,:)              = Rapprox(row(k),:)/Rapprox(row(k),col(k)); % [1] L(3) (kth Iteration)
        i             = buildMoMblock(Const,num,Solver_setup, zMatrices, freq, ObservRWGs, SourceRWGs(col(k)));
        u(:,k)              = i.values(:,col(k)); % [1] L(5) (kth Iteration)

        Rapprox(:,col(k))   = u(:,k)-ZBlock(:,col(k));  % [1] L(4) (kth Iteration)
        u(:,k)              = Rapprox(:,col(k)); % [1] L(5) (kth Iteration)
        
        ZBlock              = u*v; % Slow
    end%for
end

% -----------------------------------------
% Some useful debugging output
%if (true)
%    message(Const,sprintf('    (ACA debug) Processing Z%d%d k=%d size(U):(%d,%d) size(V):(%d,%d)',...
%        m,n,k,size(u,1),size(u,2),size(v,1),size(v,2)));
%end%if

if (false)
    % Also some simplified output for further processing in e.g. POSTFEKO
    check_row = 1;
    if (m == check_row)
        message(Const,sprintf('%d,%d',n,k));
    end %if (m == check_row)
end%if

% -----------------------------------------
% Set the return value here (note, for Alg = 3 we only return U and V)
if (Const.ACAalg < 3)
    ZBlock = u*v; %Slow
else
    ZBlock = [];
end

% norm(ZBlock - Zoff,'fro')./norm(Zoff,'fro')*100
return