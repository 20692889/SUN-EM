function [Uvector, Vvector] = AdapadtiveCrossApproximation(Zmn, m, n)


r = min(m, n);
I_row = [];
I_col = [];
%Uvector = zeros(m,r);
%Vvector= zeros(r,n);
Uvector = [];
Vvector = [];
I_row(1) = 1;

Rerr = zeros(m,n);
Rerr(I_row(1),:) = Zmn(I_row(1),:);
Z_old = zeros(m,n);
Z_new = zeros(m,n);
[M,j] = max(abs(Rerr(I_row(1),:)));
I_col(1) = j;
Vvector(1,:) = Rerr(I_row(1),:)/Rerr(I_row(1),I_col(1));
Rerr(:,I_col(1))= Zmn(:,I_col(1));
Uvector(:,1) = Rerr(:,I_col(1));
Z_old = ((norm(Uvector(:,1),'fro')).^2)*((norm(Vvector(1,:),'fro')).^2);
[M,j] = max(abs(Rerr(:,I_col(1))));
I_row(2) = j;
k=2;
sum_uv = 0;
sum_vu =0;
sum_d = 0;
UV_frob = 0;
flag = 'true';

   % while(flag)
   for i=1:r
    
       for i=1:k-1
           sum_uv = sum_uv + Uvector(I_row(k), i)*Vvector(i,:);
       end
       
       Rerr(I_row(k),:) = Zmn(I_row(k),:) - sum_uv;
       [M,j] = max(abs(Rerr(I_row(k),:)));
       I_col(k) = j;
       Vvector(k, :) = Rerr(I_row(k),:)/Rerr(I_row(k), I_col(k));
       for i=1:k-1
           sum_vu = sum_vu + Vvector(i, I_col(k))*Uvector(:,i);
       end
       Rerr(:,I_col(k))= Zmn(:,I_col(k))-sum_vu;
       Uvector(:,k) = Rerr(:,I_col(k));
      % for j = 1:k-1
      %    a = abs((Uvector(:,j)')*(Uvector(:,k))); 
      %    b = abs((Vvector(j,:)')*(Vvector(k,:)));
          %d= dot(a,b);
      %    sum_d = sum_d + (a*b);
      % end
      % Z_new = Z_old + (2*sum_d) + ((norm(Uvector(:,k),'fro')).^2)*((norm(Vvector(k,:),'fro')).^2);
      % UV_frob = norm(Uvector(:,k),'fro')*norm(Vvector(k,:),'fro');
      % if(UV_frob <= (0.01*sqrt(Z_new)))
      %     flag = 'false';
      % else
           [M,j] = max(abs(Rerr(:,I_col(k)))); 
           k=k+1;
           I_row(k) = j;
       %end
   end   
end    


% r = min(m, n);
% Uapprox = zeros(m,r);
% Vapprox = zeros(n,r);
% Zapproxk
% Ifind = [];
% Jfind= [];
% totK = 0;
% eps_err = 10^-3; 
% Rx = zeros(m,n);
% iter = 'true' ;
% sum_uv = 0;
% dotSum = 0;
% 
% %for I_itter = Ifind
% % for J_itter = Jfind
%     while(iter)
%         for k = k + 1
%             Ifind(1) = 1;
%             Zapprox(0) = 0;
%             for i = 1:k
%                 sum_uv = Uapprox(i,Ifind(k))*Vapprox(i,:);
%             end
%             
%             Rx(Ifind(k), :) = Zmnx(Ifind(k),:) - sum_uv;
%             [M,K] = max(Rx(Ifind(k), :));
%             Jfind(k) = K;
%             Vapprox(k,:) = Rx(Ifind(k),:)./Rx(Ifind(k), Jfind(k));
%             Rx(:, Jfind(k)) = Zmnx(:, Jfind(k));
%             Uapprox(:,k) = Rx(:, Jfind(k));
%             [M,K] = max(Rx(: , Jfind(k)));
%             Ifind(k) =  K;
% 
%            % totk = totK +1;
%             
%             for j = 1:k-1
%                 
%                 Zapprox(k-1) = Zapprox(k-1) + Uapprox(:,i)*Vapprox(i,:);
%                 summ1 = abs((Uapprox(:,j).')*(Uapprox(:,k)));
%                 summ2 = abs((Vapprox(j,:).')*(Vapprox(k,:)));
%                 dotSum = dotSum + dot(summ1,summ2);   
%             end
%             
%             dotSum = dotSum*2;
%             Zk_norm2 = (norm(Zapprox(k-1),'fro'))^2 + dotSum + ((norm(Uapprox(:,k),'fro'))^2)*((norm(Vapprox(k,:),'fro'))^2) ;    
% 
%             Zapprox_norm(k) = sqrt(Zk_norm2);
%             uk = norm(Uapprox(:,k),'fro');
%             vk = norm(Vapprox(k,:),'fro');
%             
%                 if((uk*vk) <= (eps_err*Zapprox_norm(k)))
%                     iter = 'false';
%                 end
%                
%         end
%         
%      end
% 
%  
   
   %syms i
   %Zapprox = symsum((Uapprox(n,Ifind(k))*Vapprox(n)),i,0,r);
 %  end
 



















 