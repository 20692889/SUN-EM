function [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, p,q,r_cp,k,quad_pts,sing)
% INT_PQ returns the four integrals for faces (=elements) p and q,
% corresponding to [eqs.34a-34d,RWG82]. Note that these integrals are
% evaluated in normalized coordinates, and are not scaled by the area.
% The overall matrix element is assembled from these elsewhere.
% Note also that p is the field point, q the source point. The integrals are
% performed over triangle q.
% Note further that xi, eta and zeta are equivalent to
% lambda_1, lambda_2, and lambda_3 respectively.
% If flag sing is set, then the singular terms are evaluated using a
% special integration rule.

% Author: D B Davidson, Dec 2009.
% Corrections for singular intergral evaluation: 1 June 2010 DBD.

%global ELEMENTS NODE_COORD

qnodes = elements(q,:);
n1 = node_coord(qnodes(1),:);
n2 = node_coord(qnodes(2),:);
n3 = node_coord(qnodes(3),:);
area = tri_area3D(n1,n2,n3); 

Ipq=0;
Ipq_xi=0;
Ipq_eta=0;
Ipq_zeta=0;

if p==q && sing
    [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = intg_sing_SGF(k,r_cp,n1,n2,n3,3,4);
    Ipq = Ipq/(2*area); 
    % The factor of 2 here is required to get good results, but it is not
    % clear why this is correct. If removed, then the w/2 below should also be
    % changed to w, but the results are then out by 2. 
    Ipq_xi = Ipq_xi/(2*area);
    Ipq_eta = Ipq_eta/(2*area);
    Ipq_zeta = Ipq_zeta/(2*area);
    % 2A factor above required to give result in simplex coordinates - see
    % [eq 31. RWG82]. (Function intg_sing_SGF returns the LHS thereof).

else
    [w,lambda] = tri_quad(quad_pts);
    w=w/2; % Rule must be correctly normalized.
    
    %r_cp
    for nn=1:quad_pts
        r_prime = lambda(nn,1)*n1 + lambda(nn,2)*n2 + lambda(nn,3)*n3; % [eq.30,RWG82]
        R_p = norm(r_cp-r_prime); % [eq.27,RWG82]
        
        % DJdbg --> plot r_cp and r_prime
        %plot(r_prime(1,1),r_prime(1,2),'or','MarkerFaceColor','g','MarkerSize',10);        
        %plot(r_cp(1,1),r_cp(1,2),'or','MarkerFaceColor','k','MarkerSize',10);   
        
        GF  = exp(-1i*k*R_p)/R_p; % Green's function
        %w(nn)
        Ipq     = Ipq+w(nn)*GF;
        %lambda(nn,1)
        Ipq_xi  = Ipq_xi+w(nn)*lambda(nn,1)*GF;
        %lambda(nn,2)
        Ipq_eta = Ipq_eta+w(nn)*lambda(nn,2)*GF;
        %Ipq_zeta = Ipq_zeta+w(nn)*lambda(nn,3)*GF; % Code gives same
        %answers as below.
    end
    Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
end
end



function [area]=tri_area3D(a,b,c)
% TRI_AREA_3d returns the unsigned area of the triangle located in 3D space
% with vertices a, b, c. 

% The algorithm (but not code) comes from
% [Press] Press et al, "Numerical Recipes: the Art of Scientific
% Computing", 3rd ed, CUP 2007, p.1115.
vec_area = 0.5 *cross((a-c),(b-c));
area = norm(vec_area);
end



function [w,lambda] = tri_quad (n)
% This function returns the weights and evaluation points (in simplex
% coordinates lambda) from D.A.Dunavant, "Gaussian quadrature formulas for
% triangles", IJNME, vol 21,  1985, pp. 1129-1148.

switch n
    case 1 % 1 point rule, degree of precision 1
        w=zeros(1,1);
        alpha = zeros(1,1);
        beta = zeros(1,1);
        gamma = zeros(1,1);
        w(1) = 1;
        alpha(1) = 1/3;
        beta(1) = 1/3;
        gamma(1) = 1/3;
    case 3 % 3 point rule, degree of precision 2
        w=zeros(3,1);
        alpha = zeros(3,1);
        beta = zeros(3,1);
        gamma = zeros(3,1);
        w(1) = 1/3; 
        w(2) = w(1);
        w(3) = w(1);
        alpha(1) = 2/3;
        beta(1) = 1/6;
        gamma(1) = 1/6;
        alpha(2) = gamma(1);
        beta(2) = alpha(1);
        gamma(2) = beta(1);
        alpha(3) = beta(1);
        beta(3) = gamma(1);
        gamma(3) = alpha(1);                
    case 6 % symmetric 6 point rule, degree of precision 4
        w = zeros(6,1);
        alpha = zeros(6,1);
        beta = zeros(6,1);
        gamma = zeros(6,1);
        w(1) =     0.223381589678011;
        w(2) = w(1);
        w(3) = w(1);
        w(4) =     0.109951743655322;
        w(5) = w(4);
        w(6) = w(4);
        alpha(1) = 0.108103018168070;
        beta(1)  = 0.445948490915965;
        gamma(1) = beta(1);
        alpha(2) = beta(1);
        alpha(3) = beta(1);
        beta(2)  = alpha(1);
        beta(3)  = beta(1);
        gamma(2) = beta(1);
        gamma(3) = alpha(1);
        alpha(4) = 0.816847572980459;
        beta(4)  = 0.091576213509771;
        gamma(4) = beta(4);
        alpha(5) = beta(4);
        alpha(6) = beta(4);
        beta(5)  = alpha(4);
        beta(6)  = beta(4);
        gamma(5) = beta(4);
        gamma(6) = alpha(4);
    case 12 % symmetric 12 point rule, degree of precision 6
        w = zeros(12,1);
        alpha = zeros(12,1);
        beta = zeros(12,1);
        gamma = zeros(12,1);
        w1 = 0.116786275726379;
        w2 = 0.050844906370207;
        w3 = 0.082851075618374;
        a1 = 0.501426509658179;
        b1 = 0.249286745170910;
        a2 = 0.873821971016996;
        b2 = 0.063089014491502;
        a3 = 0.053145049844817;
        b3 = 0.310352451033784;
        g3 = 0.636502499121399;
        w(1) = w1;
        w(2) = w1;
        w(3) = w1;
        w(4) = w2;
        w(5) = w2;
        w(6) = w2;
        w(7) = w3;
        w(8) = w3;
        w(9) = w3;
        w(10) = w3;
        w(11) = w3;
        w(12) = w3;
        alpha(1)  = a1; beta(1) = b1;  gamma(1) = b1;
        alpha(2)  = b1; beta(2) = a1;  gamma(2) = b1;
        alpha(3)  = b1; beta(3) = b1;  gamma(3) = a1;
        alpha(4)  = a2; beta(4) = b2;  gamma(4) = b2;
        alpha(5)  = b2; beta(5) = a2;  gamma(5) = b2;
        alpha(6)  = b2; beta(6) = b2;  gamma(6) = a2;
        alpha(7)  = a3; beta(7) = b3;  gamma(7) = g3;
        alpha(8)  = a3; beta(8) = g3;  gamma(8) = b3;
        alpha(9)  = b3; beta(9) = g3;  gamma(9) = a3;
        alpha(10) = b3; beta(10) = a3; gamma(10) = g3;
        alpha(11) = g3; beta(11) = a3; gamma(11) = b3;
        alpha(12) = g3; beta(12) = b3; gamma(12) = a3;        
    otherwise
        error 'Unimplemented rule'
end
lambda(1:n,1) = alpha;
lambda(1:n,2) = beta;
lambda(1:n,3) = gamma;
end