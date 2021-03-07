function [Zmn,U,V] = calcZmn(Const,num,Solver_setup, zMatrices, freq, m, n, ObservRWGs, SourceRWGs)
    %calcZmn
    %   Usage:
    %       [zMatrices] = calcZmn(Const, zMatrices, m, n, ObservRWGs, SourceRWGs)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       zMatrices
    %           The Z-matrices data
    %       freq
    %           The frequency index for which this MoM matrix should be extracted
    %       m,n
    %           The matrix row and column subscripts
    %
    %   Output Arguments:
    %       Zmn
    %           The coupling matrix Zmn that was calculated with the MoM (fast or slow algorithm) /
    %           the ACA
    %
    %   Description:
    %       Calculates the coupling matrix Zmn using either the ACA or the
    %       MoM (using either a fast or slow implementation. For the ACA, the calculation is based
    %       on the partially pivoted
    %       ACA algorithm ([1], [2]).
    %
    %   =======================
    %   Written by Danie Ludick on March 24, 2014.
    %   Last updated on June 25, 2017.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    %   References:
    % [1] Zhao, K., Vouvakis, M.N. and Lee, J.-F.: The adaptive cross approximation algorithm for
    %     accelerated method of moments computations of emc problems. Electromagnetic Compati-
    %     bility, IEEE Transactions on, vol. 47, no. 4, pp. 763{773, 2005.
    %
    % [2] Rjasanow, Sergej, and Olaf Steinbach. The fast solution of boundary integral equations.
    %     Springer, 2007. (See pp. 126 - 127, Partially pivoted ACA algorithm).

    error(nargchk(7,9,nargin));
    %narginchk(7,7);

    % Initialise U and V (only set for ACA)
    U = [];
    V = [];

    if (Const.useACA)

        [Zmn,U,V] = buildMoMblockACA (Const,num,Solver_setup, zMatrices,freq, m, n, ObservRWGs, SourceRWGs); % Rob Maaskant implementation
        
        % Debug output below:
        if (true)
            check_row = 1;
            if (m == check_row)
                % Some debug output
               % message(Const,sprintf('\n  Calculating the Z%d%d matrix using the ACA alg. %d',m,n,Const.ACAalg));

                if (Const.ACAalg == 3)
                    % Construct here Zmn as this is empty when returned by the buildMoMblocACA for this efficient 
                    % version of the ACA algorithm in order to calcualte the accuracy and memory usage below.
                    Zmn = U*V;
                end%if

                % Calculate here the memory usage of the coupling submatrices, U and V (consider that we will eventually
                % be using a robuust storage functionality for the H-matrix)
                Zmn_memUsage = byteSize(Zmn);
                U_memUsage = byteSize(U);
                V_memUsage = byteSize(V);
                
               % message(Const,sprintf('    Memory usage of Z%d%d=%s , U=%s, V=%s',m,n,Zmn_memUsage,U_memUsage,V_memUsage));

                err=norm(zMatrices.values(ObservRWGs,SourceRWGs)-Zmn,'fro')/norm(zMatrices.values(ObservRWGs,SourceRWGs),'fro');%*100;
                %message(Const,sprintf('    ACA approximation error = %5.5f', err));

            end %if (m == check_row)
        end%if (false/true)
    else
        if (Const.debug)
            % Some debug output
            message(Const,sprintf('  Calculating the Z%d%d matrix using the MoM',m,n));
        end%if
        Zmn = extractZmnfromFEKOmatfile(Const, zMatrices, freq, ObservRWGs, SourceRWGs);
        %Zmn = buildMoMblock (Const, Solver_setup, zMatrices, freq, ObservRWGs, SourceRWGs);
        %Zmn = extractZmnfromFEKOmatfile(Const, zMatrices, freq, ObservRWGs, SourceRWGs);
        if (Const.debug)
            % Some debug output
            % None required, just get the same IF test here as for the above.
        end%if
    end%if

    % TO-DO: Danie, format the outpout (any timing?) or flop count?. Store
    % to array for each of the submatrices and plot afterwards.
    %message(Const,sprintf('Finished DGFM solver in %f sec.',dgfm.solTime));

