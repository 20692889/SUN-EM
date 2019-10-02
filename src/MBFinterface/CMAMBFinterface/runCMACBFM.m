function [cbfm] = runCMACBFM(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs)
    cbfm = [];
    Nmom = Solver_setup.num_mom_basis_functions;
    numArrEls = Solver_setup.num_finite_array_elements;
    cbfm.numSols = 1;
    cbfm.Isol = complex(zeros(Nmom, 1));

    numR = (mbfs.numMBFs)*(Solver_setup.mom_basis_functions_per_array_element);
    numC = (mbfs.numMBFs)*(Solver_setup.mom_basis_functions_per_array_element);

    cbfm.Zred = complex(zeros(numR, numC));
    cbfm.Vred = complex(zeros(numR, 1));

    for p=1:numArrEls
        PindxEnd = 0;
        PindxStart = 0;
        numMBFsP = 0;
        domain_p_bfs = Solver_setup.rwg_basis_functions_domains{p};
        N_p = length(domain_p_bfs);

        domainPCbfs = complex(zeros(N_p, mbfs.numMBFs));%3));%mbfs.numRedMBFs(p, 1)));
        for(red_ii=1:mbfs.numMBFs)%3)%mbfs.numRedMBFs(p, 1))
            domainPCbfs(:, red_ii) = mbfs.orthSol(:, red_ii);%RedIsol(domain_p_bfs, red_ii, p, 1);
        end
        Vrwg = yVectors.values(domain_p_bfs, 1);
        for q = 1:numArrEls
            domain_q_bfs = Solver_setup.rwg_basis_functions_domains{q};
            N_q = length(domain_q_bfs);
            domainQCbfs = complex(zeros(N_q, mbfs.numMBFs));%3));%mbfs.numRedMBFs(q, 1)));
            for red_ii=1:mbfs.numMBFs%3%mbfs.numRedMBFs(q, 1)
                domainQCbfs(:, red_ii) = mbfs.orthSol(:, red_ii);%RedIsol(domain_q_bfs, red_ii, q,1);
            end
            [Z] = zMatrices.values(domain_p_bfs, domain_q_bfs, 1);

            numMBFsP = 0;
            for domain = 1:(p-1)
                 numMBFsP = numMBFsP + mbfs.numMBFs;%3;%mbfs.numRedMBFs(domain,1);
            end
            PindxStart = numMBFsP+1;
            PindxEnd   = (PindxStart - 1) + mbfs.numMBFs;%3; %mbfs.numRedMBFs(p,1);

            numMBFsQ = 0;
            for domain = 1:(q-1)
                 numMBFsQ = numMBFsQ + mbfs.numMBFs;%3;%mbfs.numRedMBFs(domain,1);
            end
            QindxStart = numMBFsQ+1;
            QindxEnd   = (QindxStart - 1) + mbfs.numMBFs;%3;%mbfs.numRedMBFs(q,1);
            cbfm.Zred(PindxStart:PindxEnd, QindxStart:QindxEnd) = (domainPCbfs)'*Z*domainQCbfs;
        end

        cbfm.Vred(PindxStart:PindxEnd) = (domainPCbfs)' * Vrwg;
    end
    [L,U] = lu(cbfm.Zred);
    b = L\cbfm.Vred;
    cbfm.Ired = U\b;
    for p=1:numArrEls
        domain_p_bfs = Solver_setup.rwg_basis_functions_domains{p};
        if (p > 1)
            offset = offset + mbfs.numMBFs;%3;%mbfs.numRedMBFs(p-1, 1);
        else
              offset = 0;
        end

        for red_ii=1:mbfs.numMBFs%3%mbfs.numRedMBFs(p, 1)
             fak = cbfm.Ired(offset + red_ii);
            cbfm.Isol(domain_p_bfs, 1) = cbfm.Isol(domain_p_bfs, 1) + fak.*mbfs.orthSol(:, red_ii);%mbfs.RedIsol(domain_p_bfs,red_ii,p,1);
        end
    end
   cbfm.relError = calculateErrorNormPercentage(xVectors.Isol(:,1), cbfm.Isol(:,1));
   %message_fc(Const, sprintf('Rel. error norm. for Sol. %d of %d compared to FEKO sol. %f percent.',1, 1, cbfm.relError));
   message_fc(Const, sprintf('Using %d MBFs there is an error norm of %f percent compared to FEKO sol.',  mbfs.numMBFs,cbfm.relError));
end

