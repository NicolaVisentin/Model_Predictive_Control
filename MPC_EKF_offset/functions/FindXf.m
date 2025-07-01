function [Xf,Z]=FindXf(A,B,K,xlb,xub,ulb,uub,terminal,solver)
%
% Calculates the various sets for the given system and parameters.
%
% Most inputs are self-explanatory. terminal should be either the string
% 'termeq' to use a terminal equality constraint, or 'lqr' to use the
% maximum LQR-invariant set.
% Solver is the one used to find Olqr: 'linpr' for linprog, anything else
% for YALMIP with quadprog.
% K is the optimal Riccati gain found with [-K,P]=dlqr(A,B,Q,R).
%
% Output Xn is a cell array defining each of the \bbX_n sets. Xn{1} is the
% terminal set, Xn{2} is \bbX_1, etc. V is a cell array with each entry giving
% the extreme vertices of the corresponding Xn (as a 2 by N matrix). Z is a
% structure defining the feasible space.


% Define constraints for Z = {[x; u] | G x + H u + \psi <= 0}

Nx=size(A, 1);
[Az,bz]=hyperrectangle([xlb;ulb],[xub;uub]);
Z=struct('G',Az(:,1:Nx),'H',Az(:,(Nx+1):end),'psi',bz);

% Decide terminal constraint

switch terminal
    case 'termeq'
        % Equality constraint at the origin
        Xf = [0; 0; 0];
    
    case 'lqr'
        % Build feasible region considering x \in X and Kx \in U
        % Control input
        [A_U, b_U] = hyperrectangle(ulb, uub);
        A_lqr = A_U*K;
        b_lqr = b_U;
        
        % State input
        [A_X, b_X] = hyperrectangle(xlb, xub);
        Acon = [A_lqr; A_X];
        bcon = [b_lqr; b_X];
    
        % Use LQR-invariant set.
        Xf = struct();
        ApBK = A + B*K; % LQR evolution matrix.
        [Xf.A, Xf.b] = calcOinf(ApBK, Acon, bcon, solver);
        [~, Xf.A, Xf.b] = removeredundantcon(Xf.A, Xf.b);

end

end