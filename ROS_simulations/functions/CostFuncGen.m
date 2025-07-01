function [H_sim,h_sim]=CostFuncGen(T_sim,S_sim,Q,R,P)
%
% Cost function generation for LTV systems. Computes cost functions V(U) for
% dense MPC formulation all along the reference trajectory. Each "single 
% cost function" starting from a GENERIC point x0 along the reference 
% trajectory is:
%   V(U) = 0.5*U'*H*U + (h*x0)'*U + const
% from a typical quadratic cost function:
%   V(Y,U) = 0.5 * sum{k=0,N-1}( yk'*Q*yk + uk'*R*uk ) + 0.5 * yN'*P*yN 
% where:
%  -Y=[y(0)' y(1)' y(2)' ... y(N)']: outputs sequence from time k=0 to
%                                    k=N | dim=((N+1)*ny,1)
%  -U=[u(0)' u(1)' u(2)' ... u(N-1)']: input sequence from time k=0 to time
%                                      k=N-1 | dim=(N*nu,1)
% for a the following predicted dynamics from x0 (dense formulation):
%   Y=T*x0+S*U
% 
% Syntax:
%
%   [H_sim,h_sim]=CostFuncGen(T_sim,S_sim,Q,R,P)
%
% Inputs
%  -T_sim.........T prediction matrices (from initial state) | dim=((N+1)*ny,nx,N_sim+1)
%  -S_sim.........S prediction matrices (from input) | dim=((N+1)*ny,N*nu,N_sim+1) 
%  -Q,R,P.........weighting matrices of the cost function | dimQ=(nx,nx), 
%                 dimR=(nu,nu), dimP=(nx,nx)
% Outputs:
%  -H_sim,h_sim...matrices for the dense formulation of the cost function.
%                 They are collected in a 3D matrix, and each "slice" (:,:,k0)
%                 corresponds to a certain point along the reference
%                 trajectory


% Extract dimensions

nx=length(P);
nu=length(R);
N=size(S_sim,2)/nu;
N_sim=size(T_sim,3)-1;

% Define "dense" weight matrices (weight matrices are constant here)

Q_bar=kron(eye(N),Q);
Q_bar=blkdiag(Q_bar,P);
R_bar=kron(eye(N),R);

% Build cost function's matrices from each k0

H_sim=zeros(N*nu,N*nu,N_sim+1);
h_sim=zeros(N*nu,nx,N_sim+1);
for k0=0:N_sim    
    H_sim(:,:,k0+1)=S_sim(:,:,k0+1)'*Q_bar*S_sim(:,:,k0+1)+R_bar;
    h_sim(:,:,k0+1)=S_sim(:,:,k0+1)'*Q_bar'*T_sim(:,:,k0+1);

    H_sim(:,:,k0+1)=(H_sim(:,:,k0+1)+H_sim(:,:,k0+1)')/2;   % enforce symmetry on H

    clc
    fprintf('Cost matrices generation: %.1f%%',k0/N_sim*100)
end
clc

% NOTE: only H and h matrices are given, while "const" is not needed. In fact
%    const=0.5*x0'*(T'*Q_bar*T)*x0;
% however, x0 is given (measured state), so it's constant, it's not
% something we need to minimize in our optimisation, so we can not account
% for it.

end