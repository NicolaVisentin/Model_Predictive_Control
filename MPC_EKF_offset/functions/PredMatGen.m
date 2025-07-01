function [T_sim,S_sim]=PredMatGen(LTV,N)
%
% Prediction matrices generation for LTV system. Computes prediction 
% matrices T and S for dense MPC formulation all along the reference 
% trajectory. Each "single prediction" starting from a GENERIC point x0 
% along the reference trajectory is:
%   Y = T*x0 + S*U
% where:
%  -Y=[y(0)' y(1)' y(2)' ... y(N)']: outputs sequence from time k=0 to
%                                    k=N | dim=((N+1)*ny,1);
%  -x0: initial condition vector | dim=(nx,1);
%  -U=[u(0)' u(1)' u(2)' ... u(N-1)']: input sequence from time k=0 to time
%                                      k=N-1 | dim=(N*nu,1)
% 
% Syntax:
%
%   [T_sim,S_sim]=PredMatGen(LTV,N)
%
% Inputs
%  -LTV.....structure containing matrices A_ref (nx,nx,N_ref), B_ref 
%           (nx,nu,N_ref), C (ny,nx) that define the discrete time linear 
%           time varying system. N_ref is the number of timesteps composing 
%           the reference trajectory (i.e. the trajectory of the "virtual 
%           robot"). The MPC will be able to work from time 0 to time 
%           N_sim=N_ref-N, where N is the receiding horizon.
%  -N.......number of discrete time steps of the receiding horizon.
%
% Outputs:
%  ! TO HAVE S AND T AS STATE PREDICTION MATRICES, JUST USE C=I, OBVIOUSLY
%  -T_sim...Ty output prediction matrices (initial state) |
%           dim=((N+1)*ny,nx,N_sim+1), where N_sim=N_ref-N is the time our
%           MPC is able to work.
%  -S_sim...Sy output prediction matrices (input) | dim=((N+1)*ny,N*nu,N_sim+1), 
%           where N_sim=N_ref-N is the time our MPC is able to work.
%  ! Basically, each "slice" (:,:,k0) of T_sim or S_sim matrices contains 
%  the T or S dense prediction matrix that starts from time k0 and
%  "predicts" up to time k0+N. This of course can be done from time k0=0 to
%  time k0=N_ref-N=N_sim.


% Extract matrices

A_ref=LTV.A_ref;
B_ref=LTV.B_ref;
C=LTV.C;

% Compute dimensions

nx=size(A_ref,1);
nu=size(B_ref,2);
ny=size(C,1);
N_ref=size(A_ref,3);
N_sim=N_ref-N;

% Create matrices T and S

T_sim=zeros((N+1)*ny,nx,N_sim+1);
S_sim=zeros((N+1)*ny,N*nu,N_sim+1);

for k0=0:N_sim
    T_sim(1:ny,:,k0+1)=C*eye(nx);
    
    for k=1:N    
        % update matrix T
        T_sim(k*ny+1:(k+1)*ny,:,k0+1)=A_ref(:,:,k0+k)*T_sim((k-1)*ny+1:k*ny,:,k0+1);
    
        % update matrix S
        if k==1
            S_sim(k*ny+1:(k+1)*ny,1:nu,k0+1)=C*B_ref(:,:,k0+k);
        else
            for ii=0:k-2
                S_sim(k*ny+1:(k+1)*ny,ii*nu+1:(ii+1)*nu,k0+1)=A_ref(:,:,k0+k)*S_sim((k-1)*ny+1:k*ny,ii*nu+1:(ii+1)*nu,k0+1);
                S_sim(k*ny+1:(k+1)*ny,(ii+1)*nu+1:(ii+2)*nu,k0+1)=B_ref(:,:,k0+k);
            end
        end
    end
    clc
    fprintf('Prediction matrices generation: %.1f%%',k0/N_sim*100)
end
clc

end