function [Uvect]=Umat2Uvect(Umat)
%
%   [Uvect]=Umat2Uvect(Umat)
%
% Converts a (N*nu,1) vector in form:
%   Uvect=[u(0)
%          u(1) 
%          u(2) 
%           .
%           .
%           .
%          u(N-1)]
% into its corresponding (N,nu) matrix form:
%   Umat=[u(0)'
%         u(1)'
%         u(2)'
%          .
%          .
%          .
%         u(N-1)']


N=size(Umat,1);
nu=size(Umat,2);
Uvect=reshape(Umat',nu*N,1);

end