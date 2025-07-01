function [Umat]=Uvect2Umat(Uvect,nu)

%   [Umat]=Uvect2Umat(Uvect,nu)
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


N=length(Uvect)/nu;
Umat=reshape(Uvect,nu,N)';

end