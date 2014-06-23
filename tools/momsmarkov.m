function [mmean,msdev] = momsmarkov(X,pie) 

% SAMPLEMOMS Computes stationary moments of  controlled Markov process
%
% Usage: 
%    [mmean,msdev] = momsmarkov(X,pie)
%
% Inputs:
%  X      : n by N matrix of N column-vectors sized n by 1 each, where
%           n is the number of discrete states and N is the number of
%           variables
%  pie    : stationary density of controlled Markov process
% Ouputs:
%  mmean  : Nx1 vector of mean of N time series vectors
%  msdev  : Nx1 vector of standard deviations of N time series vectors

[n,N] = size(X);
mmean = X'*pie;
msdev = zeros(N,1);

for i=1:N
    msdev(i)=sqrt(((X(:,i)-mmean(i)).^2)'*pie);
end