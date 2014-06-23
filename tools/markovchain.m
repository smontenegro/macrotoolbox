function [Z,PI] = markovchain(N,p,q,e,m,a)

% MARKOVCHAIN Markov process and transition matrix 
% 
%Purpose:    Follows Rouwenhorst to find a Markov chain whose sample paths 
%            approximate a family of random walk processes with serial 
%            correlation p+q-1. Depending on p and q the generated process 
%            could be AR(1), unit root, or ARCH.
%            See: Cooley's RBC book, chapter 10 by K. G. Rouwenhorst
%
%Usage:     [Z, PI] = markovchain(N,p,q,e)
%
%Input:      N   scalar, number of nodes for Z, usually odd number
%            p   scalar, governs cond. prob of moving up after low state
%            q   scalar, governs cond. prob of moving down after low state
%            e   scalar, size of epsilons
%            m   scalar, value of mean for random variable
%            a   scalar, value of downside asymmetry, a>1
%Output:     Z   N*1 vector, nodes for Z
%            PI  N*N matrix, transition probabilities
%
%
% NOTE: For an evenly spaced grid of Z:
%        - if 0<p=q<1 then Z follows a standard AR(1) process.
%        - if p=q=0.5 then Z follows has no persistence.
%        - if p>q recessions are longer than expansions (ARCH process).
%
%       Also, note that in the last case, altering p and q does not 
%       preserve the mean of Z (i.e. is not a mean-preserving shift).

% Code by Franz Hamann, Banco de la Republica, December 2013

 PI  = [p 1-p; 1-q q];

 for n = 3:N
    PI = p*[PI zeros(n-1,1); zeros(1,n)] + ...
         (1-p)*[zeros(n-1,1) PI; zeros(1,n)] + ...
         (1-q)*[zeros(1,n); PI zeros(n-1,1)] + ...
         q*[zeros(1,n); zeros(n-1,1) PI];
    PI(2:end-1,:) = PI(2:end-1,:)/2;
 end

 fi = sqrt(N-1)*e;
 Z = linspace(-fi,fi,N)';
 Z = Z + m;
 
 if nargin>5; Z(1) = -a*fi+m; end