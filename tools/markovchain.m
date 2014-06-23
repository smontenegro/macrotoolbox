function [X,Q] = markovchain(N,p,q,e,m,a)

% MARKOVCHAIN Markov process and transition matrix 
% 
% Purpose:   Follows Rouwenhorst to find a Markov chain whose sample paths
%            approximate a family of random walk processes with mean m, std
%            dev e and serial correlation p+q-1. Depending on p and q, the
%             process could be AR(1) or ARCH.
%            See: Cooley's RBC book, chapter 10 by K. G. Rouwenhorst
%
% Usage:    [X, Q] = markovchain(N,p,q,e)
%
% Input:     N   scalar, number of nodes for X, usually odd number
%            p   scalar, governs cond. prob of moving up after low state
%            q   scalar, governs cond. prob of moving down after low state
%            e   scalar, size of epsilons
%            m   scalar, value of mean for random variable
%            a   scalar, value of downside asymmetry, a>1
% Output:    X   N*1 vector, nodes for X
%            Q   N*N matrix, transition probabilities
%
%
% NOTE: For an evenly spaced grid of X:
%        - if 0<p=q<1 then X follows a standard AR(1) process.
%        - if p=q=0.5 then X follows has no persistence.
%        - if p>q bad times are longer than good ones (ARCH process).
%
%       Also, note that in the last case, altering p and q does not 
%       preserve the mean of X (i.e. is not a mean-preserving shift).

% Code by Franz Hamann, Banco de la Republica, December 2013

 Q  = [p 1-p; 1-q q];

 for n = 3:N
    Q = p*[Q zeros(n-1,1); zeros(1,n)] + ...
         (1-p)*[zeros(n-1,1) Q; zeros(1,n)] + ...
         (1-q)*[zeros(1,n); Q zeros(n-1,1)] + ...
         q*[zeros(1,n); zeros(n-1,1) Q];
    Q(2:end-1,:) = Q(2:end-1,:)/2;
 end

 fi = sqrt(N-1)*e;
 X = linspace(-fi,fi,N)';
 X = X + m;
 
 if nargin>5; X(1) = -a*fi+m; end