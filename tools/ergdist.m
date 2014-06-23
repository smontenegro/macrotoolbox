function pi = ergdist(P,pi)

% ERGDIST Ergodic distribution of a n-state Markov Chain 
%
% Usage:
%            pi = ergdist(P,pi)
%
% INPUTS:
%          P   : n by n transition probability matrix
%          pi  : n by 1 initial guess of invariant distribution
%
% where n is the total number of states in the Markov Chain.
%
% OUTPUTS:
%          pi  : n by 1 ergodic distribution
%
% NOTE: Does not work for Markov Chains with transient states.  
%       Follows Ljundqvist and Sargent (2000). 

if sum(P')~=1; warning('Not a transition probability matrix'); end;

n = length(P);

if nargin<2
  pi = (1/n)*ones(n,1);  % Start with uniform dist., if no seed provided.
end

maxit = 500;
tol   = sqrt(eps);

for it = 1:maxit
  pinext = (pi'*P)';
  change = max(abs(pinext-pi));
  %fprintf ('%5i %10.1e\n',it,change) 
  if change<tol; break; end;
  pi  = pinext;
end  
  
