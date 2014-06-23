function [spath,xpath] = imprespmarkov(P,S,x,T,s)

% IMPRESPMARKOV Impulse response of optimal controls of a Markov Chain 
%
% Usage:
%        [Spath,xpath] = imprespmarkov(P,S,x,T,s)
%
% INPUTS:
%          P     : n by n transition probability matrix
%          S     : n by k state vector 
%          x     : n by 1 policy function (optimal controls)
%          T     : number of simulated time periods
%          s     : initial state index
%
% where n is the total number of states in the Markov Chain
%       k is the total number of state variables
%
% OUTPUTS:
%          spath  : T by 1 path of states conditional on starting at s
%          xpath  : T by 1 path of corresponding optimal controls 
%

 if sum(P')~=1; warning('Not a transition probability matrix'); end;

 spath = zeros(2,T);
 n     = length(P);
 pi    = zeros(n,1); pi(s)=1;        % initial distribution
 
 for t = 1:T;
  spath(:,t) = S'*pi;                % expected state at t
  pi         = (pi'*P)';             % update the distribution pi
 end;
 
 xpath   = x(getindex(spath',S));    % optimal control at t