function [vbad,xbad,vgood,xgood,default] = solvedpos(f,fd,prob,beta,beta2,lambda,azero)

% SOLVEDP  Solves stochastic dynamic program by state-space discretization
% Usage:
%         [v,x,pstar] = solvedpos(f,P,beta,alg,v);

% Syntax: Let m = # actions, n = # of ALL possible states, then
%   INPUTS
%      f       : n by m reward function
%      P       : mn by n transition probability matrix
%      alg     : algorithm used = 'policy' for policy function iteration
%                               = 'value'  for value fucntion iteration
%      v       : initial guess for value function (default is zero) 
%    OUTPUTS
%      v       : value function, v(s)
%      x       : optimal controls, x(s)
%      pstar   : optimal transition probability for states s, P(x)


% SET CONVERGENCE PARAMETER DEFAULTS
  maxit    = 10000;                    % maximum number of iterations
  tol      = 10e-6;                    % convergence tolerance, usually tol=sqrt(eps)
  prtiters = 1;                        % print iterations (1) or not (0)

[n,m] = size(f);
vbad=ones(n,1);
vgood=ones(n,1);
vbg = reshape(vgood,m,n/m);
vbg = vbg(azero,:);
vbg = repmat(vbg,m,1);
vbg = reshape(vbg,n,1);

    if prtiters, disp('Solving Bellman equation by Value function iteration'); end
    for it=1:maxit
      vbadold=vbad;
      vgoodold=vgood;
      Evgood = prob*(reshape(vgood,m,n/m)');Evgood = kron(Evgood,ones(m,1));
      Evbad = prob*(reshape(vbad,m,n/m)');Evbad = kron(Evbad,ones(m,1));
      Evbg = prob*(reshape(vbg,m,n/m)');Evbg = kron(Evbg,ones(m,1));
      [vbad,xbad]   = max(fd+lambda*(beta.*Evbg)+(1-lambda)*(beta.*Evbad),[],2);
      [vgood,xgood] = max(f+beta.*Evgood,[],2);
      default=vbad>vgood | isnan(vgood)==1;
      vgood(find(vbad>vgood | isnan(vgood)==1))=vbad(find(vbad>vgood| isnan(vgood)==1));
      vbg = reshape(vgood,m,n/m);
      vbg = vbg(azero,:);
      vbg = repmat(vbg,m,1);
      vbg = reshape(vbg,n,1);
      v=[vbad vgood];vold=[vbadold vgoodold];
      change = norm(v-vold);                % compute change
      if prtiters 
         fprintf ('%5i %10.1e\n',it,change) % print progress
      end
      if change<tol, break, end;            % convergence check
    end
    if change>tol, warning('Failure to converge in solvedp'), end;

    

