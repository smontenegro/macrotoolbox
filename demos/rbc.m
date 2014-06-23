% RBC.M Solves standard real business cycle model with inelastic labor
%       supply. 
%
% See any macro textbook.
%
% Written by F. Hamann. Feel free to copy, change and distribute

% Model parameters
  beta  =  0.98;                % discount factor
  mu    =  2.0;                 % utility parameter
  alpha =  0.3;                 % production parameter
  delta =  0.1;                 % depreciation rate
  sigma =  0.35;                 % productivity shock volatility
  rho   = 0.5;                  % productivity shock persistence

  n1 =  100;                    % number of endogenous states
  n2 =  9;                      % number of exogenous states
  m  =  n1;                     % number of actions
  n  = n1*n2;
  
% Approximate productivity shocks with n2 states Markov chain
  [a, prob] = rouwenhorst(n2,0,rho,sigma)
  a = exp(a);
  
% Construct state space
  kmin  = 0.1;                       % minimum state
  kmax  = 20;                      % maximum state
  k = linspace(kmin,kmax,n1);
  [A,K,Knext] = gridmake(a',k',k');
  
% Construct the reward function, f
  c = A.*K.^alpha + (1-delta)*K - Knext;
  f = c.^(1-mu)./(1-mu);
  i = find(c<0);
  f(i) = NaN;
  f = reshape(f,n,m);
  
% Construct the transition matrix
  PI = repmat(prob,n1,1);                         % stacks [prob;prob; ... prob] = PI
  I  = speye(n1,n1);                              % sparse identity (avoid mem troubles)
  P  = kron(I,PI);                                % stack moving PI one place columnwise

% Solve the DP problem using policy iteration
  [v,x,Popt] = solvedp(f,P,beta,'policy');

% Solve infinite-horizon model via function iteration
  %[v,x,pstar] = solvedp(f,P,delta,'value');
 
 piss = ergdist(Popt);

  v = reshape(v,n2,n1);
  figure(1); surf(k,a,v);
  xlabel('k'); ylabel('a'); zlabel('v');

  figure(2); bar(piss);
  title('Ergodic density of K')
