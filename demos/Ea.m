function [amean, phi,pi] = Ea(rate,beta,mu,delta,A,alpha,b,tau,smean,rho,sigma)

% Ea.m Aiyagari's asset demand function under incomplete markets

% Set number of actions and states 
% NOTE: follows Fackler & Miranda's notation!
n1    = 7;                                      % # gridpoints exog. states
n2    = 100;                                    % # gridpoints endo. states
n     = n1*n2;                                  % # total states gridpoints
m     = n2;                                     % # total action gridpoints

% Construct State-Space

% 1. Approximate shocks with n1 states Markov chain
  [s, prob] = rouwenhorst(n1,smean,rho,sigma)
  s = exp(s);

% 2. Asset Grid (and setting the borrowing limit)
  wage = (1-alpha)*(A*(alpha/(rate+delta))^alpha)^(1/(1-alpha));

  if rate<=0
    phi = b;                                    % ad-hoc borrowing limit
  else
    phi  = min(b, wage*s(1)/rate);              % natural borrowing limit
  end

  amax = 40;                                    % max value of debt grid  
  amin = -phi;                                  % borrowing constraint
  a    = linspace(amin,amax,n2);                % linearly spaced vector
  inca = (max(a)-min(a))/(n2-1);                % width of the asset grid

[S,A] = gridmake(s',a');

% Construct the reward function, f
c  = zeros(n,m);
f  = zeros(n,m);

  for i=1:m
    c(:,i) = (wage*S+(1+rate)*A-a(i)-tau);
  end

j = find(c<=0);
f = (c.^(1-mu)-1)/(1-mu); f(j) = NaN;

% Construct the transition matrix
PI = repmat(prob,n2,1);                         
I  = speye(n2,n2);                              % sparse identity
P  = kron(I,PI);                                % stack moving PI diagonally


% Solve the Model
[v,x,pstar] = solvedp(f,P,beta,'policy');         % can also use 'value'


% Analysis of the Model
  ax = a(x);                                      % Store policy in ax vector
  v = reshape(v,n1,n2)';                          % Rearrange v in a matrix
  x = reshape(x,n1,n2)';                          % Rearrange x in a matrix


% Ergodic Distribution
  pi = ergdist(pstar);
  amean =ax*pi;
 
