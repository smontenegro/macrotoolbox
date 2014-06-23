% AGUIAR-GOPINATH Aguiar & Gopinath (2006) model by state space discretization
% Based on the code by Aguiar & Gopinath, website:
% http://www.economics.harvard.edu/faculty/gopinath/papers_gopinath;
% publication number 8.
% The present code solves the A&G model for trend shocks and exogenous labor
% supply. The notation here is not the same used in the A&G code.

clear,clc

% Model Parameters
mu     = 2;                                    % risk aversion              
beta   = 0.8;                                  % subjective discount factor 
delta  = 0.02;                                 % additional lost of output in autarky
alpha  = 0.32;                                 % capital's share of income
lambda = 0.1;                                  % prob of redemption 0.1
rate   = 0.01;                                 % exogenous interest rate
n1 = 7;                                        % # gridpoints exog. states for the growth rate
n2 = 400;                                      % # gridpoints endo. states for assets (400)
n = n1*n2;                                     % # total states gridpoints
m = n2;                                        % # total action gridpoints

% Construct State-Space
% 1. Approximate 
%   Ln(g_t)= (1-rho)*(Ln(meg)-cg) + rho*Ln(g_t-1) + e_t
% with n1 states Markov chain

rho    = 0.17;                                 % autocorrelation coefficient
sigma  = 0.03;                                 % standard deviation of shocks for growth rate
meg    = 1.006;                                % long run mean for trend income
cg     = .5*((sigma^2)/(1-rho^2));             % AR(1) parameter
cg     = (1-rho)*(log(meg)-cg);                % constant term of the AR(1)
width  = 4.1458;                               % width of state shocks, needs to be larger becasue default ocurrs in extreme states

[prob,g,pi] = appmarkov2(cg,lambda,sigma,width,n1);
g=exp(g);

%[prob,g,pi] = appmarkov1(cg,rho,meg,sigma,width,n1);

% 2. Asset Grid 

amax = 0;                                     % max value of asset grid  
amin = -0.22;                                 % min value of asset grid  
a    = linspace(amin,amax,n2);                % linearly spaced vector
inca = (max(a)-min(a))/(n2-1);                % width of the asset grid
azero= find(a==0);                            % state where assets is zero

[A,G] = gridmake(a',g');                   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration to find the "q" that satisfies the capital market equilibrium %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3. Initial guess for q

q0    = (1/(1+rate))*ones(n,m);
maxit = 1000;

for it=1:maxit
 q0old = q0;    

% 4. Construct the reward function without default and solve

 y  = G./meg;
 y  = repmat(y,1,m);  

 c  = zeros(n,m);
 f  = zeros(n,m);

  for i=1:m
    c(:,i) = y(:,i)+A-(G*a(i)).*q0(:,i); % in detrend form
  end

 f = (c.^(1-mu))/(1-mu); f(c<=0) = NaN;

% 5. Construct the reward function with default 
 cdefault  = zeros(n,m);
 fdefault  = zeros(n,m);

  for i=1:m
    cdefault(:,i) = (1-delta)*y(:,i);
  end

 fdefault = (cdefault.^(1-mu))/(1-mu); fdefault(cdefault<=0) = NaN;

% 6. Comupte new discount factor to account for the detrend form
 betah = (G.^(1-mu))*beta;
 betah = repmat(betah,1,m);

% 7. Solve the system of Bellman eq

 [vbad,xbad,vgood,xgood,default] = solvedpos(f,fdefault,prob,betah,lambda,azero);

% 8. Compute the expected value of default
 Edef = prob*(reshape(default,m,n1)');
 ind = reshape(default,m,n1)';
 ind = sum(ind)==n1; %ind finds a' in which default happens with prob 1 next period

for j=1:m;
    if ind(j)==1;
        Edef(:,j)=1;
    end;
end;

 Edef = kron(Edef,ones(m,1)); %add current a to rows to make same size as q

% 9. Update the price "q"
 q0 = (1/(1+rate))*(1-Edef);
 q0 = max(q0,0);

 if abs(q0old-q0)<10e-6,break,end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of iteration to find the "q"                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10. Optimal Policy Function

policy=(1-default).*xgood+default.*azero;
