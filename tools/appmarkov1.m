function [Tran,s,probst]=appmarkov1(const,lambda,meg,sigma,m,N)

% APPMARKOV approximates an AR(1) process with Markov chain
%
% Syntax: Let   y_t = lambda * y_(t-1) + u_t where
%         u_t: Gaussian white noise process with standard deviation sigma
% then
%        [Tran,s,probst,alambda,asigmay] = appmarkov(lambda,sigma,m,N)
% 
%    INPUTS:
%       const   = constant term
%       lambda  = autocorrelation coefficient
%       sigma   = standard deviation of shocks
%       m       = "width" of discretized state space s
%       N       = number of possible states to approximate {y_t} (N=9 usually)
%
%    OUTPUTS:
%       Tran    = transition matrix of the Markov chain
%       s       = discretized state space of y_t
%       probst  = invariant distribution of Markov Chain
%       alambda = theoretical first order autoregression coeff. for Markov chain
%       asigma  = theoretical standard deviation for Markov chain Y_t
%
% NOTE: Inside the code the following options can be set:
%   tol      : convergence tolerance, default is sqrt(eps)
%   maxit    : maximum number of iterations, default is 100


maxit = 100;
tol   = sqrt(eps);

% Discretize the state space

stvy = sqrt(sigma^2/(1-lambda^2)); % standard deviation of y_t

% Equally space grid of points
% ymax=m*vary,ymin=-m*vary, ymax and ymin are two boundary points
ymax = meg+m*stvy;                 % upper boundary of state space
ymin = meg-m*stvy;                 % lower boundary of state space
w=2*m/(N-1)*stvy;                  %length of intervals between g in discretizations
s=meg-m*stvy:w:meg+m*stvy;         % the discretized state space



% Calculate the transition matrix

for i=1:N;
    for j=1:N;
        Tran(i,j)=logncdf(s(j)+w/2,(1-lambda)*(log(meg)-0.5*stvy^2)+lambda*log(s(i)),sigma)-...
            logncdf(s(j)-w/2,(1-lambda)*(log(meg)-0.5*stvy^2)+lambda*log(s(i)),sigma);
    end;
end;


if sum(Tran') ~= ones(1,N)
   str = find(abs(sum(Tran')-ones(1,N))>1e-6);  % find rows not adding up to one
   disp('Error in transition matrix');
   disp(['rows ',num2str(str),' does not sum to one']);
   disp('Normalizing the probabilities');
   d=sum(Tran');
   for i=1:length(str)
       Tran(str(i),:)=(1/d(str(i)))*Tran(str(i),:);
   end    
end


% Calculate the invariant distribution of Markov chain
Trans= Tran';
probst = (1/N)*ones(N,1);            % initial distribution of states

for i=1:maxit
   probst1 = Trans*probst;
   if max(abs(probst1-probst))<tol, break, end;
   probst = probst1; 
end
   

   
   






