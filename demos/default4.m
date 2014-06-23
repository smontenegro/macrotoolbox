% DEFAULT4 Solves sovereign's default problem (with fixed reentry cost)
clear all

% Parameters
beta    = 0.97;             % Discount factor 0.97, literature 0.96-0.98
sigma   = 2.0;              % Coeff risk aversion 2.0, literature 2-5
R       = 1.0295;           % Riskfree rate 1.02
zmean   = 0.0;              % Mean of output level 0
rho     = 0.8;              % autocorrelation coefficient 0.9
stddev  = 0.02;             % std deviation of shocks 0.125
tau     = 0.5;              % permanent output loss if default (%) 0.7 
ec      = 0.3*exp(zmean);   % one-time entry cost (in units of c) 3
xc      = 2*ec;             % one-time default cost (in units of c) 6

% Build the state-space grid
% 1. Approximate output-shocks with nz states Markov chain
nz   = 2;                       % number of skill levels
[z, prob] = rouwenhorst(nz,zmean,rho,stddev);
pzs = markov(prob);
y = exp(z)
% 2. Approximate aseets with na discrete states
na   = 500;                     % number of asset grid points
amax = 0;                       % max value of debt grid 
amin = -2*exp(zmean);
%amin = -min((1-tau)*y)/(R-1);   % natural borrowing constraint
a    = linspace(amin,amax,na);  % linearly spaced vector
d    = [0 1];                   % in default (d = 0); in repayment (d=1)
nd   = length(d);                       

[Y,A,D] = gridmake(y',a',d');
[aa,dd] = gridmake(a',d');         % action-space grid

n = nz*na*nd;                        % number state gridpoints
m = na*nd;                           % number action gridpoints  

% Construct the reward function 
c  = zeros(n,m);
  for i=1:m
      c(:,i) = ((1-tau)*Y).*(1-D)*(1-dd(i))+ ...
               ((1-tau)*Y+R*A-aa(i)-ec).*(1-D)*dd(i)+ ...
               (Y+R*A-aa(i)).*D*dd(i) + ...
               (Y+R*A-aa(i)-xc).*D*(1-dd(i));
   end

u = (c.^(1-sigma)-1)/(1-sigma); u(c<=0) = NaN;  

% Construct the transition matrix
I  = speye(m,m);                            % sparse identity
P  = kron(I,repmat(prob,m,1));              % stack moving PI diagonally

% Solve the Model  
[v,x,pstar] = solvedp(u,P,beta,'policy');  clear P;      
anext = aa(x);
%pie = ergdist(pstar);
[pie,F] = markov(pstar);

% Calculations
ameanr = A(na*nz+1:end)'*pie(na*nz+1:end)
a2yr   = ameanr/(y*pzs)

% Reshaping
vr     = reshape(v,nz*na,nd);
anextr = reshape(anext,nz*na,nd);
pier   = reshape(pie,nz*na,nd);

vd = vr(:,1); anextd = anextr(:,1); pied = pier(:,1);
vr = vr(:,2); anextr = anextr(:,2); pief = pier(:,2);

vdr = reshape(vd,nz,na)'; anextdr = reshape(anextd,nz,na)';
vrr = reshape(vr,nz,na)'; anextrr = reshape(anextr,nz,na)';

pieir = reshape(pie(1:na*nz)',nz,na)';
piefr = reshape(pie(na*nz+1:end)',nz,na)';

Fi=cumsum(pieir); 
Ff=cumsum(piefr);

