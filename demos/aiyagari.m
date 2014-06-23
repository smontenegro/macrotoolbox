% AIYAGARI Aiyagari (1993) model
%
% Written by F. Hamann. Feel free to copy, change and distribute

clear all; close all; clc

%% Model Parameters
mu        = 3;                                  % risk aversion              
discrate  = 0.05;                               % discount rate
beta      = 1/(1+discrate);                     % subjective discount factor 
delta     = 0.08;                               % depreciation
A         = 0.80;                               % production technology
alpha     = 0.36;                               % capital's share of income
b         = 3;                                  % ad-hoc debt limit
smean     = 1;                                  % mean of shocks
rho       = 0.2;                                % autocorrelation coefficient
sigma     = 0.4*sqrt(1-rho*rho);                % standard deviation of shocks
tau       = 0.1;                                % tax rate

n1 = 7;                                         % number of nodes for s   
[s, prob] = rouwenhorst(n1,smean,rho,sigma);        
s = exp(s);

disp('Mininimum value of y_i'),s(1)

% Natural debt limit

rate=(0.02:0.02:0.1)';

for i=1:length(s)
    for j=1:length(rate)
        y(i,j)=s(i)*(1-alpha)*(A*(alpha/(rate(j)+delta))^alpha)^(1/(1-alpha));
    end
end

my=y(1,:);

for i=1:length(rate)
    for j=1:length(my)
        limnat(i,j)=my(j)/rate(i);
    end
end

plot(limnat,rate)

% Numeric solution, graphs r & Ea (equilibrium assets)

rate       = [0:0.0025:discrate-0.0025];           % interest rates
amean      = zeros(1,length(rate));
phi        = zeros(1,length(rate));

for i=1:length(rate)
    [amean(i), phi(i)] = Ea(rate(i),beta,mu,delta,A,alpha,b,tau,smean,rho,sigma);
end

figure(1)
plot(amean,rate); hold on;
plot(tau./rate,rate); hold on;
plot(-phi,rate)

% Equilibrium r y Ea

% Bisection method:
% 1.) Define initial r, r0, and compute Ea(r0)
% 2.) Find r1 such that bhat(r1)=Ea(r0)
% 3.) Define r2=(r0+r1)/2 and compute Ea(r2) (bisection)
% 4.) If Ea(r2)>bhat(r2), then replace r0 by r2;
%     if Ea(r2)<bhat(r2), then replace r1 by r2,
%     and use bisection.

r0=0.04 % r0 < discount rate, but close
[am0,phi0]=Ea(r0,beta,mu,delta,A,alpha,b,tau,smean,rho,sigma);
r1=tau/am0;

for i=1:100
    r2=(r0+r1)/2;
    [am2,phi2]=Ea(r2,beta,mu,delta,A,alpha,b,tau,smean,rho,sigma);
    bhat2=tau/r2;i,am2,bhat2
    if abs(am2-bhat2)<1.e-10,break,end
    if am2>bhat2
        r0=r2;
    else am2<bhat2
        r1=r2;
    end 
end    

disp('Equilibrium interest rate'),r2
disp('Equilibrium stock of assets'),am2


% Stationary distribution of assets at equilibrium rate

[Eae,phi,pi]=Ea(r2,beta,mu,delta,A,alpha,b,tau,smean,rho,sigma);

hold off
figure(2)
bar(pi)

