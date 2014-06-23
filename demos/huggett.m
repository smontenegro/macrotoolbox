% HUGGETT.M  PARTIAL EQUILIBRIUM version of Huggett (1993)
%            See Huggett (1993) The importance of the risk-free rate
%            in heterogeneous-agent incomplete insurance economies
%            Journal of Economic Dynamics and Control 17 953-969
%
% Written by F. Hamann. Feel free to copy, change and distribute.

sigma = 1.5;
beta  = 0.993;
p1	  = 0.5;
p2    = 0.8;
PI	  = [1-p1 p1;1-p2 p2];
em	  = 0.1;
es	  = 1.0;
q     = 0.997;

e     = [em;es];
ne    = length(e);
amin  = -8;
amax  = 3;
na    = 250;
a     = linspace(amin,amax,na)';
C     = zeros(na*ne,na);
u     = zeros(na*ne,na);

[E,A] = gridmake(e,a);
anext    = a;

for i=1:na
    C(:,i) = E+A-q*anext(i);
end

j=find(C<=0); C(j)=NaN;
u = (C.^(1-sigma)-1)./(1-sigma); u(j)=-Inf;
clear j C;

I = speye(na,na);
P = kron(I,repmat(PI,na,1));

[v,x,pstar] = solvedp(u,P,beta,'policy');
clear P u;

ax = a(x);
pi    = markov(pstar);
F     = cumsum(pi);

cons = E+A-q.*ax;
amean = pi'*A;
cmean = pi'*cons;

plotdp(v,x,pstar,E,A,e,a)


%% Summary statistics

emean = E'*pi;
sd.c  = sqrt(((E+A-q*ax-cmean).^2)'*pi);
sd.e  = sqrt(((E-emean).^2)'*pi);


fprintf('\nModel Statistics \n')

rowname = ['y mean       ';
           'c mean       ';
           'b mean       ';
           'b/y          ';
           'Y stdv       ';
           'C stdv       ';
           'Csdv/Ysdv    '];
         

values=[emean;cmean;amean;amean/emean;sd.e;sd.c;sd.c/sd.e];

[rowname num2str(values)]