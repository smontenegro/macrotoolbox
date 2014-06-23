%% ENDOWMENT.M  Small Open Endowment Economy with Incomplete Markets
%
% Written by F. Hamann. Feel free to copy, change and distribute.

%% Parametros

sigma	= 2;         % risk aversion
beta	= 0.985;     % discount factor
R       = 1.0125;    % gross asset return rate (vs. R=1)

if R*beta>=1; display('Set beta*R<1 for convergence'); break; end;

%% Markov chain for y

ymean   = 1;       % meean
ystdv   = .3;      % std.dev. [0.7 0.3] if too low/high, pdf hits limit
rho     = .1;      % autocorrelation
ny      = 2;       % # of discrete states

[y, prob] = rouwenhorst(ny,ymean,rho,ystdv)

%% B = {bmin<b2<b3<...<bmax}

nb    = 500;
bmin  = -15;
bmax  =  10;
b     = linspace(bmin,bmax,nb)';

%% State-space S = YxB

[Y,B] = gridmake(y',b);
bb    = b;

%% Utility function and feasible consumption C>=0

C     = zeros(nb*ny,nb);
u     = zeros(nb*ny,nb);

for i=1:nb    
    C(:,i)=Y+R*B-bb(i);   % Y nonstorable, B noncontingent
end

j    = find(C<=0); C(j) = NaN;
u    = (C.^(1-sigma)-1)./(1-sigma); u(j)=-Inf;
clear j C;

%% Transition probability matrix (see Sargent and Ljundqvist)

I = speye(nb,nb);
P = kron(I,repmat(prob,nb,1));

%% Bellman equation

[v,x,pstar] = solvedp(u,P,beta,'policy'); 
bnext       = b(x);
clear P u;

%% Steady State Distribution

pie   = markov(pstar);
F     = cumsum(pie);

cons  = Y+R*B-bnext;

bmean = pie'*bnext;
cmean = pie'*cons;

%% Model properties

%plotdp(v,x,pstar,Y,B,y,b);


%% Simulation

T      = 500;      
s0     = findnearest(bmean,B);   
spath  = simulmarkov(pstar,T,s0);

ypath  = Y(spath);
cpath  = ypath + R*B(spath)-b(x(spath));
CApath = b(x(spath))-B(spath);

sd.y = std(ypath);
sd.c = std(cpath);

figure(2)
plot([ypath cpath])

[sdev,corrcont,corr,acov] = samplemoms([ypath cpath CApath],1,3)

%% Model steady state statistics
display(' ')
display('Steady State Statistics')

rowname = ['y mean       ';
           'c mean       ';
           'b mean       ';
           'b/y          ';
           'Y stdv       ';
           'C stdv       ';
           'Csdv/Ysdv    '];
         

values=[ymean;cmean;bmean;bmean/ymean;sd.y;sd.c;sd.c/sd.y];

[rowname num2str(values)]

