% BROCK_MIRMAN Brock and Mirman (1972) optimal stochastic growth model 
 clear all;

%% Parameters
 gamma	= 1;        % coeficiente de aversion al riesgo 2
 beta	= 0.98;     % factor de descuento
 alpha  = 0.3;      % capital share

%% Cadena de Markov de z
 zmean   = 1;       % media
 zstdv   = .1;      % desviacion estandar 0.8 y 0.3
 rho     = 0;       % autocorrelacion
 nz      = 5;       % # estados discretos
 
 [z, PI] = rouwenhorst(nz,zmean,rho,zstdv);

%% Grilla de K = {kmin<k2<k3<...<kmax}
 nk    = 500;
 kmin  = 1.0000e-03;
 kmax  = 0.5;
 k     = linspace(kmin,kmax,nk)';

%% Estado-espacio S = ZxK
 [Z,K] = gridmake(z',k);

%% Solucion exacta
 kex = alpha*beta*zmean*K.^alpha; % evaluada en zmean

%% Conjunto de consumo factible C>0 y funcion de utilidad
 C     = zeros(nk*nz,nk);
 u     = zeros(nk*nz,nk);

 for i=1:nk    
    C(:,i)=Z.*K.^alpha - k(i);   
 end
 if gamma==1;
    u = log(C);
 else
    u = (C.^(1-gamma))./(1-gamma); 
 end
 
 u(C<=0)=-Inf; clear C;

%% Matriz de probabilidad de transicion
 P = kron(speye(nk,nk),repmat(PI,nk,1));

%% Ecuacion de Bellman
 [v,x,pstar] = solvedp(u,P,beta,'policy');  clear P u;

%% Calculo de la distribucion estacionaria
 pie   = markov(pstar);
 figure(1); bar(pie)
 kmean = k(x)'*pie
 cmean = zmean*kmean^alpha-kmean

%% Grafica de comparacion solucion numerica vs exacta
 figure(2); scatter(K,k(x)); hold on; scatter(K,kex); hold off