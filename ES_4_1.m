% Programa que grafica la FIR, FAS y FAP de un AR(p)
% Elaborado por: Juan Tenorio 
% Referencia: Hamilton J. (1994) Time Series Analysis pp. 155
% Curso: Series de Tiempo

clc; clear;

%%%%%%%%%%%%%%%% Proceso AR 
c     = 0;
phi   = [0.01 0.03 0.05];
sigma = 1;
T     = 500;
yinic = [0 0 0];
J     = 30;

p = length(yinic);
w = sigma*randn(T,1);
y = zeros(T,1);
y(1:p) = yinic(:);
for t = p+1:T
  y(t) = c + phi(:)'*y(t-1:-1:t-p) + w(t);
end

%%%%%%%%%%%%%%%% Calculo de la FIR

% Orden de la EeD
p = length(phi);

% Matriz F
F = [ phi(:)' ; ...
      eye(p-1) zeros(p-1,1)];

% Valores propios
L = eig(F);

% Obtencion de los coeficientes c (Proposición 1.2)
p = length(L);
c = zeros(p,1);
for ii=1:p
    den = 1;
    for kk=1:p
        if kk~=ii
            den = den*( L(ii) - L(kk) );
        end
    end
    c(ii) = L(ii)^(p-1)/den;
end

% Calculo de la FIR (funcion impulso respuesta o multiplicador dinamico)
IRF = zeros(J+1,1);
for jj = 0:J
    IRF(jj+1) = c.'*L.^jj;
end


%%%%%%% Graficas

subplot(2, 1, 1); plot( y );
strphi = num2str(phi(1));
for k=2:length(phi)
   strphi = [strphi ', ' num2str(phi(k))];
end
title(['AR(' num2str(p) '): \phi = [ ' num2str(phi) ' ]' ]);
xlim([1 T]);
set(gca, 'FontSize', 8);


subplot(2,1,2); bar( IRF );
xlim([1 J]);
title('FIR');
set(gca, 'FontSize', 8);
xlabel('J (Lag)');
