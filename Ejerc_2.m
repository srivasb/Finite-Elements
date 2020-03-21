% Este programa calcula una función de aproximación de la función 
%                         f(x)=exp(-x) *x^3.
% La aproximación se construye a través de la combinación lineal de unas
% funciones base Nm, es decir,
%            faprox(x)= a1.N1 + a1.N2 + a3.N3 + ... + aM . NM.
% A diferencia del Ejercicio 1, acá se emplean polinomios de Legendre como
% funciones base. El objetivo del código es encontrar los coeficientes de 
% la combinación lineal que minimizan el residuo de la aproximación. 
%
% Aplicaciones de elementos finitos - I-2020
% Universidad Nacional de Colombia
%
% Profesor: Carlos Galeano
%
close all;clear all; clc;
%% Definición de las variables simbólicas y las funciones 
syms x m;                 % Se definen las variables simbólicas
fi=@(x) exp(-x)*x^3;  % Se define la función que se desea aproximar
Nm=@(x,m) legendreP(m,x); % Se define la familia de funciones base a emplear
M=3;                      % Se define la cantidad de términos a emplear en la aproximación

%% LLenado de la matriz [K] y el vector [F] 
%Se incializa la matriz [K] y el vector [F]
K=zeros(M,M);
F=zeros(M,1);

% Se calcula cada término Kij de la matriz [K] y cada término Fi del vector
% [F]. Este cálculo está basado en el método de residuos ponderados -
% Galerkin
for i=1:M
  for j=1:M
    K(i,j)=int(Nm(x,i-1)*Nm(x,j-1),x,-1,1);  %Se calcula el término (i,j) de la matriz K
  end
  F(i)=int(fi(x)*Nm(x,i-1),x,-1,1);          %Se calcula el i-ésimo término del vector f
end

%% Solución del sistema de ecuaciones
%a=inv(K)*F;          % Invirtiendo la matriz [K]
a=linsolve(K,F);      % Utilizando un solucionador del sistema de ecuaciones

%% Grafica de los resultados obtenidos
np=100;  %Número de puntos con los que se representan las curvas
X=linspace(-1,1,np)';  %Conjunto de valores x en los que se evalúan las funciones
fexacta=zeros(np,1);  %Se inicializa el vector con los valores de la funciòn original
faprox=zeros(np,1);   %Se inicializa el vector con los valores de la función aproximada
R=zeros(np,1);        %Se inicializa el vector que almacena el residuo
%Se calcula el valor de las funciones para cada punto X
for i=1:np
  fexacta(i)=fi(X(i));
  for j=1:M
     faprox(i)=faprox(i)+a(j)*Nm(X(i),j-1);
  end
end
%Se grafica la función original
plot(X,fexacta,'--r','LineWidth',1.3);
grid on; hold on;
%Se grafica la función aproximada
plot(X,faprox,'-b');
xlabel('x');
legend('f(x)','f_{aprox} (x)');
%Se grafica el residuo de la aproximación a lo largo del dominio
figure
plot(X,fexacta-faprox);
grid on;
xlabel('x');
ylabel('R_\Omega')