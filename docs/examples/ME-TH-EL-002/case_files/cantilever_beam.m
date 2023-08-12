%% Author:Jacob D. R. Bordón
%% ANALYTICAL SOLUTION TO EULER-BERNOULLI CANTILEVER BEAM UNDER BASE LATERAL 
%% UNIT DISPLACEMENT 
clear
clc
close all

%% INPUT DATA
L   = 120;
h   = 10;
b   = 80;
xi  = 0.01;
E   = 19599.92e6;
rho = 2300;
f   = 0.01:0.01:6.9;

%% CALCULATION
E = E*(1+1i*2*xi);
A = b*h;
I = b*h^3/12;
m = rho*A;
nf = length(f);
u = zeros(1,nf);
for kf = 1:nf
  w = 2*pi*f(kf);
  k = m*w^2/(E*I);
  lambda(1) = sqrt(sqrt(k));
  lambda(2) = -lambda(1);
  lambda(3) = sqrt(-sqrt(k));
  lambda(4) = -lambda(3);
  expl = exp(lambda*L);
  A = [1,1,1,1;lambda;lambda.^2.*expl;lambda.^3.*expl];
  b = [1;0;0;0];
  x = linsolve(A,b);
  u(kf) = expl*x;
end

%% EXPORT SOLUTION
% 1st column: frequency [Hz]
% 2nd column: real part of the tip displacement [m]
% 3rd column: imaginary part of the tip displacement [m]
sol = [f;real(u);imag(u)]';
save("-ascii","analytical_solution.txt","sol");
figure
subplot(3,1,1)
semilogy(f,abs(u))
xlabel 'f [Hz]'
ylabel '|u(x=L)| [m]'
subplot(3,1,2)
plot(f,real(u))
xlabel 'f [Hz]'
ylabel 'Re(u(x=L)) [m]'
subplot(3,1,3)
plot(f,imag(u))
xlabel 'f [Hz]'
ylabel 'Im(u(x=L)) [m]'
