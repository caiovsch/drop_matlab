close all; clear variables; clc;

tic

% Diretorio de trabalho
cd 'C:\Users\caios\Desktop\PG\Matlab\periodico';


N = 100;                          % número de pontos (espacial)
M = 200;                         % número de pontos (temporal)
r_inicial = 0;                   % r inicial
r_final = 1;                     % r final
dr = (r_final-r_inicial)/(N-1);  % passo em r
r = (r_inicial:dr:r_final)';     % vetor de r
dt = 0.025;                      % passo em t 
sigma = -dt/(12*dr^2);           % parâmetro auxiliar
i = (1:N)';                      % vetor dos passos espaciais
f = 1;                           % fator pra desligar um termo da equação

h0 = 1-r.^2;                     % perfil da superfícies livre em t = 0
h_anterior = h0;                 % primeiro h para entrar no loop

fid = fopen('perfil_0000.dat','w');
for k = 1:numel(r)
    fprintf(fid, '%f\t%f\n', r(k), h0(k));
end
fclose(fid);

fid = fopen('tempo.dat', 'w');
for k = 0:M
   fprintf(fid, '%d\t%d\n', k, k*dt); 
end
fclose(fid);

% Construção da matriz B_j
alpha = sigma*(1-f./(2*(i(2:N-1)-1)));
beta = -2*sigma;
gamma = sigma*(1+f./(2*(i(2:N-1)-1)));

BdiagS = diag([2*sigma; gamma], 1);
BdiagP = diag(ones(N,1)*beta);
BdiagI = diag([alpha; 2*sigma], -1);
B = BdiagS + BdiagP + BdiagI;


for j =1:M

% Construção da matriz A_j
a = 2*sigma*(1-f./(2*(i(2:N-1)-1))).*h_anterior(1:N-2).^3;
b = 1-4*sigma*h_anterior(1:N).^3;
c = 2*sigma*(1+f./(2*(i(2:N-1)-1))).*h_anterior(3:N).^3;

AdiagS = diag([4*sigma*h_anterior(2)^3; c], 1);     % diagonal superior
AdiagP = diag(b);                                   % diagonal principal    
AdiagI = diag([a; 4*sigma*h_anterior(N-1)^3], -1);  % diagonal inferior
A = AdiagS+AdiagP+AdiagI;                           % matriz A

% Construção do vetor delta_j
delta = h_anterior;

% Cálculo no passo temporal seguinte H_(j+1)
s = A^-1*(B*h_anterior.^4 + delta);
h_anterior = s;              % redefine o h para o próximo loop

if j<10
    filename = ['perfil_000', num2str(j), '.dat'];
elseif j<100
    filename = ['perfil_00', num2str(j), '.dat'];
elseif j<1000
    filename = ['perfil_0', num2str(j), '.dat'];
end
fid = fopen(filename,'w');
for k = 1:numel(r)
    fprintf(fid, '%d\t%d\n', r(k), h_anterior(k));
end
fclose(fid);

end


h_final = h_anterior;

figure
plot(r, h0, 'b')
hold on
plot(r, h_final, 'r*')
title('Perfil da superficie livre')
legend('Condicao inicial', 't final')
xlabel('r'); ylabel('h(r,t)')
hold off

trapezio_inicial = trapz(r, r.*h0);
trapezio_final = trapz(r, r.*h_final);


toc
