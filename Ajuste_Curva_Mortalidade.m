% =========================================================
% AJUSTE POLINOMIAL DAS TAXAS DE MORTALIDADE (D. saccharalis)
% Estagio 1 e Estagio 2 (Temperatura), Estagio 3 (Temperatura x Senescencia)
% =========================================================
clear all
clc
close all

% 1. Dados de Entrada (Variavel Independente: Temperatura em oC)
T = [18, 22, 25, 28, 32];

% 2. Coordenadas Y (Taxas de mortalidade calculadas)
d1 = [0.01005, 0.00054, 0.00711, 0.00354, 0.01705];  % Estagio 1 (Ovos/L1 a L3)
d2 = [0.01177, 0.00293, 0.00484, 0.00373, 0.01178];  % Estagio 2 (L4 a L7/Pupas)
d3 = [0.37594, 0.17452, 0.17986, 0.19048, 0.32895];  % Estagio 3 (Adultos —> d = 1 / tau)

% 3. Ajuste dos Polinomios de Grau 2 (polyfit)
p1 = polyfit(T, d1, 2);
p2 = polyfit(T, d2, 2);
p3 = polyfit(T, d3, 2);

% Exibindo os parametros no ecra (Command Window)
fprintf(‘— Parametros do Estagio 1 —\n’);
fprintf(‘a1 = %e, b1 = %e, c1 = %e\n\n’, p1(1), p1(2), p1(3));

fprintf(‘— Parametros do Estagio 2 —\n’);
fprintf(‘a2 = %e, b2 = %e, c2 = %e\n\n’, p2(1), p2(2), p2(3));

fprintf(‘— Parametros do Estagio Adulto (Base Termica) —\n’);
fprintf(‘a3 = %e, b3 = %e, c3 = %e\n\n’, p3(1), p3(2), p3(3));

% 4. Geracao de pontos para suavizar as curvas no grafico
T_curva = linspace(min(T), max(T), 100);
d1_curva = polyval(p1, T_curva);
d2_curva = polyval(p2, T_curva);
d3_curva = polyval(p3, T_curva);
view(3)
