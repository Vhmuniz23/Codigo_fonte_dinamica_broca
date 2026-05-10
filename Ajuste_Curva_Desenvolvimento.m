% =========================================================
% AJUSTE DO MODELO DE BRIERE-1 (D. saccharalis)
% Estagio 1 e Estagio 2 (Taxas de Desenvolvimento)
% =========================================================
clear all
clc
close all

% 1. Dados de Entrada (Temperatura em C)
T_dados = [18, 22, 25, 28, 32];

% 2. Taxas de Desenvolvimento Calculadas (R = 1/tau)
% Valores corrigidos conforme o apendice metodologico
R1_dados = [0.02039, 0.03833, 0.05414, 0.06382, 0.06557]; % Estagio 1
R2_dados = [0.00954, 0.01972, 0.02829, 0.03269, 0.03643]; % Estagio 2

% 3. Definicao do Modelo de Briere-1
% Equacao: R(T) = a * T * (T - Tmin) * sqrt(Tmax - T)
briere1 = @(p, T) (p(1) * T .* (T - p(2)) .* sqrt(max(0, p(3) - T))) .* (T > p(2) & T < p(3));

% 4. Ajuste Nao-Linear (lsqcurvefit)
% Parametros iniciais: [a, Tmin, Tmax]
chute_inicial = [1e-4, 10, 40];
opcoes = optimoptions(‘lsqcurvefit’, ‘Display’, ‘off’, ‘FunctionTolerance’, 1e-10);
lb = [0, 0, 32];  % Limites inferiores (Tmax deve ser pelo menos a maior temp testada)
ub = [1, 15, 50]; % Limites superiores

% Ajuste para o Estagio 1
p1 = lsqcurvefit(briere1, chute_inicial, T_dados, R1_dados, lb, ub, opcoes);
% Ajuste para o Estagio 2
p2 = lsqcurvefit(briere1, chute_inicial, T_dados, R2_dados, lb, ub, opcoes);

% 5. Exibicao dos Resultados
fprintf(‘— Parametros Briere-1: Estagio 1 —\n’);
fprintf(‘a = %e\nTo (Tmin) = %.2f\nTL (Tmax) = %.2f\n\n’, p1(1), p1(2), p1(3));
legend(‘Location’, ‘northwest’);
