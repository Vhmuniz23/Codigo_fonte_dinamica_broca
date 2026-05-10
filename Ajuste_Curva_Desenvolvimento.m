% =========================================================
% AJUSTE DO MODELO DE BRIÈRE-1 (D. saccharalis)
% Estágio 1 e Estágio 2 (Taxas de Desenvolvimento)
% =========================================================
clear all
clc
close all

% 1. Dados de Entrada (Temperatura em °C)
T_dados = [18, 22, 25, 28, 32];

% 2. Taxas de Desenvolvimento Calculadas (R = 1/tau)
% Valores corrigidos conforme o apêndice metodológico
R1_dados = [0.02039, 0.03833, 0.05414, 0.06382, 0.06557]; % Estágio 1
R2_dados = [0.00954, 0.01972, 0.02829, 0.03269, 0.03643]; % Estágio 2

% 3. Definição do Modelo de Brière-1
% Equação: R(T) = a * T * (T - Tmin) * sqrt(Tmax - T)
briere1 = @(p, T) (p(1) * T .* (T - p(2)) .* sqrt(max(0, p(3) - T))) .* (T > p(2) & T < p(3));

% 4. Ajuste Não-Linear (lsqcurvefit)
% Parâmetros iniciais: [a, Tmin, Tmax]
chute_inicial = [1e-4, 10, 40];
opcoes = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-10);
lb = [0, 0, 32];  % Limites inferiores (Tmax deve ser pelo menos a maior temp testada)
ub = [1, 15, 50]; % Limites superiores

% Ajuste para o Estágio 1
p1 = lsqcurvefit(briere1, chute_inicial, T_dados, R1_dados, lb, ub, opcoes);
% Ajuste para o Estágio 2
p2 = lsqcurvefit(briere1, chute_inicial, T_dados, R2_dados, lb, ub, opcoes);

% 5. Exibição dos Resultados
fprintf('--- Parâmetros Brière-1: Estágio 1 ---\n');
fprintf('a = %e\nTo (Tmin) = %.2f\nTL (Tmax) = %.2f\n\n', p1(1), p1(2), p1(3));

fprintf('--- Parâmetros Brière-1: Estágio 2 ---\n');
fprintf('a = %e\nTo (Tmin) = %.2f\nTL (Tmax) = %.2f\n\n', p2(1), p2(2), p2(3));

% 6. Geração de Curvas e Gráficos
T_ajuste = linspace(5, 45, 200);
R1_ajuste = briere1(p1, T_ajuste);
R2_ajuste = briere1(p2, T_ajuste);

figure('Name', 'Ajuste de Desenvolvimento - Brière-1', 'Color', 'w');
hold on;

% Plot Estágio 1
plot(T_dados, R1_dados, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Dados Estágio 1');
plot(T_ajuste, R1_ajuste, 'b-', 'LineWidth', 2, 'DisplayName', 'Brière-1 Estágio 1');

% Plot Estágio 2
plot(T_dados, R2_dados, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Dados Estágio 2');
plot(T_ajuste, R2_ajuste, 'r-', 'LineWidth', 2, 'DisplayName', 'Brière-1 Estágio 2');

grid on;
xlabel('Temperatura (°C)', 'FontSize', 12);
ylabel('Taxa de Desenvolvimento (dia^{-1})', 'FontSize', 12);
title('Ajuste das Taxas de Desenvolvimento (D. saccharalis)', 'FontSize', 14);
legend('Location', 'northwest');