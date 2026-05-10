% =========================================================
% AJUSTE POLINOMIAL DAS TAXAS DE MORTALIDADE (D. saccharalis)
% Estágio 1 e Estágio 2 (Temperatura), Estágio 3 (Temperatura x Senescência)
% =========================================================
clear all
clc
close all

% 1. Dados de Entrada (Variável Independente: Temperatura em °C)
T = [18, 22, 25, 28, 32];

% 2. Coordenadas Y (Taxas de mortalidade calculadas e corrigidas)
d1 = [0.01005, 0.00054, 0.00711, 0.00354, 0.01705];  % Estágio 1 (Ovos/L1 a L3)
d2 = [0.01177, 0.00293, 0.00484, 0.00373, 0.01178];  % Estágio 2 (L4 a L7/Pupas)
d3 = [0.37594, 0.17452, 0.17986, 0.19048, 0.32895];  % Estágio 3 (Adultos --> d = 1 / tau)

% 3. Ajuste dos Polinômios de Grau 2 (polyfit)
p1 = polyfit(T, d1, 2);
p2 = polyfit(T, d2, 2);
p3 = polyfit(T, d3, 2);

% Exibindo os parâmetros no ecrã (Command Window)
fprintf('--- Parâmetros do Estágio 1 ---\n');
fprintf('a1 = %e, b1 = %e, c1 = %e\n\n', p1(1), p1(2), p1(3));

fprintf('--- Parâmetros do Estágio 2 ---\n');
fprintf('a2 = %e, b2 = %e, c2 = %e\n\n', p2(1), p2(2), p2(3));

fprintf('--- Parâmetros do Estágio Adulto (Base Térmica) ---\n');
fprintf('a3 = %e, b3 = %e, c3 = %e\n\n', p3(1), p3(2), p3(3));

% 4. Geração de pontos para suavizar as curvas no gráfico
T_curva = linspace(min(T), max(T), 100);
d1_curva = polyval(p1, T_curva);
d2_curva = polyval(p2, T_curva);
d3_curva = polyval(p3, T_curva);

% =============================================================
% 5. AJUSTE DO TERMO DE SENESCÊNCIA DO ESTÁGIO 3 - ADULTOS (d3)
% =============================================================
A3 = 10;     % Idade máxima atingível em dias
T_alvo = 25; % Exemplo: Temperatura de estudo para a senescência

% Extração automática de a3, b3 e c3 a partir do polinômio ajustado acima
a3 = p3(1);
b3 = p3(2);
c3 = p3(3);

% ---> DADOS OBSERVADOS
idades_observadas = [1, 3, 5, 7, 9]; 
mortalidade_obs = [0.05, 0.15, 0.40, 1.20, 5.00];

% Definição da função matemática de senescência
func_senescencia = @(coef, a) (a3*(T_alvo^2) + b3*T_alvo + c3) .* coef(1) .* a ./ (A3 - a).^2;

% Encontrando o parâmetro d3
chute_d3 = 0.1;
opcoes = optimoptions('lsqcurvefit', 'Display', 'off');
d3_estimado = lsqcurvefit(func_senescencia, chute_d3, idades_observadas, mortalidade_obs, 0, 100, opcoes);
fprintf('--- Parâmetro de Senescência (Adultos) ---\n');
fprintf('d3 estimado = %e\n\n', d3_estimado);

% =========================================================
% 6. SUPERFÍCIE 3D: MORTALIDADE TÉRMICA + ETÁRIA
% =========================================================
% Malha de temperatura e idade
T_3D = linspace(min(T), max(T), 80);
a_3D = linspace(0.1, 9.5, 80);
[TT, AA] = meshgrid(T_3D, a_3D);

% Mortalidade térmica dos adultos
mu_termica = a3.*TT.^2 + b3.*TT + c3;

% Mortalidade etária (senescência)
mu_etaria = d3_estimado .* (AA ./ (A3 - AA)).^2;

% Mortalidade total adulto dependente de T e idade
mu_total = mu_termica .* mu_etaria;

% =========================================================
% 7. PLOTAGEM GRÁFICA
% =========================================================
figure('Name', 'Taxas de Mortalidade vs Temperatura', 'Color', 'w', 'Position', [50, 50, 1400, 700]);

% ---------------------------------------------------------
% Subplot 1: Imaturos
% ---------------------------------------------------------
subplot(2,2,1);
plot(T, d1, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
plot(T_curva, d1_curva, 'b-', 'LineWidth', 2);
plot(T, d2, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(T_curva, d2_curva, 'r-', 'LineWidth', 2);
grid on;
xlabel('Temperatura T (°C)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mortalidade Térmica \mu_{1,2}(T)', 'FontSize', 11, 'FontWeight', 'bold');
title('Ajuste: Estágios Imaturos - Temperatura', 'FontSize', 13);
legend('Estágio 1 (Dados)', 'Estágio 1 (Ajuste)', 'Estágio 2 (Dados)', 'Estágio 2 (Ajuste)', 'Location', 'north', 'FontSize', 8);

% ---------------------------------------------------------
% Subplot 2: Adultos térmico
% ---------------------------------------------------------
subplot(2,2,2);
plot(T, d3, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
plot(T_curva, d3_curva, 'k-', 'LineWidth', 2);
grid on;
xlabel('Temperatura T (°C)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mortalidade Térmica \mu_{3}(T)', 'FontSize', 11, 'FontWeight', 'bold');
title('Ajuste: Estágio Adulto - Temperatura', 'FontSize', 13);
legend('Estágio 3 (Dados)', 'Estágio 3 (Ajuste)', 'Location', 'northwest', 'FontSize', 8);

% ---------------------------------------------------------
% Subplot 3: Senescência
% ---------------------------------------------------------
subplot(2,2,3);
plot(idades_observadas, mortalidade_obs, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); hold on;
a_curva = linspace(0.1, 9.5, 100);
mu_senes_curva = func_senescencia(d3_estimado, a_curva);
plot(a_curva, mu_senes_curva, 'g-', 'LineWidth', 2);
grid on;
xlabel('Idade Fisiológica do Adulto a (dias)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Mortalidade Etária \mu_3(a)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Estimativa: Estágio Adulto - Senescência (T = %d °C)', T_alvo), 'FontSize', 13);
legend('Dados Etários', 'Curva de Senescência', 'Location', 'northwest', 'FontSize', 8);

% ---------------------------------------------------------
% Subplot 4: Superfície 3D
% ---------------------------------------------------------
subplot(2,2,4);
surf(TT, AA, mu_total);
shading interp
colorbar
colormap(jet)
grid on
xlabel('Temperatura T (°C)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Idade Fisiológica a (dias)', 'FontSize', 11, 'FontWeight', 'bold');
zlabel('Mortalidade Total \mu_3(T,a)', 'FontSize', 11, 'FontWeight', 'bold');
title('Mortalidade Adulta Térmica x Etária', 'FontSize', 13);
view(3)