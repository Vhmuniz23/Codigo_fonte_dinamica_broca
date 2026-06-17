% =========================================================================
% MODELO POPULACIONAL ESTRUTURADO EM IDADE (EDP) - DIATRAEA SACCHARALIS
% =========================================================================
% TIPO DE MODELO: 
% Sistema de Equações Diferenciais Parciais (EDP) de Transporte 
% (Baseado na equação de McKendrick-von Foerster).
% 
% EQUAÇÃO GOVERNANTE GERAL (Para cada estágio):
% dn/dt + d(g(T)*n)/da = -mu(T,a)*n
%
% ONDE:
%    n(t,a): Densidade de indivíduos (indivíduos por dia de idade)
%         t: Tempo Cronológico (dias)
%         a: Idade fisiológica (adimensional ou normalizada)
%      g(T): Taxa de advecção (velocidade de envelhecimento fisiológico 
%            dependente da temperatura T)
%   mu(T,a): Taxa de mortalidade instantânea (depende da temperatura T e 
%            da idade a, representando estresse térmico e senescência)
%
% MÉTODO DE RESOLUÇÃO NUMÉRICA:
% Método das Linhas (Method of Lines - MOL). A dimensão espacial da idade 
% 'a' é discretizada em compartimentos (caixas numéricas) usando diferenças 
% finitas (esquema Upwind). Isso transforma a EDP num grande sistema acoplado 
% de Equações Diferenciais Ordinárias (EDOs), que é resolvido ao longo do 
% tempo 't' utilizando o integrador 'ode15s' (adequado para sistemas rígidos).
% =========================================================================
% [MAT] Limpeza completa do ambiente para evitar conflitos de memória
clearvars; 
clc; 
close all;

% --- PADRONIZAÇÃO GRÁFICA GLOBAL ---
set(groot, 'defaultAxesFontSize', 12);                   % Eixos e números (Tick labels) tamanho 12
set(groot, 'defaultTextFontSize', 12);                   % Textos gerais e legendas tamanho 12
set(groot, 'defaultAxesTitleFontSizeMultiplier', 14/12); % Força os títulos para tamanho 14
% -----------------------------------

% =========================================================================
% 1. PARAMETRIZAÇÃO BIOLÓGICA E AMBIENTAL DO MODELO
% =========================================================================
% --- 1.1 Forçante Climática (Sazonalidade) ---
p.period = 365; % Período do ciclo térmico anual (dias)

% --- 1.2 Taxas de Desenvolvimento (Modelo Teórico de Brière-1) ---
% O modelo de Brière descreve o desenvolvimento assimétrico de insetos.
% a: constante empírica; Tmin: limiar térmico inferior; Tmax: limiar letal superior.
p.briere_stg1 = struct('gamma', 5.460229e-5, 'Tmin', 13.22, 'Tmax', 36.00); % Estágio 1 (Ovos/L1 a L3)
p.briere_stg2 = struct('gamma', 2.759856e-5, 'Tmin', 13.52, 'Tmax', 36.95); % Estágio 2 (L4 a L7/Pupas)

% --- 1.3 Parâmetros de Reprodução / Fecundidade (Curva Gaussiana) ---
  p.b_max = 93.63;  % Fecundidade máxima potencial (ovos/fêmea/dia) no pico térmico
p.T_opt_3 = 24.27;  % Temperatura ótima exata para a máxima oviposição (°C)
p.sigma_3 = 4.46;   % Desvio padrão da curva (largura da janela térmica reprodutiva)
     p.nu = 0.5843; % Razão sexual (Proporção de fêmeas na população de adultos)

% --- 1.4 Parâmetros de Mortalidade ---
% Limites térmicos absolutos para sobrevivência (fora deste intervalo, a mortalidade é total)
 p.mu_Tmax = 37.0;  % Limiar térmico superior para a função de mortalidade (°C)
 p.mu_Tmin = 12.0;  % Limiar térmico inferior para a função de mortalidade (°C)
p.mu_death = 15.0;  % Taxa letal padronizada (usada como teto para evitar infinitos numéricos)

% Imaturos: Mortalidade térmica modelada por um polinômio de 2º grau (Formato em 'U')
p.mu_stg1_coeffs = [2.095618e-4, -9.978089e-3, 1.212723e-1]; % Estágio 1 (Ovos/L1 a L3)
p.mu_stg2_coeffs = [1.763156e-4, -8.794485e-3, 1.125844e-1]; % Estágio 2 (L4 a L7/Pupas)

% Adultos: A mortalidade é dividida em dois componentes: Térmico e Senescência
    p.mu_stg3_coeffs = [3.885256e-3, -1.966856e-1, 2.648668];  % Coeficientes do estresse térmico
p.mu_stg3_senescence = 3.527091;                               % Fator de peso para a morte por velhice

% --- 1.5 Competição Intraespecífica (Densidade-Dependência) ---
% Controla o crescimento populacional explosivo (Modelo tipo Beverton-Holt)
p.q_3 = 2e-5; % Intensidade da competição entre adultos por locais de oviposição

% =========================================================================
% 2. DISCRETIZAÇÃO DA MALHA (GRID) PARA O MÉTODO DAS LINHAS (MOL)
% =========================================================================
% Define-se uma temperatura de referência padrão para normalizar o envelhecimento
T_ref = 25.0;
p.rate_stg1_ref = local_calc_briere(T_ref, p.briere_stg1);
p.rate_stg2_ref = local_calc_briere(T_ref, p.briere_stg2);

% Domínios de Idade Fisiológica Máxima de cada estágio na T_ref
p.A1 = 1 / p.rate_stg1_ref;  % Idade fisiológica onde o inseto muda de Estágio 1 para Estágio 2
p.A2 = 1 / p.rate_stg2_ref;  % Idade fisiologica onde o inseto muda de Estágio 2 para Estágio 3
p.A3 = 10.0;                 % Longevidade fisiológica máxima do inseto no Estágio 3

% 1. DEFINE A RESOLUÇÃO E CALCULA O PASSO PRIMEIRO
res_stg1 = 80;              % Número de divisões/caixas desejadas no Estágio 1
    p.da = p.A1 / res_stg1; % Passo de integração 'da' (Fica totalmente dinâmico)

% 2. CONSTRUÇÃO DA MALHA ESPACIAL DE IDADE
% Usamos linspace com "res_stg1 + 1" para garantir que gere os nós (bordas), 
% cravando matematicamente o 0 absoluto no início e o p.A1 absoluto no final.
p.ages_stg1 = linspace(0, p.A1, res_stg1 + 1)'; 

% Mantendo a sua regra original para os Estágios 2 e 3 (Pontos centrais)
p.ages_stg2 = (p.da/2 : p.da : p.A2)'; 
p.ages_stg3 = (p.da/2 : p.da : p.A3)'; 

% 3. EXTRAI O NÚMERO TOTAL FINAL DE COMPARTIMENTOS 
% Agora as variáveis p.N capturam o tamanho real e seguro para as matrizes do MOL
p.N1 = length(p.ages_stg1); 
p.N2 = length(p.ages_stg2);
p.N3 = length(p.ages_stg3);

% [MAT] OTIMIZAÇÃO DE PERFORMANCE: Pré-alocação dos índices vetoriais para o integrador
p.idx_stg1 = 1 : p.N1;
p.idx_stg2 = (p.N1 + 1) : (p.N1 + p.N2);
p.idx_stg3 = (p.N1 + p.N2 + 1) : (p.N1 + p.N2 + p.N3);

% Pré-calcula a curva estática de envelhecimento natural (senescência) dos adultos
p.mu_stg3_age_vec = local_calc_senescence(p.ages_stg3, p);

% Vetores deslocados linearmente (apenas para plotagem contínua no eixo X dos gráficos)
ages_plot_stg1 = p.ages_stg1;
ages_plot_stg2 = p.ages_stg2 + p.A1;
ages_plot_stg3 = p.ages_stg3 + p.A1 + p.A2;

% =========================================================================
% 3. CONDIÇÕES INICIAIS DO SISTEMA NUMÉRICO (t = 0)
% =========================================================================
% Utiliza-se uma distribuição gaussiana deslocada para representar a entrada gradual 
% da população inicial no sistema, garantindo estabilidade no esquema numérico.
n0_stg1 = gaussian_dist(p.ages_stg1, p.A1, 510, p.da);  % Condição inicial: 510 Ovos/L1 a L3 por hectare
n0_stg2 = gaussian_dist(p.ages_stg2, p.A2, 270, p.da);  % Condição inicial: 270 L4 a L7/Pupas por hectare
n0_stg3 = gaussian_dist(p.ages_stg3, p.A3, 220, p.da);  % Condição inicial: 220 Adultos por hectare

% Vetor de estado inicial global (Concatenação de todos os estágios)
y0 = [n0_stg1; n0_stg2; n0_stg3];

% Paleta de cores global para padronização das figuras
    colors = [0, 0.85, 0; 0, 0, 1; 1, 0, 0; 0, 0, 0; 0.6, 0.2, 0.8]; 
color_temp = [0.8500 0.3250 0.0980];
    c_vals = ['b'; 'r'];

% =========================================================================
% 4. ANÁLISE TEÓRICA DAS FUNÇÕES BIOLÓGICAS (Gráficos Preliminares 1 a 3)
% =========================================================================
disp('Gerando Gráficos das Funções Biológicas Globais (10°C a 40°C)...');
T_analysis = linspace(5, 45, 500); % Malha fina de temperatura para curvas contínuas

% 4.1 Calcula as taxas relativas de advecção (Brière)
g_stg1 = local_calc_briere(T_analysis, p.briere_stg1) / p.rate_stg1_ref;
g_stg2 = local_calc_briere(T_analysis, p.briere_stg2) / p.rate_stg2_ref;

% 4.2 Fecundidade (com e sem penalização por competição de Beverton-Holt)
     fec_pure = p.b_max * exp(-0.5 * ((T_analysis - p.T_opt_3)/p.sigma_3).^2);
N_competition = 1e+5; % Densidade arbitrária apenas para teste da curva
     fec_comp = fec_pure ./ (1 + p.q_3 * N_competition);

% 4.3 Funções de penalidade por mortalidade térmica (Separado por estágios)
        mu_stg1_analysis = local_calc_mu_imature(T_analysis, p.mu_stg1_coeffs);
        mu_stg2_analysis = local_calc_mu_imature(T_analysis, p.mu_stg2_coeffs);
mu_stg3_thermal_analysis = local_calc_mu_adult_thermal(T_analysis, p);

% 4.4 Superfície 2D de mortalidade do adulto (Combina Temperatura e Idade)
       ages_analysis = linspace(0, p.A3-0.1, 300);
mu_stg3_age_analysis = local_calc_senescence(ages_analysis, p);
            [TX, AX] = meshgrid(T_analysis, ages_analysis);
       mu_total_grid = local_calc_mu_adult_thermal(TX, p) .* local_calc_senescence(AX, p);

% --- PLOT FIG 1: Curvas de Brière (Envelhecimento) ---
figure(1); 
clf; 
set(gcf, 'Name', 'Taxa de Desenvolvimento (Brière-1)');
g_vals = [g_stg1; g_stg2];
hold on; 
for i = 1:2
    plot(T_analysis, g_vals(i,:), 'Color', c_vals(i,:), 'LineWidth', 2);
end
title('Taxa Desenvolvimento Imaturos (Brière-1)'); 
set_ax('Temperatura (°C)', 'Taxa Desenvolvimento Relativa g(T)', [10 40]);
xticks(10:2:40);
ylim([0 max(g_vals(:)) * 1.1]);
legend('g_1(T) (Estágio 1)', 'g_2(T) (Estágio 2)', 'Location', 'northeast');
grid on; 
box on;
hold off;

% --- PLOT FIG 2: Distribuição Gaussiana de Fecundidade ---
figure(2); 
clf; 
set(gcf, 'Name', 'Taxa de Fecundidade vs Temperatura');
f_vals = [fec_pure; fec_comp];
hold on; 
for i = 1:2
    plot(T_analysis, f_vals(i,:), 'Color', c_vals(i,:), 'LineWidth', 2);
end
title('Taxa Fecundidade Fêmeas Reprodutoras (Oviposição)'); 
set_ax('Temperatura (°C)', '\beta_{max} (Ovos/Fêmea/Dia)', [10 40]);
xticks(10:2:40);
legenda_1 = '\beta(T) sem Comp. (Densidade Baixa)';
legenda_2 = sprintf('B(T,N_3) com Comp. (N_3 = %d Indiv.)', N_competition);
legend({legenda_1, legenda_2}, 'Location', 'northeast');
grid on; 
box on;
hold off;

% --- PLOT FIG 3: Superfícies e Curvas de Decaimento (Mortalidade) ---
figure(3); 
clf; 
set(gcf, 'Name', 'Taxa de Mortalidade Térmica vs Senescência');

subplot(2,2,1);
plot(T_analysis, mu_stg1_analysis, 'Color', colors(1,:), 'LineWidth', 2); 
hold on;
plot(T_analysis, mu_stg2_analysis, 'Color', colors(2,:), 'LineWidth', 2); 
hold off;
title('Mortalidade Térmica Imaturos'); 
set_ax('Temperatura (°C)', 'Taxa Térmica \mu(T) (1/Dia)', [10 40]);
xticks(10:5:40);
legend('\mu_1(T)', '\mu_2(T)', 'Location', 'northwest', 'Orientation','vertical');

subplot(2,2,2);
plot(T_analysis, mu_stg3_thermal_analysis, 'Color', colors(3,:), 'LineWidth', 2); 
title('Mortalidade Térmica Adultos'); 
set_ax('Temperatura (°C)', 'Taxa Térmica \mu_3(T) (1/Dia)', [10 40]);
xticks(10:5:40);

subplot(2,2,3);
plot(ages_analysis, mu_stg3_age_analysis, 'Color', colors(4,:), 'LineWidth', 2); 
title('Mortalidade Senescência Adultos'); 
set_ax('Idade Fisiológica (Dias)', 'Taxa Senescência \mu_3(a) (1/Dia)', [0 p.A3]);
xticks(0:1:10);

subplot(2,2,4);
surf(TX, AX, mu_total_grid, 'EdgeColor', 'none'); 
title('Mortalidade Total Adultos'); 
set_ax('Temperatura (°C)', 'Idade Fisiológica (Dias)', [10 40]);
xticks(10:5:40); yticks(0:2:10); zticks(0:1000:4000);
zlabel('Taxa Total \mu_3(T,a) (1/Dia)'); 
ylim([0 p.A3]); zlim([0 4000]); clim([0 4000]);
shading interp; 
view(3); 
colormap('jet'); 
colorbar; 

% =========================================================================
% 5. SIMULAÇÃO EM LOTE (SWEEP): VARREDURA CLIMÁTICA E SAZONAL
% =========================================================================
% Esta seção resolve o modelo para múltiplas combinações de temperatura média 
% e amplitude sazonal, extraindo o potencial de infestação para cada cenário.
disp('Iniciando simulações de todos os cenários climáticos (ODE15s)...');

   T_averages_test = [12, 16, 20, 24, 28, 32, 36, 40]; % Malha de Temperaturas Médias Anuais
   amplitudes_test = [0, 3, 6, 9];         % Malha de Amplitudes Sazonais (Flutuação)
   time_simulation = [0 1275];             % Janela de integração (t0 a tf)[Cenários: 180d(6m), 270d(9m), 1185d(3a3m), 1275d(3a6m), 1640d(4a6m)]
day_start_analysis = 180;                  % Início do período para extrair picos (ignora transiente)[Cenários: 60d, 180d]
  day_end_analysis = time_simulation(2);   % Fim do período para extrair picos [Cenários: 180d, 270d, 1185d, 1275d, 1640d]

% Criação do espaço de combinações (Grid paramétrica)
[T_grid, Amp_grid] = ndgrid(T_averages_test, amplitudes_test);
             T_vec = T_grid(:); 
           Amp_vec = Amp_grid(:);
       N_scenarios = length(T_vec);

% Células e matrizes para armazenamento dos resultados
       t_all = cell(N_scenarios, 1); 
       y_all = cell(N_scenarios, 1);
Pop_Min_stg1 = zeros(N_scenarios, 1); Pop_Max_stg1 = zeros(N_scenarios, 1);
Pop_Min_stg2 = zeros(N_scenarios, 1); Pop_Max_stg2 = zeros(N_scenarios, 1);
Pop_Min_stg3 = zeros(N_scenarios, 1); Pop_Max_stg3 = zeros(N_scenarios, 1);

% Configurações do Integrador Numérico (Garante que populações não fiquem negativas)
opts_sweep = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NonNegative', 1:length(y0));

tic; 
% Laço paralelo (parfor) acelera enormemente as simulações da varredura
parfor k = 1:N_scenarios 
    % 1. Pega os valores específicos dessa iteração
    T_current = T_vec(k);
    Amp_current = Amp_vec(k);
    
    % 2. Chama a função encapsulada
    % Ela retorna exatamente os dados que você precisa salvar nos vetores
    [t_all{k}, y_all{k}, Pop_Min_stg1(k), Pop_Max_stg1(k), Pop_Min_stg2(k), Pop_Max_stg2(k), ...
    Pop_Min_stg3(k), Pop_Max_stg3(k)] = simulate_scenario(T_current, Amp_current, p, time_simulation, y0, ...
    opts_sweep, day_start_analysis, day_end_analysis);
end

% Reconstrução das matrizes 2D (Média x Amplitude) para os gráficos de varredura
Mat_Min_stg1 = reshape(Pop_Min_stg1, length(T_averages_test), length(amplitudes_test));
Mat_Max_stg1 = reshape(Pop_Max_stg1, length(T_averages_test), length(amplitudes_test));
Mat_Min_stg2 = reshape(Pop_Min_stg2, length(T_averages_test), length(amplitudes_test));
Mat_Max_stg2 = reshape(Pop_Max_stg2, length(T_averages_test), length(amplitudes_test));
Mat_Min_stg3 = reshape(Pop_Min_stg3, length(T_averages_test), length(amplitudes_test));
Mat_Max_stg3 = reshape(Pop_Max_stg3, length(T_averages_test), length(amplitudes_test));

fprintf('Simulações Concluídas em %.2f segundos.\n', toc);

% =========================================================================
% 6. GRÁFICOS DE VARREDURA (BIFURCAÇÃO ECOLÓGICA / BOXPLOTS)
% =========================================================================
  titles_var = {'Estágio 1', 'Estágio 2', 'Estágio 3'};
   Datas_Min = {Mat_Min_stg1, Mat_Min_stg2, Mat_Min_stg3}; 
   Datas_Max = {Mat_Max_stg1, Mat_Max_stg2, Mat_Max_stg3};
str_interval = sprintf('Período de Monitoramento: Dia %g a %g', day_start_analysis, day_end_analysis);
 offset_base = 0.2; 
           w = 0.08; 
 positions_y = 1:length(T_averages_test);

% Gráficos estilo "Tornado" ou boxplot horizontal demonstrando os mín e máx populacionais
for k = 1:3
    figure(3+k);
    clf; 
    set(gcf, 'Name', ['Amplitude Populacional vs Temperatura: ' titles_var{k}]); 
    hold on;
          d_min = round(Datas_Min{k}); 
          d_min(d_min < 0) = 0;
          d_max = round(Datas_Max{k}); 
          d_max(d_max < 0) = 0;
    max_x_geral = max(10, max(max(d_max)));
          h_leg = zeros(1, length(amplitudes_test)); 
    legend_strs = cell(1, length(amplitudes_test));
    
    for j = 1:length(amplitudes_test)
        y_plot = positions_y + ((j - (length(amplitudes_test)+1)/2) * offset_base);
        for i = 1:length(T_averages_test)
            x_min = d_min(i, j); 
            x_max = d_max(i, j); 
            y_c = y_plot(i);
            patch([x_min, x_min, x_max, x_max], [y_c - w, y_c + w, y_c + w, y_c - w], colors(j,:), 'FaceAlpha', 0.5, 'EdgeColor', colors(j,:), 'LineWidth', 2.5);
            text_val = sprintf('%d - %d', x_min, x_max); 
            if x_min == x_max 
                text_val = sprintf('%d', x_max); 
            end
            text(x_max + (max_x_geral * 0.015), y_c, text_val, 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 11, 'Color', 'k');
        end
              h_leg(j) = patch([NaN NaN NaN NaN], [NaN NaN NaN NaN], colors(j,:), 'FaceAlpha', 0.5, 'EdgeColor', colors(j,:), 'LineWidth', 2.5);
        legend_strs{j} = sprintf('Amp %d °C', amplitudes_test(j));
    end
    title(['Dinâmica Populacional vs Amplitude Térmica: ' titles_var{k}], str_interval, 'FontWeight', 'bold');
    set_ax('Nº de Indivíduos (Mín. a Máx.)', 'Temperatura Média Anual (°C)', [0, max_x_geral * 1.25]); 
    if k == 1 || k == 2
       xticklabels([0:1.0e+07:max_x_geral * 1.25, max_x_geral]);
    elseif k == 3
       xticklabels([0:0.25e+07:max_x_geral * 1.25, max_x_geral]);
    end
    ylim([0.5, length(T_averages_test) + 0.5]); 
    yticks(positions_y);
    yticklabels(string(T_averages_test));
    for y_sep = 1.5 : 1 : (length(T_averages_test) - 0.5)
        yline(y_sep, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1); 
    end
    grid on; 
    box on; 
    legend(h_leg, legend_strs, 'Location', 'northeast'); 
    hold off;
end

% =========================================================================
% 7. PAINEL DE CONTROLE: GERAÇÃO DE GRÁFICOS DO CENÁRIO ESPECÍFICO
% =========================================================================
% A partir deste ponto, o código foca na dinâmica temporal de um ÚNICO cenário.
% => DEFINIÇÃO DO CENÁRIO ESCOLHIDO PARA PLOTAR (DEVE EXISTIR NA VARREDURA):
  T_chosen = 24;   
Amp_chosen = 0;

disp(['Extraindo e plotando dados completos para T = ' num2str(T_chosen) '°C e Amp = ' num2str(Amp_chosen) '°C...']);
% Procura o cenário na memória
idx_plot = find(T_vec == T_chosen & Amp_vec == Amp_chosen);
if isempty(idx_plot)
   error('O cenário escolhido não foi simulado na varredura.'); 
end

% Extração das variáveis temporais do cenário escolhido
            t = t_all{idx_plot};
        y_sol = y_all{idx_plot};
  p.T_media_C = T_chosen; 
p.amplitude_C = Amp_chosen;
     p.T_func = @(t) p.T_media_C + p.amplitude_C * cos(2 * pi * t / p.period);

% Separação espacial das densidades e cálculo dos totais por estágio
n_sol_stg1 = y_sol(:, p.idx_stg1);   total_stg1 = sum(n_sol_stg1, 2)*p.da;
n_sol_stg2 = y_sol(:, p.idx_stg2);   total_stg2 = sum(n_sol_stg2, 2)*p.da;
n_sol_stg3 = y_sol(:, p.idx_stg3);   total_stg3 = sum(n_sol_stg3, 2)*p.da;
  
% --- RECONSTRUÇÃO DOS FLUXOS INTERNOS (Para fins de gráficos) ---
            temp_t_C = p.T_func(t); 
   v_rate_g_stg2_eff = local_calc_briere(temp_t_C, p.briere_stg2) / p.rate_stg2_ref;
        mu_stg1_base = local_calc_mu_imature(temp_t_C, p.mu_stg1_coeffs);
        mu_stg2_base = local_calc_mu_imature(temp_t_C, p.mu_stg2_coeffs);
mu_stg3_thermal_plot = max(1e-4, local_calc_mu_adult_thermal(temp_t_C, p));

% Reconstrução da dinâmica de reprodução
v_rate_b_potencial = p.b_max * exp(-0.5 * ((temp_t_C - p.T_opt_3)/p.sigma_3).^2);
     v_rate_b_real = v_rate_b_potencial ./ (1 + p.q_3 .* total_stg3);

% Reconstrução dos fluxos de transição entre fases
v_establishment_flux = (local_calc_briere(temp_t_C, p.briere_stg1)/p.rate_stg1_ref) .* n_sol_stg1(:, end);
    v_emergence_flux = v_rate_g_stg2_eff .* n_sol_stg2(:, end);

% Reconstrução do status demográfico do adulto
        n_sol_stg3_safe = max(0, n_sol_stg3);
        total_stg3_safe = max(1, sum(n_sol_stg3_safe, 2));
    mu_senescence_media = (n_sol_stg3_safe * p.mu_stg3_age_vec) ./ total_stg3_safe;
      mu_stg3_vec_total = min(p.mu_death, mu_stg3_thermal_plot .* mu_senescence_media);

% --- FLAG DE CONTROLE GRÁFICO (Sazonal vs Constante) ---
is_temp_const = max(temp_t_C) == min(temp_t_C);
if is_temp_const
    str_temp = sprintf(' (Temp. de %.1f °C)', temp_t_C(1)); 
else
    str_temp = ''; 
end

% =========================================================================
% --- FIGURA 7: POTENCIAL TÉRMICO E ABUNDÂNCIA MÁXIMA GLOBAL ---
% =========================================================================
% Une as curvas teóricas (Brière e Fecundidade) com os picos reais 
% simulados gerados pela matriz de sweep. Mostra a "Janela de Viabilidade".
figure(7); 
clf; 
set(gcf, 'Name', 'População Máxima e Janela Térmica');
true_max_pops = [max(Mat_Max_stg1(:)), max(Mat_Max_stg2(:)), max(Mat_Max_stg3(:))];
    curves_f7 = {local_calc_briere(T_analysis, p.briere_stg1), local_calc_briere(T_analysis, p.briere_stg2)};
     T_params = {p.briere_stg1, p.briere_stg2};
titles_var_f7 = {'','Curvas Populacionais: Amplitudes Máximas vs Janelas Térmicas',''};
    titles_f7 = {'Estágio 1', 'Estágio 2'};

for i = 1:2
    subplot(1,3,i);
    curve = (curves_f7{i} / max(curves_f7{i})) * true_max_pops(i);
    [~, idx] = max(curve); 
    plot(T_analysis, curve, 'Color', colors(i,:), 'LineWidth', 2); 
    title(titles_var_f7{i}, titles_f7{i}, 'FontWeight', 'bold');
    hold on;
    xline(T_params{i}.Tmin, '--', 'Color', colors(4,:), 'Label', sprintf('T_{min} = %.1f °C', T_params{i}.Tmin), 'FontWeight', 'bold', 'FontSize', 11);    
    xline(T_analysis(idx), '--', 'Color', colors(4,:), 'Label', sprintf('T_{opt} = %.1f °C', T_analysis(idx)), ...
        'FontWeight', 'bold', 'FontSize', 11, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'center');  
    xline(T_params{i}.Tmax, '--', 'Color', colors(4,:), 'Label', sprintf('T_{max} = %.1f °C', T_params{i}.Tmax), 'FontWeight', 'bold', 'FontSize', 11);   
    xticks(10:5:40);
    if i == 1
        set_ax('Temperatura (°C)', 'Nº Máximo de Indivíduos', [10 40]);
    else
        set_ax('Temperatura (°C)', '', [10 40]); 
    end
    if true_max_pops(i) > 0
       ylim([0 true_max_pops(i)*1.15]); 
    end
    hold off;
end

subplot(1,3,3);
y_at_min = p.b_max * exp(-0.5 * ((p.mu_Tmin - p.T_opt_3)/p.sigma_3).^2);
y_at_max = p.b_max * exp(-0.5 * ((p.mu_Tmax - p.T_opt_3)/p.sigma_3).^2);
curve_3_base = fec_pure - (y_at_min + ((y_at_max - y_at_min)/(p.mu_Tmax - p.mu_Tmin)) * (T_analysis - p.mu_Tmin));
curve_3_base(T_analysis <= p.mu_Tmin | T_analysis >= p.mu_Tmax | curve_3_base < 0) = 0;
if max(curve_3_base) > 0
    curve_3 = (curve_3_base / max(curve_3_base)) * true_max_pops(3);
else 
    curve_3 = curve_3_base; 
end
[~, idx3] = max(curve_3);
plot(T_analysis, curve_3, 'Color', colors(3,:), 'LineWidth', 2); 
title('Estágio 3'); 
hold on;
xline(p.mu_Tmin, '--', 'Color', colors(4,:), 'Label', sprintf('T_{min} = %.1f °C', p.mu_Tmin), 'FontWeight', 'bold', 'FontSize', 11);  
xline(T_analysis(idx3), '--', 'Color', colors(4,:), 'Label', sprintf('T_{opt} = %.1f °C', T_analysis(idx3)), ...
    'FontWeight', 'bold', 'FontSize', 11, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'center');  
xline(p.mu_Tmax, '--', 'Color', colors(4,:), 'Label', sprintf('T_{max} = %.1f °C', p.mu_Tmax), 'FontWeight', 'bold', 'FontSize', 11);
xticks(10:5:40); 
set_ax('Temperatura (°C)', '', [10 40]);
if true_max_pops(3) > 0 
   ylim([0 true_max_pops(3)*1.15]); 
end
hold off;

% =========================================================================
% --- FIGS 8 e 9: CONDIÇÃO INICIAL E EVOLUÇÃO TEMPORAL (Densidades) ---
% =========================================================================
  n0_stgs_all = {n0_stg1, n0_stg2, n0_stg3}; 
     ages_all = {ages_plot_stg1, ages_plot_stg2, ages_plot_stg3}; 
       totals = {total_stg1, total_stg2, total_stg3}; 
    titles_f8 = {'Estágio 1', 'Estágio 2', 'Estágio 3'};
titles_var_f8 = {'','Distribuições Populacionais Iniciais da Broca-da-Cana',''};
 
% --- FIGURA 8 (Distribuições t=0) ---
figure(8); 
clf; 
set(gcf, 'Name', 'Condições Iniciais da Praga');
for i = 1:3
    subplot(1,3,i);
    x_stage = ages_all{i}(:); 
    y_stage = n0_stgs_all{i}(:);
    plot(x_stage, y_stage, 'Color', colors(i,:), 'LineWidth', 2);   
    if isempty(titles_var_f8{i})
       title(titles_f8{i}, 'FontWeight', 'bold');
    else
       title({titles_var_f8{i}, titles_f8{i}}, 'FontWeight', 'bold');
    end
    set_ax('Idade Fisiológica (Dias)', sprintf('Densidade n_%d(0,a) (Indiv./Ha.)', i), [x_stage(1), x_stage(end)]);
    if i == 1
       xticks(linspace(x_stage(1),x_stage(end),5));
       xticklabels({'0', '5', '9', '14', '19'});
    elseif i == 2
       xticks(linspace(x_stage(1),x_stage(end),5));
       xticklabels({'19', '28', '37', '46', '55'});
    elseif i == 3
       xticks(linspace(x_stage(1),x_stage(end),5));
       xticklabels({'55', '58', '60', '63', '65'});
    end
    y_max = max(y_stage);
    if y_max > 0 
       ylim([0, y_max]); 
    end
    grid on;
end
    
% --- FIGURA 9 (Curvas de crescimento) ---
figure(9); 
clf; 
set(gcf, 'Name', 'Evolução Populacional Temporal');
for i = 1:3
    if i >= 3
       x_lbls_f9 = 'Tempo Cronológico (Dias)';
    else
       x_lbls_f9 = ' ';
    end
    subplot(3,1,i); 
    titles_f9 = {['Evolução Populacional vs Temperatura Constante: ' titles_f8{i}], ['Evolução Populacional vs Sazonalidade Térmica: ' titles_f8{i}]};
    if is_temp_const
       plot(t, totals{i}, 'Color', colors(i,:), 'LineWidth', 2);
       title([titles_f9{1} ' ' str_temp]); 
       set_ax(x_lbls_f9, 'Nº Total de Indivíduos', [min(t), max(t)]); 
       xticks(unique([min(t):100:max(t), max(t)]));
       xtickangle(45);
       grid on;
    else
        yyaxis left;
        plot(t, totals{i}, 'Color', colors(i,:), 'LineWidth', 2);
        title(titles_f9{2});
        set_ax(x_lbls_f9, 'Nº Total de Indivíduos', [min(t), max(t)]);
        xticks(unique([min(t):100:max(t), max(t)]));
        xtickangle(45);
        grid on;
        yyaxis right;
        plot(t, temp_t_C, '--', 'Color', color_temp, 'LineWidth', 0.8);
        ylabel('Temperatura (°C)'); 
        set(gca, 'YColor', color_temp);
    end 
end

% --- FIGS 10 a 12: SUPERFÍCIES 3D (ESPAÇO TEMPO-IDADE) ---
% Mostra graficamente a advecção (movimento da população pela matriz etária)
idx_plot_3d = 1:max(1, floor(length(t)/150)):length(t);
  t_plot_3d = t(idx_plot_3d);
       sols = {n_sol_stg1, n_sol_stg2, n_sol_stg3}; 
     figs3D = {'Estágio 1', 'Estágio 2', 'Estágio 3'};

for i = 1:3
    figure(9+i); 
    clf; 
    set(gcf, 'Name', ['Superfície Populacional-Etária 3D: ' figs3D{i}]); 
    surf(ages_all{i}, t_plot_3d, sols{i}(idx_plot_3d, :), 'EdgeColor', 'none');
    title(['Superfície Populacional-Etária: ' figs3D{i}], 'FontWeight', 'bold');
    set_ax('Idade Fisiológica (Dias)', 'Tempo Cronológico (Dias)'); 
    zlabel('Nº Total de Indivíduos');
    shading interp; 
    colormap('jet'); 
    colorbar; 
    view(3);  
end

% --- FIG 13: FORÇANTE CLIMÁTICA (Termômetro do Sistema) ---
figure(13); 
clf; 
set(gcf, 'Name', 'Monitoramento Climático (Sazonalidade Térmica)');
plot(t, temp_t_C, 'Color', color_temp, 'LineWidth', 2.0); 
title('Sazonalidade Térmica Cossenoidal T(t)'); 
set_ax('Tempo Cronológico (Dias)', 'Temperatura (°C)', [min(t) max(t)]);  
xticks(unique([min(t):100:max(t), max(t)]));
xtickangle(45);
yline(p.T_media_C, '--k');

% --- FIGS 14 e 15: FLUXOS DEMOGRÁFICOS DE CONTORNO (Limites do Domínio) ---
figure(14); 
clf; 
set(gcf, 'Name', 'Fluxos Demográficos da Praga');
x_lbls_f14 = {'', 'Tempo Cronológico (Dias)'};
y_vals_f14 = {v_establishment_flux(:), v_emergence_flux(:)}; 
y_lbls_f14 = {'Larvas/Dia'; 'Adultos/Dia'};
titles_f14 = {'Fluxo de Estabelecimento: Entrada no Colmo (Estágio 1 --> Estágio 2)'; 'Fluxo de Emergência: Maturação da Praga (Estágio 2 --> Estágio 3)'};
for i = 1:2
    subplot(2,1,i);
    if is_temp_const
       plot(t, y_vals_f14{i}, 'Color', c_vals(i,:), 'LineStyle', '-', 'LineWidth', 2);
       title([titles_f14{i} str_temp]); 
       set_ax(x_lbls_f14{i}, y_lbls_f14{i}, [min(t) max(t)]);
       xticks(unique([min(t):100:max(t), max(t)])); 
       xtickangle(45);
       grid on;
    else 
       yyaxis left; 
       plot(t, y_vals_f14{i}, 'Color', c_vals(i,:), 'LineStyle', '-', 'LineWidth', 2);
       title(titles_f14{i});
       set_ax(x_lbls_f14{i}, y_lbls_f14{i}, [min(t) max(t)]);
       xticks(unique([min(t):100:max(t), max(t)]));
       xtickangle(45); 
       grid on;
       yyaxis right; 
       plot(t, temp_t_C(:), '--', 'Color', color_temp, 'LineWidth', 0.8); 
       ylabel('Temperatura (°C)');
    end
end

figure(15); 
clf; 
y_vals_f15 = {v_rate_b_potencial(:), v_rate_b_real(:)};
set(gcf, 'Name', 'Taxa de Fecundidade e Competição Temporal');
if is_temp_const
    hold on;
    plot(t, y_vals_f15{1}, 'Color', c_vals(1,:), 'LineStyle', '-', 'LineWidth', 2);  
    plot(t, y_vals_f15{2}, 'Color', c_vals(2,:), 'LineStyle', '-', 'LineWidth', 2); 
    hold off;
    title(['Oviposição Máxima Diária (Fecundidade)' str_temp]); 
    legend('\beta(T) sem Comp.', 'B(T,N_3) com Comp.', 'Location', 'north','Orientation','horizontal');
    set_ax('Tempo Cronológico (Dias)', '\beta_{max} (Ovos/Fêmea/Dia)', [min(t) max(t)]);
    xticks(unique([min(t):100:max(t), max(t)]));
    xtickangle(45);
    ylim([0 max(y_vals_f15{1}) * 1.1]);
    grid on;
    box on;
else 
     yyaxis left; 
     hold on; 
     plot(t, y_vals_f15{1}, 'Color', c_vals(1,:), 'LineStyle', '-', 'LineWidth', 2);
     plot(t, y_vals_f15{2}, 'Color', c_vals(2,:), 'LineStyle', '-', 'LineWidth', 2); 
     hold off; 
     title('Oviposição Máxima Diária (Fecundidade)'); 
     set_ax('Tempo Cronológico (Dias)', '\beta_{max} (Ovos/Fêmea/Dia)', [min(t) max(t)]);
     xticks(unique([min(t):100:max(t), max(t)]));
     ylim([0 max(y_vals_f15{1}) * 1.1]);
     xtickangle(45);
     yyaxis right; 
     plot(t, temp_t_C, '--', 'Color', color_temp, 'LineWidth', 0.8);
     legend('\beta(T) sem Comp.', 'B(T,N_3) com Comp.', 'Temp. (˚C)', 'Location', 'north','Orientation','horizontal');
     ylabel('Temperatura (°C)');
     ylim([min(temp_t_C) - 1 max(temp_t_C) * 1.1]);
     grid on;
     box on;
end

% --- FIG 16: DINÂMICA DAS TAXAS DE MORTALIDADE TEMPORAIS ---
figure(16); 
clf; 
set(gcf, 'Name', 'Análise de Mortalidade Temporal');
y_vals_f16 = {[mu_stg1_base(:), mu_stg2_base(:)], mu_stg3_thermal_plot(:), mu_senescence_media(:), mu_stg3_vec_total(:)};
titles_f16 = {'Mortalidade Térmica Imaturos', 'Mortalidade Térmica Adultos', 'Mortalidade Senescência Adultos', 'Mortalidade Total Adultos'};
y_lbls_f16 = {'Taxa Térmica \mu_i(T) (1/Dia)', 'Taxa Térmica \mu_3(T) (1/Dia)', 'Taxa Senescência \mu_3(a) (1/Dia)', 'Taxa Total \mu_3(T,a) (1/Dia)'};
t = t(:); temp_col = temp_t_C(:);

for i = 1:4
    subplot(2,2,i);
    if i >= 3
       x_lbls_f16 = 'Tempo Cronológico (Dias)';
    else
       x_lbls_f16 = '';
    end
    max_y = max(y_vals_f16{i}(:));
    limite_y = max_y * 1.25;
    if limite_y == 0
       limite_y = 1;  
    end
    
    if is_temp_const 
        if i == 1
            hold on;
            plot(t, y_vals_f16{i}(:,1), 'Color', colors(1,:), 'LineStyle', '-', 'LineWidth', 2);
            plot(t, y_vals_f16{i}(:,2), 'Color', colors(2,:), 'LineStyle', '-', 'LineWidth', 2);
            hold off;
            title([titles_f16{i} str_temp]);
            legend('\mu_1(T)', '\mu_2(T)', 'Location', 'north', 'Orientation','horizontal'); 
            set_ax(x_lbls_f16, y_lbls_f16{i}, [min(t) max(t)]);
            xticks(unique([min(t):100:max(t), max(t)])); 
            xtickangle(45);
            ylim([0 limite_y]);    
            box on; 
        else
            plot(t, y_vals_f16{i}, 'Color', colors(i+1,:), 'LineWidth', 2); 
            title([titles_f16{i} str_temp]);
            set_ax(x_lbls_f16, y_lbls_f16{i}, [min(t) max(t)]);
            xticks(unique([min(t):100:max(t), max(t)]));
            xtickangle(45);
            ylim([0 limite_y]);
            box on; 
        end  
    else
        % --- CASO TEMPERATURA VARIÁVEL (YYAXIS) ---
        yyaxis left; 
        if i == 1
            hold on; % Ativa para plotar as duas curvas de mortalidade
            plot(t, y_vals_f16{i}(:,1), 'Color', colors(1,:), 'LineStyle', '-', 'LineWidth', 2);
            plot(t, y_vals_f16{i}(:,2), 'Color', colors(2,:), 'LineStyle', '-', 'LineWidth', 2);
            title(titles_f16{i});
            set_ax(x_lbls_f16, y_lbls_f16{i}, [min(t) max(t)]);
            xticks(unique([min(t):100:max(t), max(t)]));
            xtickangle(45);
            ylim([0 limite_y]);
        else
            plot(t, y_vals_f16{i}, 'Color', colors(i+1,:), 'LineWidth', 2);
            title(titles_f16{i});
            set_ax(x_lbls_f16, y_lbls_f16{i}, [min(t) max(t)]);
            xticks(unique([min(t):100:max(t), max(t)]));
            xtickangle(45);
            ylim([0 limite_y]);
        end      
        yyaxis right; 
        plot(t, temp_col, '--', 'Color', color_temp, 'LineWidth', 0.8); 
        ylabel('Temperatura (°C)'); 
        ylim([min(temp_col) - 1 max(temp_col) * 1.15]);
        set(gca, 'YColor', color_temp);
        if i == 1
           legend('\mu_1(T)', '\mu_2(T)', 'Temp. (˚C)', 'Location', 'north','Orientation','horizontal');
        end
        
        grid on; 
        box on;          % Ativa o contorno para o lado direito (ativo no momento)
        yyaxis left;     % Alterna brevemente para o esquerdo...
        box on;          % ...e força o contorno fechar perfeitamente no topo e laterais
    end
end

% --- FIG 17: PERFIL ETÁRIO GERAL (Agrupamento Histórico) ---
              N_total_stages = trapz(t, y_sol, 1)'; 
                      x_cont = [ages_plot_stg1; ages_plot_stg2; ages_plot_stg3];
            window_smoothing = round(length(N_total_stages) * 0.10);
                    y_smooth = smoothdata(N_total_stages, 'gaussian', window_smoothing);
                      n_tail = round(length(ages_plot_stg3) * 0.40);
y_smooth(end-n_tail+1 : end) = y_smooth(end-n_tail+1 : end) .* ((cos(linspace(0, pi, n_tail)') + 1) / 2);
                        idxs = {1:p.N1, p.N1:(p.N1 + p.N2), (p.N1 + p.N2):length(x_cont)};

figure(17); 
clf; 
set(gcf, 'Name', 'Perfil Populacional-Etário Global'); 
hold on;
h = zeros(1,3);
for i = 1:3
    x_poly = [x_cont(idxs{i})', fliplr(x_cont(idxs{i})')]; 
    y_poly = [y_smooth(idxs{i})', zeros(1, length(idxs{i}))];
    fill(x_poly, y_poly, colors(i,:), 'EdgeColor', colors(i,:), 'LineWidth', 2.5, 'FaceAlpha', 0.25);
    plot(x_cont(idxs{i}), y_smooth(idxs{i}), 'Color', colors(i,:), 'LineWidth', 3);
    h(i) = fill(NaN, NaN, colors(i,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(i,:), 'LineWidth', 2.5);
end
box on; 
set(gca, 'TickDir', 'in');
title('Estrutura Populacional-Etária Global Acumulada', sprintf('Duração de Monitoramento: %g Dias', t(end)), 'FontWeight', 'bold');
legend(h, {'Estágio 1', 'Estágio 2', 'Estágio 3'}, 'Location', 'northeast');
set_ax('Idade Fisiológica (Dias)', 'Nº Total de Indivíduos', [0, p.A1 + p.A2 + p.A3]);
xticks([0:5:max(x_cont), 66])
ylim([0, max(y_smooth) * 1.1]);
hold off;

% =========================================================================
% SALVAMENTO AUTOMÁTICO DE TODOS OS GRÁFICOS
% =========================================================================
fprintf('\n=======================================================\n');
disp('SALVANDO GRÁFICOS...');

% 1. Cria uma pasta no seu computador para não bagunçar seus arquivos
destination_folder = 'Resultados_Imagens';
if ~exist(destination_folder, 'dir')
    mkdir(destination_folder);
end

% 2. Encontra todas as figuras que estão abertas
figures = findall(groot, 'Type', 'figure');
[~, order] = sort([figures.Number]); % Pega os números das figuras e descobre a ordem crescente
figures = figures(order);            % Reorganiza o vetor do menor para o maior

% 3. Faz um loop salvando uma a uma
for i = 1:length(figures)
    fig = figures(i);
    
    % Pega o número da figura (ex: 16)
    num_fig = fig.Number; 
    
    % Cria o nome do arquivo (ex: Resultados_Imagens/Figura_16.png)
    file_name = fullfile(destination_folder, sprintf('Figura_%d.png', num_fig));
    
    % Encontra todos os eixos (gráficos) dentro desta figura específica
    figure_axes = findall(fig, 'Type', 'axes');
    for j = 1:length(figure_axes)
        % Trava os limites dos eixos para que o exportgraphics não os recalcule
        figure_axes(j).XLimMode = 'manual';
        figure_axes(j).YLimMode = 'manual';
        
        % Opcional: Garante que a caixa de fora cubra toda a linha
        figure_axes(j).Clipping = 'on';
    end
    
    drawnow;     
    pause(0.5);  % Dá meio segundo para a memória estabilizar (muito útil para 300 DPI)
    
    % Salva a figura em formato PNG com alta resolução (300 DPI)
    exportgraphics(fig, file_name, 'Resolution', 250);
    
    fprintf('Salvo com sucesso: %s\n', file_name);
end

disp('Todos os gráficos foram salvos na pasta "Resultados_Imagens"!');
fprintf('=======================================================\n');

% =========================================================================
% RELATÓRIOS DE CONSOLE E EXPORTAÇÃO DO CENÁRIO
% =========================================================================
% Imprime no terminal do MATLAB o Paradoxo de Malha e a real T_opt do modelo

fprintf('\n=======================================================\n');
fprintf('   RESULTADOS DA VARREDURA TÉRMICA E POTENCIAL TEÓRICO \n');
fprintf('=======================================================\n');

       datas_sweep = {Mat_Max_stg1, Mat_Max_stg2, Mat_Max_stg3};
      stages_names = {'Ovos/Larvas L1 a L3', 'Larvas L4 a L7/Pupas', 'Adultos'};
                c1 = local_calc_briere(T_analysis, p.briere_stg1);
                c2 = local_calc_briere(T_analysis, p.briere_stg2);
                c3 = curve_3_base; % Reaproveita a curva limpa sem bicos gerada na Fig 7
theoretical_curves = {c1, c2, c3};

for i = 1:3
    % Encontra qual coordenada (linha/coluna) da matriz de resultados contém o valor máximo
    [max_pop, linear_idx] = max(datas_sweep{i}(:));
               [row, col] = ind2sub(size(datas_sweep{i}), linear_idx);
              T_simulated = T_averages_test(row);
            Amp_simulated = amplitudes_test(col);
    
    % Encontra a T_opt teórica analítica (matriz contínua fina de teste)
       [~, idx_opt] = max(theoretical_curves{i});
    T_opt_theorical = T_analysis(idx_opt);
    
    fprintf('ESTÁGIO %d: %s\n', i, stages_names{i});
    fprintf('  -> Temperatura Ótima Verdadeira (T_opt): %.1f °C\n', T_opt_theorical);
    fprintf('  -> Maior População Encontrada: %.0f indivíduos\n', max_pop);
    fprintf('     (Gerada no cenário simulado mais próximo: Temp = %d °C | Amp = %d °C)\n', T_simulated, Amp_simulated);
    fprintf('-------------------------------------------------------\n');
end

fprintf('\nTodos os Gráficos (1 a 17) gerados e atualizados com sucesso!\n');

% =========================================================================
% 8. FUNÇÕES LOCAIS (KERNEL MATEMÁTICO DA EDP)
% =========================================================================

% --- SISTEMA DE EDOs RESULTANTE DA DISCRETIZAÇÃO (MÉTODO DAS LINHAS) ---
function dydt = odefun_3stage(t, y, p)
    % 1. Atualiza a forçante externa (Temperatura no tempo t atual do solver)
    T = p.T_func(t);
    
    % 2. Mapeia e divide o vetor de estado global (y) nas densidades locais
    n1 = y(p.idx_stg1); % Ovos/L1 a L3
    n2 = y(p.idx_stg2); % L4 a L7/Pupas
    n3 = y(p.idx_stg3); % Adultos
    
    % Cálculo populacional macro (Integral da densidade em relação à idade 'a')
    N3_total = sum(n3)*p.da;
    
    % 3. Taxas Direcionais de Advecção (Velocidade na malha)
    % Velocidade normalizada na temperatura padrão (Envelhecimento térmico)
    g_stg1 = local_calc_briere(T, p.briere_stg1) / p.rate_stg1_ref;
    g_stg2 = local_calc_briere(T, p.briere_stg2) / p.rate_stg2_ref;
    g_stg3 = 1.0; % Adultos experimentam 1 dia fisiológico a cada 1 dia cronológico
    
    % 4. Módulo de Mortalidade / Perdas (Decaimento populacional contínuo)
    mu_stg1 = local_calc_mu_imature(T, p.mu_stg1_coeffs);
    mu_stg2 = local_calc_mu_imature(T, p.mu_stg2_coeffs);
    
    % A mortalidade do adulto flutua individualmente pelo vetor de idades (senescência)
    mu_stg3_total_vec = local_calc_mu_adult_thermal(T, p) * p.mu_stg3_age_vec;
    
    % 5. Condição de Contorno Superior (Recrutamento no nó 1 do estágio 1)
    fec_real = (p.b_max * exp(-0.5 * ((T - p.T_opt_3)/p.sigma_3)^2)) / (1 + p.q_3 * N3_total);
    
    % O fluxo de novos ovos é gerado APENAS pela fração de fêmeas na população (parâmetro nu)
    Flux_in_stg1 = fec_real * (N3_total * p.nu); 
    
    % Inicializa o vetor espelho de taxas de variação temporal (dn/dt)
    dydt = zeros(size(y));
      da = p.da; % Diferencial de passo da discretização espacial
    
    % =====================================================================
    % ESQUEMA UPWIND (Diferenças Finitas de 1ª Ordem)
    % O fluxo numérico é retroativo (Nó(i) subtrai Nó(i-1)) para assegurar 
    % estabilidade direcional e evitar oscilações espúrias, uma vez que a 
    % advecção (envelhecimento) ocorre sempre num único sentido.
    % =====================================================================
    
    % --- EQUACIONAMENTO DO ESTÁGIO 1 (Ovos/L1 a L3) ---
    % Nó de Fronteira Inicial: Entrada do sistema = Ovos colocados
        dydt(p.idx_stg1(1)) = (Flux_in_stg1 - g_stg1*n1(1))/da - mu_stg1*n1(1);
    % Nós internos da malha do Estágio 1
    dydt(p.idx_stg1(2:end)) = -(g_stg1/da)*(n1(2:end) - n1(1:end-1)) - mu_stg1*n1(2:end);

    % --- EQUACIONAMENTO DO ESTÁGIO 2 (L4 a L7/Pupas) ---
    % Nó de Fronteira: O fluxo de entrada são as larvas maduras saindo do estágio 1
        dydt(p.idx_stg2(1)) = (g_stg1*n1(end) - g_stg2*n2(1))/da - mu_stg2*n2(1);
    % Nós internos da malha do Estágio 2
    dydt(p.idx_stg2(2:end)) = -(g_stg2/da)*(n2(2:end) - n2(1:end-1)) - mu_stg2*n2(2:end);

    % --- EQUACIONAMENTO DO ESTÁGIO 3 (Adultos) ---
    % Nó de Fronteira: O fluxo de entrada é a emergência das pupas do estágio 2
        dydt(p.idx_stg3(1)) = (g_stg2*n2(end) - g_stg3*n3(1))/da - mu_stg3_total_vec(1)*n3(1);
    % Nós internos da malha do Estágio 3 (Nota: A mortalidade aqui é matrizada pontualmente '.*')
    dydt(p.idx_stg3(2:end)) = -(g_stg3/da)*(n3(2:end) - n3(1:end-1)) - (mu_stg3_total_vec(2:end) .* n3(2:end));
end

% -------------------------------------------------------------------------
% FUNÇÕES BIOLÓGICAS AUXILIARES
% -------------------------------------------------------------------------

% Modelo de Brière para simulação de velocidade do metabolismo dependente da temperatura
function val = local_calc_briere(T, params)
         val = zeros(size(T));
        % Zera rigorosamente o fluxo caso a temperatura escape dos limites suportados
        mask = (T > params.Tmin) & (T < params.Tmax);
        if any(mask, 'all')
           T_m = T(mask);
           val(mask) = params.gamma .* T_m .* (T_m - params.Tmin) .* sqrt(abs(params.Tmax - T_m));
        end
end

% Função quadrática de mortalidade térmica para fases jovens
function mu = local_calc_mu_imature(T, coeffs)
         mu = max(0, polyval(coeffs, T));
end

% Função quadrática de mortalidade térmica específica para fase reprodutiva
function mu_therm = local_calc_mu_adult_thermal(T, p)
         mu_therm = max(0, polyval(p.mu_stg3_coeffs, T));
end

% A senescência cria um "paredão" de mortalidade ao fim da vida fisiológica, 
% convergindo artificialmente a densidade para zero no limiar do A3.
function mu_age = local_calc_senescence(age, p)
         mu_age = (p.mu_stg3_senescence .* age) ./ max(1e-6, (p.A3 - age).^2);
end

% -------------------------------------------------------------------------
% GERADOR DE CURVAS PARA O VETOR DE CONDIÇÃO INICIAL (t=0)
% -------------------------------------------------------------------------
function n0 = gaussian_dist(ages, A_max, Pop_Ini, da)
      sigma = A_max / 3;
    % Desvio base para a distribuição suave na malha
    
    % Distribui os indivíduos pelo espaço fisiológico
    n0 = exp(-((ages - ages(1)).^2)/(2*sigma^2));
    
    % Ajuste fino (Shift): Evita descontinuidades abruptas puxando a cauda para exatamente zero
            n0 = n0 - n0(end);
    n0(n0 < 0) = 0; % Salvaguarda numérica
    
    % Fator K de normalização da integral numérica (Garante que a área total seja = Pop_Ini)
    n0 = n0 * (Pop_Ini / (sum(n0)*da));
end

% -------------------------------------------------------------------------
% SIMULAÇÃO ISOLADA PARA PROCESSAMENTO EM PARALELO (PARFOR)
% -------------------------------------------------------------------------
% Encapsula a configuração do cenário, a resolução das EDOs (ode15s) e o 
% pós-processamento temporal. Isso evita vazamento de escopo, otimiza o uso 
% de memória nos workers do 'parfor' e retorna apenas os limites populacionais 
% (mín/máx) dentro da janela de análise.
function [t_var, y_var, min1, max1, min2, max2, min3, max3] = ...
    simulate_scenario(T_avg, Amp, p_local, time_sim, y_inicial, opts, day_start, day_end)
    
    % 1. Configura os parâmetros específicos do cenário
    p_local.T_average_C = T_avg;
    p_local.amplitude_C = Amp;
    
    % Construção da forçante térmica sazonal (Cosseno)
    p_local.T_func = @(t) p_local.T_average_C + p_local.amplitude_C * cos(2 * pi * t / p_local.period);
    
    % 2. Resolução do sistema de EDOs usando o ode15s
    [t_var, y_var] = ode15s(@(t,y) odefun_3stage(t, y, p_local), time_sim, y_inicial, opts);
    
    % 3. Integração numérica espacial para calcular a população total por estágio
    tot_stg1_var = sum(y_var(:, p_local.idx_stg1), 2) * p_local.da;
    tot_stg2_var = sum(y_var(:, p_local.idx_stg2), 2) * p_local.da;
    tot_stg3_var = sum(y_var(:, p_local.idx_stg3), 2) * p_local.da;
    
    % 4. Filtro temporal para recolha dos picos populacionais
    idx_analysis = (t_var >= day_start) & (t_var <= day_end);
    if ~any(idx_analysis)
       idx_analysis = true(size(t_var)); 
    end
    
    % 5. Captura do limite inferior e superior (mín e máx) de indivíduos
    min1 = min(tot_stg1_var(idx_analysis)); max1 = max(tot_stg1_var(idx_analysis));
    min2 = min(tot_stg2_var(idx_analysis)); max2 = max(tot_stg2_var(idx_analysis));
    min3 = min(tot_stg3_var(idx_analysis)); max3 = max(tot_stg3_var(idx_analysis));
end

% -------------------------------------------------------------------------
% UTILITÁRIO GRÁFICO GLOBAL (VERSÃO INTELIGENTE ANTI-CORTE)
% -------------------------------------------------------------------------
function set_ax(xlbl, ylbl, x_lim)
    xlabel(xlbl, 'FontSize', 12);
    ylabel(ylbl, 'FontSize', 12);
    
    if nargin < 3
        x_lim = xlim; 
    end
    
    % Calcula a folga (0.1% para a esquerda e 1% para a direita)
    Break1 = (x_lim(2) - x_lim(1)) * 0.001; % Folga da esquerda
    Break2 = (x_lim(2) - x_lim(1)) * 0.001;  % Folga da direita
    
    % Aplica Break1 no limite inferior (esquerda) e Break2 no limite superior (direita)
    xlim([x_lim(1) - Break1, x_lim(2) + Break2]); 
    
    grid on;
    set(gca, 'FontSize', 12, 'GridAlpha', 0.15, 'LineWidth', 1.1);
    set(gcf, 'Color', 'w');
end