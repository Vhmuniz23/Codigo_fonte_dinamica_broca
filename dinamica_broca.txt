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
%   n(t,a)  : Densidade de indivíduos (indivíduos por dia de idade)
%   t       : Tempo Cronológico (dias)
%   a       : Idade fisiológica (adimensional ou normalizada)
%   g(T)    : Taxa de advecção (velocidade de envelhecimento fisiológico 
%             dependente da temperatura T)
%   mu(T,a) : Taxa de mortalidade instantânea (depende da temperatura T e 
%             da idade a, representando estresse térmico e senescência)
%
% MÉTODO DE RESOLUÇÃO NUMÉRICA:
% Método das Linhas (Method of Lines - MOL). A dimensão espacial da idade 
% 'a' é discretizada em compartimentos (caixas numéricas) usando diferenças 
% finitas (esquema Upwind). Isso transforma a EDP num grande sistema acoplado 
% de Equações Diferenciais Ordinárias (EDOs), que é resolvido ao longo do 
% tempo 't' utilizando o integrador 'ode15s' (adequado para sistemas rígidos).
% =========================================================================

% [MAT] Limpeza completa do ambiente para evitar conflitos de memória
clearvars; clc; close all;  

% =========================================================================
% 1. PARAMETRIZAÇÃO BIOLÓGICA E AMBIENTAL DO MODELO
% =========================================================================

% --- 1.1 Forçante Climática (Sazonalidade) ---
p.periodo = 365;          % Período do ciclo térmico anual (dias)

% --- 1.2 Taxas de Desenvolvimento (Modelo Teórico de Brière-1) ---
% O modelo de Brière descreve o desenvolvimento assimétrico de insetos.
% a: constante empírica; Tmin: limiar térmico inferior; Tmax: limiar letal superior.
p.briere_egg = struct('a', 1.85e-4, 'Tmin', 12.8, 'Tmax', 36.8); % Ovos e L1
p.briere_larva = struct('a', 3.60e-5, 'Tmin', 14.3, 'Tmax', 36.3); % L2 a Pupas

% --- 1.3 Parâmetros de Reprodução / Fecundidade (Curva Gaussiana) ---
p.b_max = 93.63;   % Fecundidade máxima potencial (ovos/fêmea/dia) no pico térmico
p.T_opt_b = 24.27; % Temperatura ótima exata para a máxima oviposição (°C)
p.sigma_b = 4.46;  % Desvio padrão da curva (largura da janela térmica reprodutiva)

% --- 1.4 Parâmetros de Mortalidade ---
% Limites térmicos absolutos para sobrevivência (fora deste intervalo, a mortalidade é total)
p.mu_Tmax = 37.0;   
p.mu_Tmin = 12.0;   
p.mu_death = 15.0;  % Taxa letal padronizada (usada como teto para evitar infinitos numéricos)

% Imaturos: Mortalidade térmica modelada por um polinômio de 2º grau (Formato em 'U')
p.mu_im_coeffs = [0.0002, -0.0104, 0.1362]; 

% Adultos: A mortalidade é dividida em dois componentes: Térmico e Senescência
p.mu_ad_coeffs = [0.011, -0.572, 7.586];    % Coeficientes do estresse térmico
p.mu_ad_senescencia = 0.748123;             % Fator de peso para a morte por velhice

% --- 1.5 Competição Intraespecífica (Densidade-Dependência) ---
% Controla o crescimento populacional explosivo (Modelo tipo Beverton-Holt)
p.q_b = 2e-5;  % Intensidade da competição entre adultos por locais de oviposição

% =========================================================================
% 2. DISCRETIZAÇÃO DA MALHA (GRID) PARA O MÉTODO DAS LINHAS (MOL)
% =========================================================================
% Define-se uma temperatura de referência padrão para normalizar o envelhecimento
T_ref = 25.0; 
p.rate_egg_ref = local_calc_briere(T_ref, p.briere_egg);
p.rate_larva_ref = local_calc_briere(T_ref, p.briere_larva);

% Domínios de Idade Fisiológica Máxima de cada estágio na T_ref
p.A1 = 1 / p.rate_egg_ref;      % Idade fisiológica onde o inseto muda de Ovos para L2
p.A2 = 1 / p.rate_larva_ref;    % Idade fisiologica onde muda de Pupa para Adulto
p.A3 = 10.0;                    % Longevidade fisiológica máxima do estágio adulto

% Construção da Malha Espacial de Idade
p.da = p.A1 / 80;                   % Passo de integração 'da' (Tamanho de cada caixa etária)
p.ages_1 = (p.da/2 : p.da : p.A1)'; % Vetor de idades do Estágio 1 (Ponto central das caixas)
p.ages_2 = (p.da/2 : p.da : p.A2)'; % Vetor de idades do Estágio 2
p.ages_3 = (p.da/2 : p.da : p.A3)'; % Vetor de idades do Estágio 3

% Número total de compartimentos em cada estágio
p.N1 = length(p.ages_1); 
p.N2 = length(p.ages_2); 
p.N3 = length(p.ages_3);

% [MAT] OTIMIZAÇÃO DE PERFORMANCE: Pré-alocação dos índices vetoriais para o integrador
p.idx_O = 1 : p.N1;
p.idx_L = (p.N1 + 1) : (p.N1 + p.N2);
p.idx_A = (p.N1 + p.N2 + 1) : (p.N1 + p.N2 + p.N3);

% Pré-calcula a curva estática de envelhecimento natural (senescência) dos adultos
p.mu_ad_age_vec = local_calc_senescence(p.ages_3, p);

% Vetores deslocados linearmente (apenas para plotagem contínua no eixo X dos gráficos)
idades_plot_O = p.ages_1;
idades_plot_L = p.ages_2 + p.A1;      
idades_plot_A = p.ages_3 + p.A1 + p.A2; 

% =========================================================================
% 3. CONDIÇÕES INICIAIS DO SISTEMA NUMÉRICO (t = 0)
% =========================================================================
% Utiliza-se uma distribuição gaussiana deslocada para representar a entrada gradual 
% da população inicial no sistema, garantindo estabilidade no esquema numérico.
n0_O = gaussian_dist(p.ages_1, p.A1, 25, p.da);  % Condição inicial: 25 Ovos
n0_L = gaussian_dist(p.ages_2, p.A2, 45, p.da);  % Condição inicial: 45 Larvas
n0_A = gaussian_dist(p.ages_3, p.A3, 30, p.da);  % Condição inicial: 30 Adultos

% Vetor de estado inicial global (Concatenação de todos os estágios)
y0 = [n0_O; n0_L; n0_A];

% Paleta de cores global para padronização das figuras
colors = [0, 1, 0.498; 1, 0.871, 0.0129; 1, 0.173, 0.173; 0.6 0.2 0.8]; 
color_temp = [0.8500 0.3250 0.0980];
c_vals = {'b', 'r'};

% =========================================================================
% 4. ANÁLISE TEÓRICA DAS FUNÇÕES BIOLÓGICAS (Gráficos Preliminares 1 a 3)
% =========================================================================
disp('Gerando Gráficos das Funções Biológicas Globais (10°C a 40°C)...');

T_analise = linspace(10, 40, 500); % Malha fina de temperatura para curvas contínuas

% 4.1 Calcula as taxas relativas de advecção (Brière)
g_egg = local_calc_briere(T_analise, p.briere_egg) / p.rate_egg_ref;
g_larva = local_calc_briere(T_analise, p.briere_larva) / p.rate_larva_ref;

% 4.2 Fecundidade (com e sem penalização por competição de Beverton-Holt)
fec_pura = p.b_max * exp(-0.5 * ((T_analise - p.T_opt_b)/p.sigma_b).^2);
N_competicao = 1e+5; % Densidade arbitrária apenas para teste da curva
fec_comp = fec_pura ./ (1 + p.q_b * N_competicao);

% 4.3 Funções de penalidade por mortalidade térmica
mu_im_analise = local_calc_mu_imature(T_analise, p);
mu_ad_termica_analise = local_calc_mu_adult_thermal(T_analise, p);

% 4.4 Superfície 2D de mortalidade do adulto (Combina Temperatura e Idade)
idades_analise = linspace(0, p.A3-0.1, 300);
mu_ad_idade_analise = local_calc_senescence(idades_analise, p);
[TX, AX] = meshgrid(T_analise, idades_analise);
mu_total_grid = local_calc_mu_adult_thermal(TX, p) .* local_calc_senescence(AX, p); 

% --- PLOT FIG 1: Curvas de Brière (Envelhecimento) ---
figure(1); clf; set(gcf, 'Name', 'Desenvolvimento (Brière)');
g_vals = {g_egg, g_larva}; titles_f1 = {'Estágio 1 - Ovos + L1', 'Estágio 2 - L2 + Pupa'};
for i=1:2
    subplot(1,2,i); plot(T_analise, g_vals{i}, c_vals{i}, 'LineWidth', 2);
    title(titles_f1{i}); set_ax('Temperatura °C', sprintf('Taxa Relativa (g_%d) (1/dia)', i), [10 40]);
end

% --- PLOT FIG 2: Distribuição Gaussiana de Fecundidade ---
figure(2); clf; set(gcf, 'Name', 'Fecundidade (Oviposição)');
f_vals = {fec_pura, fec_comp}; 
titles_f2 = {'Fecundidade s/ Competição', sprintf('Fecundidade c/ Competição (N_3=%d)', N_competicao)};
y_labels_f2 = {'\beta_{max} (Ovos/Fêmea/Dia)', ''}; 
for i=1:2
    subplot(1,2,i); plot(T_analise, f_vals{i}, c_vals{i}, 'LineWidth', 2);
    title(titles_f2{i}); set_ax('Temperatura °C', y_labels_f2{i}, [10 40]);
end

% --- PLOT FIG 3: Superfícies e Curvas de Decaimento (Mortalidade) ---
figure(3); clf; set(gcf, 'Name', 'Mortalidade Imaturos e Adultos');
subplot(2,2,1); plot(T_analise, mu_im_analise, 'b', 'LineWidth', 2); title('Mortalidade Imaturos (Térmica)'); set_ax('Temperatura °C', 'Taxa Térmica - \mu_{1,2} (1/dia)', [10 40]);
subplot(2,2,2); plot(T_analise, mu_ad_termica_analise, 'r', 'LineWidth', 2); title('Mortalidade Adultos (Térmica)'); set_ax('Temperatura °C', 'Taxa Térmica - \mu_3 (1/dia)', [10 40]);
subplot(2,2,3); plot(idades_analise, mu_ad_idade_analise, 'Color', colors(4,:), 'LineWidth', 2); title('Mortalidade Adultos (Senescência)'); set_ax('Idade Fisiológica (dias)', 'Taxa Etária - \mu_3 (1/dia)', [0 p.A3]);
subplot(2,2,4); surf(TX, AX, mu_total_grid, 'EdgeColor', 'none'); shading interp; view(3); colormap('jet'); colorbar; 
title('Mortalidade Total Adultos'); set_ax('Temperatura °C', 'Idade Fisiológica (dias)', [10 40]); zlabel('Taxa Total - \mu_3 (1/dia)'); ylim([0 p.A3]);

% =========================================================================
% 5. SIMULAÇÃO EM LOTE (SWEEP): VARREDURA CLIMÁTICA E SAZONAL
% =========================================================================
% Esta seção resolve o modelo para múltiplas combinações de temperatura média 
% e amplitude sazonal, extraindo o potencial de infestação para cada cenário.
disp('Iniciando simulações de todos os cenários climáticos (ODE15s)...');

T_medias_teste = [12, 16, 20, 24, 28, 32, 36, 40]; % Malha de Temperaturas Médias Anuais
amplitudes_teste = [0, 3, 6, 9];                   % Malha de Amplitudes Sazonais (Flutuação)
tempo_simulacao = [0 1278];                        % Janela de integração (t0 a tf)
dia_inicio_analise = 180;                          % Início do período para extrair picos (ignora transiente)
dia_fim_analise = 730; 

% Criação do espaço de combinações (Grid paramétrica)
[T_grid, Amp_grid] = ndgrid(T_medias_teste, amplitudes_teste);
T_vec = T_grid(:); Amp_vec = Amp_grid(:);
N_cenarios = length(T_vec);

% Células e matrizes para armazenamento dos resultados
t_all = cell(N_cenarios, 1); y_all = cell(N_cenarios, 1);
Pop_Min_O = zeros(N_cenarios, 1); Pop_Max_O = zeros(N_cenarios, 1);
Pop_Min_L = zeros(N_cenarios, 1); Pop_Max_L = zeros(N_cenarios, 1);
Pop_Min_A = zeros(N_cenarios, 1); Pop_Max_A = zeros(N_cenarios, 1);

% Configurações do Integrador Numérico (Garante que populações não fiquem negativas)
opts_sweep = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NonNegative', 1:length(y0));

tic; 
% Laço paralelo (parfor) acelera enormemente as simulações da varredura
parfor k = 1:N_cenarios 
    p_local = p;
    p_local.T_media_C = T_vec(k);
    p_local.amplitude_C = Amp_vec(k);
    
    % Construção da forçante térmica sazonal (Cosseno)
    p_local.T_func = @(t) p_local.T_media_C + p_local.amplitude_C * cos(2 * pi * t / p_local.periodo);
    
    % Resolução do sistema de EDOs usando o ode15s
    [t_var, y_var] = ode15s(@(t,y) odefun_3stage(t, y, p_local), tempo_simulacao, y0, opts_sweep);
    
    t_all{k} = t_var; y_all{k} = y_var;
    
    % Integração numérica espacial para calcular a população total por estágio
    tot_O_var = sum(y_var(:, p_local.idx_O), 2) * p_local.da;
    tot_L_var = sum(y_var(:, p_local.idx_L), 2) * p_local.da;
    tot_A_var = sum(y_var(:, p_local.idx_A), 2) * p_local.da;
    
    % Filtro temporal para recolha dos picos populacionais
    idx_analise = (t_var >= dia_inicio_analise) & (t_var <= dia_fim_analise);
    if ~any(idx_analise), idx_analise = true(size(t_var)); end
    
    % Captura do limite inferior e superior (mín e máx) de indivíduos
    Pop_Min_O(k) = min(tot_O_var(idx_analise)); Pop_Max_O(k) = max(tot_O_var(idx_analise));
    Pop_Min_L(k) = min(tot_L_var(idx_analise)); Pop_Max_L(k) = max(tot_L_var(idx_analise));
    Pop_Min_A(k) = min(tot_A_var(idx_analise)); Pop_Max_A(k) = max(tot_A_var(idx_analise));
end
fprintf('Simulações Concluídas em %.2f segundos.\n', toc);

% Reconstrução das matrizes 2D (Média x Amplitude) para os gráficos de varredura
Mat_Min_O = reshape(Pop_Min_O, length(T_medias_teste), length(amplitudes_teste));
Mat_Max_O = reshape(Pop_Max_O, length(T_medias_teste), length(amplitudes_teste));
Mat_Min_L = reshape(Pop_Min_L, length(T_medias_teste), length(amplitudes_teste));
Mat_Max_L = reshape(Pop_Max_L, length(T_medias_teste), length(amplitudes_teste));
Mat_Min_A = reshape(Pop_Min_A, length(T_medias_teste), length(amplitudes_teste));
Mat_Max_A = reshape(Pop_Max_A, length(T_medias_teste), length(amplitudes_teste));

% =========================================================================
% 6. GRÁFICOS DE VARREDURA (BIFURCAÇÃO ECOLÓGICA / BOXPLOTS)
% =========================================================================
% Exporta resultados consolidados do sweep para análise externa
Tabela_Varredura = table(T_vec, Amp_vec, Pop_Min_O, Pop_Max_O, Pop_Min_L, Pop_Max_L, Pop_Min_A, Pop_Max_A);
writetable(Tabela_Varredura, 'Varredura_Temp_Amp_Broca.xlsx');

titulos_var = {'Estágio 1 - Ovos + L1', 'Estágio 2 - L2 + Pupas', 'Estágio 3 - Adultos (Machos + Fêmeas)'};
Dados_Min = {Mat_Min_O, Mat_Min_L, Mat_Min_A}; Dados_Max = {Mat_Max_O, Mat_Max_L, Mat_Max_A};
str_intervalo = sprintf('(Dias %g a %g)', dia_inicio_analise, dia_fim_analise);
offset_base = 0.2; w = 0.08; posicoes_y = 1:length(T_medias_teste);

% Gráficos estilo "Tornado" ou boxplot horizontal demonstrando os mín e máx populacionais
for k = 1:3
    figure(3+k); clf; set(gcf, 'Name', ['Amplitude Populacional - ' titulos_var{k}], 'Color', 'w'); hold on;
    d_min = round(Dados_Min{k}); d_min(d_min < 0) = 0;
    d_max = round(Dados_Max{k}); d_max(d_max < 0) = 0;
    max_x_geral = max(10, max(max(d_max)));
    h_leg = zeros(1, length(amplitudes_teste)); legend_strs = cell(1, length(amplitudes_teste));
    
    for j = 1:length(amplitudes_teste)
        y_plot = posicoes_y + ((j - (length(amplitudes_teste)+1)/2) * offset_base);
        for i = 1:length(T_medias_teste)
            x_min = d_min(i, j); x_max = d_max(i, j); y_c = y_plot(i);
            % Desenha a barra horizontal de amplitude da população
            patch([x_min, x_min, x_max, x_max], [y_c - w, y_c + w, y_c + w, y_c - w], colors(j,:), 'FaceAlpha', 0.65, 'EdgeColor', colors(j,:), 'LineWidth', 1.0);
            texto_val = sprintf('%d - %d', x_min, x_max); if x_min == x_max, texto_val = sprintf('%d', x_max); end
            text(x_max + (max_x_geral * 0.015), y_c, texto_val, 'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 8.5, 'Color', 'k');
        end
        h_leg(j) = patch([NaN NaN NaN NaN], [NaN NaN NaN NaN], colors(j,:), 'FaceAlpha', 0.65, 'EdgeColor', colors(j,:));
        legend_strs{j} = sprintf('Amp %d °C', amplitudes_teste(j));
    end
    title(sprintf('%s\nAmplitude Populacional %s', titulos_var{k}, str_intervalo), 'FontWeight', 'bold', 'FontSize', 12);
    xlabel('Nº de Indivíduos (Mín a Máx)'); ylabel('Temperatura Média Anual (°C)'); 
    ylim([0.5, length(T_medias_teste) + 0.5]); xlim([0, max_x_geral * 1.25]); 
    yticks(posicoes_y); yticklabels(string(T_medias_teste)); 
    for y_sep = 1.5 : 1 : (length(T_medias_teste) - 0.5), yline(y_sep, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1); end
    grid on; box on; legend(h_leg, legend_strs, 'Location', 'northeast'); hold off;
end

% =========================================================================
% 7. PAINEL DE CONTROLE: GERAÇÃO DE GRÁFICOS DO CENÁRIO ESPECÍFICO
% =========================================================================
% A partir deste ponto, o código foca na dinâmica temporal de um ÚNICO cenário.
% => ESCOLHA AQUI O CENÁRIO QUE DESEJA PLOTAR (DEVE EXISTIR NA VARREDURA):
T_escolhida = 32;   
Amp_escolhida = 9;  
disp(['Extraindo e plotando dados completos para T = ' num2str(T_escolhida) '°C e Amp = ' num2str(Amp_escolhida) '°C...']);

% Procura o cenário na memória
idx_plot = find(T_vec == T_escolhida & Amp_vec == Amp_escolhida);
if isempty(idx_plot), error('O cenário escolhido não foi simulado na varredura.'); end

% Extração das variáveis temporais do cenário escolhido
t = t_all{idx_plot};
y_sol = y_all{idx_plot};
p.T_media_C = T_escolhida; p.amplitude_C = Amp_escolhida;
p.T_func = @(t) p.T_media_C + p.amplitude_C * cos(2 * pi * t / p.periodo);

% Separação espacial das densidades e cálculo dos totais por estágio
n_sol_O = y_sol(:, p.idx_O); n_sol_L = y_sol(:, p.idx_L); n_sol_A = y_sol(:, p.idx_A);
total_O = sum(n_sol_O, 2)*p.da; total_L = sum(n_sol_L, 2)*p.da; total_A = sum(n_sol_A, 2)*p.da;

% --- RECONSTRUÇÃO DOS FLUXOS INTERNOS (Para fins de gráficos) ---
temp_t_C = p.T_func(t); 
v_taxa_g_L_eff = local_calc_briere(temp_t_C, p.briere_larva) / p.rate_larva_ref; 
mu_im_base = local_calc_mu_imature(temp_t_C, p);
mu_ad_thermal_plot = max(1e-4, local_calc_mu_adult_thermal(temp_t_C, p));

% Reconstrução da dinâmica de reprodução
v_taxa_b_potencial = p.b_max * exp(-0.5 * ((temp_t_C - p.T_opt_b)/p.sigma_b).^2); 
v_taxa_b_real = v_taxa_b_potencial ./ (1 + p.q_b .* total_A);

% Reconstrução dos fluxos de transição entre fases
v_fluxo_estabelecimento = (local_calc_briere(temp_t_C, p.briere_egg)/p.rate_egg_ref) .* n_sol_O(:, end);
v_fluxo_emergencia = v_taxa_g_L_eff .* n_sol_L(:, end);

% Reconstrução do status demográfico do adulto
n_sol_A_safe = max(0, n_sol_A); total_adultos_safe = max(1, sum(n_sol_A_safe, 2));
idade_media_t = (n_sol_A_safe * p.ages_3) ./ total_adultos_safe;
mu_senescencia_media = (n_sol_A_safe * p.mu_ad_age_vec) ./ total_adultos_safe;
mu_ad_vec_total = min(p.mu_death, mu_ad_thermal_plot .* mu_senescencia_media);

% --- FLAG DE CONTROLE GRÁFICO (Sazonal vs Constante) ---
is_temp_const = max(temp_t_C) == min(temp_t_C);
if is_temp_const, str_temp = sprintf(' (Temp: %.1f °C)', temp_t_C(1)); else, str_temp = ''; end

% =========================================================================
% --- FIGURA 7: POTENCIAL TÉRMICO E ABUNDÂNCIA MÁXIMA GLOBAL ---
% =========================================================================
% Une as curvas teóricas (Brière e Fecundidade) com os picos reais 
% simulados gerados pela matriz de sweep. Mostra a "Janela de Viabilidade".
figure(7); clf; set(gcf, 'Name', 'Potencial Térmico e Abundância Máxima');

% Achatamento das matrizes com (:) encontra o maior pico global simulado
true_max_pops = [max(Mat_Max_O(:)), max(Mat_Max_L(:)), max(Mat_Max_A(:))];
curves_f4 = {local_calc_briere(T_analise, p.briere_egg), local_calc_briere(T_analise, p.briere_larva)};
T_params = {p.briere_egg, p.briere_larva}; 
titles_f4 = {'Estágio 1 - Ovos + L1', 'Estágio 2 - L2 + Pupas'};

for i=1:2
    subplot(1,3,i); 
    curve = (curves_f4{i} / max(curves_f4{i})) * true_max_pops(i);
    [~, idx] = max(curve); 
    
    plot(T_analise, curve, 'Color', colors(i,:), 'LineWidth', 2); hold on;
    xline(T_params{i}.Tmin, '--k', 'Label', 'T_{min}'); 
    xline(T_analise(idx), '--k', 'Label', sprintf('T_{opt} = %.1f °C', T_analise(idx)), ...
          'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'center'); 
    xline(T_params{i}.Tmax, '--k', 'Label', 'T_{max}');
    
    title(titles_f4{i}, 'FontSize', 10); 
    if i == 1, set_ax('Temperatura °C', 'Potencial Máximo (Indivíduos)', [10 40]); 
    else, set_ax('Temperatura °C', '', [10 40]); end
    if true_max_pops(i) > 0, ylim([0 true_max_pops(i)*1.15]); end
    hold off;
end

subplot(1,3,3);
% Suavização das extremidades da Gaussiana de adultos para criar forma de sino perfeita
y_at_min = p.b_max * exp(-0.5 * ((p.mu_Tmin - p.T_opt_b)/p.sigma_b).^2); 
y_at_max = p.b_max * exp(-0.5 * ((p.mu_Tmax - p.T_opt_b)/p.sigma_b).^2);
curve_3_base = fec_pura - (y_at_min + ((y_at_max - y_at_min)/(p.mu_Tmax - p.mu_Tmin)) * (T_analise - p.mu_Tmin));
curve_3_base(T_analise <= p.mu_Tmin | T_analise >= p.mu_Tmax | curve_3_base < 0) = 0;

if max(curve_3_base) > 0, curve_3 = (curve_3_base / max(curve_3_base)) * true_max_pops(3);
else, curve_3 = curve_3_base; end

[~, idx3] = max(curve_3); 
plot(T_analise, curve_3, 'Color', colors(3,:), 'LineWidth', 2); hold on;
xline(p.mu_Tmin, '--k', 'Label', 'T_{min}'); 
xline(T_analise(idx3), '--k', 'Label', sprintf('T_{opt} = %.1f °C', T_analise(idx3)), ...
      'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'center'); 
xline(p.mu_Tmax, '--k', 'Label', 'T_{max}');

title('Estágio 3 - Adultos', 'FontSize', 10); set_ax('Temperatura °C', '', [10 40]); 
if true_max_pops(3) > 0, ylim([0 true_max_pops(3)*1.15]); end
hold off;

% =========================================================================
% --- FIGS 8 e 9: CONDIÇÃO INICIAL E EVOLUÇÃO TEMPORAL (Densidades) ---
% =========================================================================
figure(8); clf; set(gcf, 'Name', 'Distribuição Etária Inicial n(0,a)', 'Color', 'w');
figure(9); clf; set(gcf, 'Name', 'Evolução Populacional Temporal', 'Color', 'w');

n0_all = {n0_O, n0_L, n0_A}; 
idades_all = {idades_plot_O, idades_plot_L, idades_plot_A}; 
totals = {total_O, total_L, total_A}; 
titles_f5 = {'Cond. Inicial: Ovos + L1', 'Cond. Inicial: L2 + Pupas', 'Cond. Inicial: Adultos'}; 
titles_f6 = {'Evolução Estágio 1 - Ovos + L1', 'Evolução Estágio 2 - L2 + Pupas', 'Evolução Estágio 3 - Adultos'};

global_x_max = max([max(idades_plot_O), max(idades_plot_L), max(idades_plot_A)]);

for i = 1:3
    % --- FIGURA 8 (Distribuições t=0) ---
    figure(8); subplot(1,3,i);
    x_stage = idades_all{i}(:); y_stage = n0_all{i}(:);
    
    % Plota estritamente os limites originais do estágio
    plot(x_stage, y_stage, 'Color', colors(i,:), 'LineWidth', 2);
    
    title(titles_f5{i}); set_ax('Idade Fisiológica (dias)', sprintf('Densidade n_%d(0,a)', i));
    
    % Trava o eixo X exatamente na primeira e na última idade da caixa
    xlim([x_stage(1), x_stage(end)]);
    
    y_max = max(y_stage);
    if y_max > 0, ylim([0, y_max * 1.15]); end
    grid on;

    % --- FIGURA 9 (Curvas de crescimento) ---
    figure(9); subplot(3,1,i);
    if is_temp_const
        plot(t, totals{i}, 'Color', colors(i,:), 'LineWidth', 2);
        title([titles_f6{i} str_temp]); xlim([min(t), max(t)]); ylabel('Total de Indivíduos'); grid on;
    else
        % Em caso de sazonalidade, cria duplo eixo (yyaxis)
        yyaxis left; plot(t, totals{i}, 'Color', colors(i,:), 'LineWidth', 2);
        ylabel('Total de Indivíduos'); title(titles_f6{i}); xlim([min(t), max(t)]); grid on;
        yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp, 'LineWidth', 0.8);
        ylabel('Temp. °C'); set(gca, 'YColor', color_temp);
    end
    if i == 3, xlabel('Tempo Cronológico (dias)'); end
end

% --- FIGS 10 a 12: SUPERFÍCIES 3D (ESPAÇO TEMPO-IDADE) ---
% Mostra graficamente a advecção (movimento da população pela matriz etária)
idx_plot_3d = 1:max(1, floor(length(t)/150)):length(t); t_plot_3d = t(idx_plot_3d);
sols = {n_sol_O, n_sol_L, n_sol_A}; figs3D = {'Estágio 1 - Ovos + L1', 'Estágio 2 - L2 + Pupas', 'Estágio 3 - Adultos (Machos + Fêmeas)'};
for i=1:3
    figure(9+i); clf; set(gcf, 'Name', ['Evolução Etária 3D: ' figs3D{i}]); 
    surf(idades_all{i}, t_plot_3d, sols{i}(idx_plot_3d, :), 'EdgeColor', 'none'); 
    shading interp; colormap('jet'); colorbar; view(3); title(titles_f6{i}); set_ax('Idade Fisiológica (dias)', 'Tempo Cronológico (dias)'); zlabel('Densidade');
end

% --- FIG 13: FORÇANTE CLIMÁTICA (Termômetro do Sistema) ---
figure(13); clf; set(gcf, 'Name', 'Monitoramento Climático'); 
plot(t, temp_t_C, 'Color', color_temp, 'LineWidth', 2.0); 
title('Função Temperatura (Sazonalidade)'); set_ax('Tempo Cronológico (dias)', 'Temperatura °C'); yline(p.T_media_C, '--k');

% --- FIGS 14 e 15: FLUXOS DEMOGRÁFICOS DE CONTORNO (Limites do Domínio) ---
figure(14); clf; set(gcf, 'Name', 'Fluxos Demográficos');
y_vals_f11 = {v_fluxo_estabelecimento, v_fluxo_emergencia}; y_lbls_f11 = {'Larvas/dia', 'Adultos/dia'};
titles_f11 = {'Estabelecimento (Entrada no colmo)', 'Emergência (Novos Adultos)'};
for i=1:2
    subplot(2,1,i);
    if is_temp_const
        plot(t, y_vals_f11{i}, c_vals{i}, 'LineWidth', 2); ylabel(y_lbls_f11{i}); title([titles_f11{i} str_temp]); grid on;
    else
        yyaxis left; plot(t, y_vals_f11{i}, c_vals{i}, 'LineWidth', 2); ylabel(y_lbls_f11{i}); title(titles_f11{i}); grid on;
        yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp, 'LineWidth', 0.8); ylabel('Temperatura °C');
    end
    if i==2, xlabel('Tempo Cronológico (dias)'); end
end

figure(15); clf; set(gcf, 'Name', 'Fecundidade e Competição');
if is_temp_const
    plot(t, v_taxa_b_potencial, [c_vals{1} '-'], 'LineWidth', 2); hold on; plot(t, v_taxa_b_real, [c_vals{2} '--'], 'LineWidth', 2); hold off;
    ylabel('\beta_{max} (Ovos/Fêmea/Dia)'); title(['Fecundidade s/ Competição vs c/ Competição' str_temp]); legend('Fecundidade s/ Comp.', 'Fecundidade c/ Comp.', 'Location', 'northeast');
else
    yyaxis left; plot(t, v_taxa_b_potencial, [c_vals{1} '-'], 'LineWidth', 2); hold on; plot(t, v_taxa_b_real, [c_vals{2} '--'], 'LineWidth', 2); hold off; ylabel('\beta_{max} (Ovos/Fêmea/Dia)');
    yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp, 'LineWidth', 0.8); ylabel('Temperatura °C');
    title('Fecundidade s/ Competição vs c/ Competição'); legend('Fecundidade s/ Comp.', 'Fecundidade c/ Comp.', 'Temperatura', 'Location', 'northeast');
end
xlabel('Tempo Cronológico (dias)'); grid on;

% --- FIG 16: DINÂMICA DAS TAXAS DE MORTALIDADE TEMPORAIS ---
figure(16); clf; set(gcf, 'Name', 'Análise de Mortalidade', 'Color', 'w');
y_vals_f13 = {mu_im_base, mu_ad_thermal_plot, mu_senescencia_media, mu_ad_vec_total};
titles_f13 = {'Mortalidade Imaturos (Térmica)', 'Mortalidade Adultos (Térmica)', 'Mortalidade Adultos (Senescência)', 'Mortalidade Total Adultos'};
y_lbls_f13 = {'Taxa Térmica - \mu_{1,2} (1/dia)', 'Taxa Térmica - \mu_3 (1/dia)', 'Taxa Etária - \mu_3 (1/dia)', 'Taxa Total - \mu_3 (1/dia)'}; 
c_f13 = {colors(1,:), colors(2,:), colors(4,:), colors(3,:)};

for i=1:4
    subplot(2,2,i); 
    
    % Define a legenda do eixo X apenas para os gráficos da linha de baixo (3 e 4)
    if i >= 3
        x_lbl = 'Tempo Cronológico (dias)';
    else
        x_lbl = '';
    end
    
    if is_temp_const 
        plot(t, y_vals_f13{i}, 'Color', c_f13{i}, 'LineWidth', 2); title([titles_f13{i} str_temp]); set_ax(x_lbl, y_lbls_f13{i});
    else
        yyaxis left; plot(t, y_vals_f13{i}, 'Color', c_f13{i}, 'LineWidth', 2); title(titles_f13{i}); set_ax(x_lbl, y_lbls_f13{i});
        yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp, 'LineWidth', 0.8); ylabel('Temperatura °C'); set(gca, 'YColor', color_temp);
    end
end

% --- FIG 17: PERFIL ETÁRIO GERAL (Agrupamento Histórico) ---
N_total_age = trapz(t, y_sol, 1)'; 
x_cont = [idades_plot_O; idades_plot_L; idades_plot_A];

janela_suavizacao = round(length(N_total_age) * 0.10);
y_smooth = smoothdata(N_total_age, 'gaussian', janela_suavizacao);
n_tail = round(length(idades_plot_A) * 0.40); 
y_smooth(end-n_tail+1 : end) = y_smooth(end-n_tail+1 : end) .* ((cos(linspace(0, pi, n_tail)') + 1) / 2);
idxs = {1:p.N1, p.N1:(p.N1 + p.N2), (p.N1 + p.N2):length(x_cont)};

figure(17); clf; set(gcf, 'Name', 'Perfil Etário Total', 'Color', 'w'); hold on;
h = zeros(1,3);
for i=1:3
    x_poly = [x_cont(idxs{i})', fliplr(x_cont(idxs{i})')]; y_poly = [y_smooth(idxs{i})', zeros(1, length(idxs{i}))];
    fill(x_poly, y_poly, colors(i,:), 'EdgeColor', colors(i,:), 'LineWidth', 2.5, 'FaceAlpha', 0.25);
    plot(x_cont(idxs{i}), y_smooth(idxs{i}), 'Color', colors(i,:), 'LineWidth', 3);
    h(i) = fill(NaN, NaN, colors(i,:), 'FaceAlpha', 0.25, 'EdgeColor', colors(i,:), 'LineWidth', 2.5);
end
box on; set(gca, 'TickDir', 'in'); set_ax('Idade Fisiológica (dias)', 'Nº Total de Indivíduos', [0, p.A1 + p.A2 + p.A3]); ylim([0, max(y_smooth) * 1.1]);
title('Estrutura Etária Total da Geração por Estágio'); legend(h, {'Estágio 1 - Ovos + L1', 'Estágio 2 - L2 + Pupas', 'Estágio 3 - Adultos'}, 'Location', 'northeast'); hold off;

% =========================================================================
% RELATÓRIOS DE CONSOLE E EXPORTAÇÃO DO CENÁRIO
% =========================================================================
% Imprime no terminal do MATLAB o Paradoxo de Malha e a real T_opt do modelo

fprintf('\n=======================================================\n');
fprintf('   RESULTADOS DA VARREDURA TÉRMICA E POTENCIAL TEÓRICO \n');
fprintf('=======================================================\n');

dados_sweep = {Mat_Max_O, Mat_Max_L, Mat_Max_A};
nomes_estagios = {'Ovos + L1', 'L2 + Pupas', 'Adultos'};

c1 = local_calc_briere(T_analise, p.briere_egg);
c2 = local_calc_briere(T_analise, p.briere_larva);
c3 = curve_3_base; % Reaproveita a curva limpa sem bicos gerada na Fig 7
curvas_teoricas = {c1, c2, c3};

for i = 1:3
    % Encontra qual coordenada (linha/coluna) da matriz de resultados contém o valor máximo
    [max_pop, linear_idx] = max(dados_sweep{i}(:));
    [row, col] = ind2sub(size(dados_sweep{i}), linear_idx);
    T_simulada = T_medias_teste(row);
    Amp_simulada = amplitudes_teste(col);
    
    % Encontra a T_opt teórica analítica (matriz contínua fina de teste)
    [~, idx_opt] = max(curvas_teoricas{i});
    T_opt_teorica = T_analise(idx_opt);
    
    fprintf('ESTÁGIO %d: %s\n', i, nomes_estagios{i});
    fprintf('  -> Temperatura Ótima Verdadeira (T_opt): %.1f °C\n', T_opt_teorica);
    fprintf('  -> Maior População Encontrada: %.0f indivíduos\n', max_pop);
    fprintf('     (Gerada no cenário simulado mais próximo: Temp = %d °C | Amp = %d °C)\n', T_simulada, Amp_simulada);
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
    n1 = y(p.idx_O); % Ovos + L1
    n2 = y(p.idx_L); % L2 + Pupas
    n3 = y(p.idx_A); % Adultos
    
    % Cálculo populacional macro (Integral da densidade em relação à idade 'a')
    N3_total = sum(n3)*p.da; 
    
    % 3. Taxas Direcionais de Advecção (Velocidade na malha)
    % Velocidade normalizada na temperatura padrão (Envelhecimento térmico)
    g1 = local_calc_briere(T, p.briere_egg) / p.rate_egg_ref;
    g2 = local_calc_briere(T, p.briere_larva) / p.rate_larva_ref;
    g3 = 1.0; % Adultos experimentam 1 dia fisiológico a cada 1 dia cronológico
    
    % 4. Módulo de Mortalidade / Perdas (Decaimento populacional contínuo)
    mu_im = local_calc_mu_imature(T, p); 
    % A mortalidade do adulto flutua individualmente pelo vetor de idades (senescência)
    mu_ad_total_vec = local_calc_mu_adult_thermal(T, p) * p.mu_ad_age_vec; 
    
    % 5. Condição de Contorno Superior (Recrutamento no nó 1 do estágio 1)
    fec_real = (p.b_max * exp(-0.5 * ((T - p.T_opt_b)/p.sigma_b)^2)) / (1 + p.q_b * N3_total);
    Flux_in_1 = fec_real * N3_total; 
    
    % Inicializa o vetor espelho de taxas de variação temporal (dn/dt)
    dydt = zeros(size(y)); 
    da = p.da; % Diferencial de passo da discretização espacial
    
    % =====================================================================
    % ESQUEMA UPWIND (Diferenças Finitas de 1ª Ordem)
    % O fluxo numérico é retroativo (Nó(i) subtrai Nó(i-1)) para assegurar 
    % estabilidade direcional e evitar oscilações espúrias, uma vez que a 
    % advecção (envelhecimento) ocorre sempre num único sentido.
    % =====================================================================
    
    % --- EQUACIONAMENTO DO ESTÁGIO 1 (Ovos + L1) ---
    % Nó de Fronteira Inicial: Entrada do sistema = Ovos colocados
    dydt(p.idx_O(1)) = (Flux_in_1 - g1*n1(1))/da - mu_im*n1(1); 
    % Nós internos da malha do Estágio 1
    dydt(p.idx_O(2:end)) = -(g1/da)*(n1(2:end) - n1(1:end-1)) - mu_im*n1(2:end); 
    
    % --- EQUACIONAMENTO DO ESTÁGIO 2 (L2 + Pupas) ---
    % Nó de Fronteira: O fluxo de entrada são as larvas maduras saindo do estágio 1
    dydt(p.idx_L(1)) = (g1*n1(end) - g2*n2(1))/da - mu_im*n2(1); 
    % Nós internos da malha do Estágio 2
    dydt(p.idx_L(2:end)) = -(g2/da)*(n2(2:end) - n2(1:end-1)) - mu_im*n2(2:end);
    
    % --- EQUACIONAMENTO DO ESTÁGIO 3 (Adultos) ---
    % Nó de Fronteira: O fluxo de entrada é a emergência das pupas do estágio 2
    dydt(p.idx_A(1)) = (g2*n2(end) - g3*n3(1))/da - mu_ad_total_vec(1)*n3(1); 
    % Nós internos da malha do Estágio 3 (Nota: A mortalidade aqui é matrizada pontualmente '.*')
    dydt(p.idx_A(2:end)) = -(g3/da)*(n3(2:end) - n3(1:end-1)) - (mu_ad_total_vec(2:end) .* n3(2:end));
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
        val(mask) = params.a .* T_m .* (T_m - params.Tmin) .* sqrt(abs(params.Tmax - T_m));
    end
end

% Função quadrática de mortalidade térmica para fases jovens
function mu = local_calc_mu_imature(T, p)
    mu = max(0, polyval(p.mu_im_coeffs, T)); 
end

% Função quadrática de mortalidade térmica específica para fase reprodutiva
function mu_therm = local_calc_mu_adult_thermal(T, p)
    mu_therm = max(0, polyval(p.mu_ad_coeffs, T)); 
end

% A senescência cria um "paredão" de mortalidade ao fim da vida fisiológica, 
% convergindo artificialmente a densidade para zero no limiar do A3.
function mu_age = local_calc_senescence(age, p)
    mu_age = (p.mu_ad_senescencia .* age) ./ max(1e-6, (p.A3 - age).^2); 
end

% -------------------------------------------------------------------------
% GERADOR DE CURVAS PARA O VETOR DE CONDIÇÃO INICIAL (t=0)
% -------------------------------------------------------------------------
function n0 = gaussian_dist(ages, A_max, Pop_Ini, da)
    sigma = A_max / 3; % Desvio base para a distribuição suave na malha
    
    % Distribui os indivíduos pelo espaço fisiológico
    n0 = exp(-((ages - ages(1)).^2)/(2*sigma^2));
    
    % Ajuste fino (Shift): Evita descontinuidades abruptas puxando a cauda para exatamente zero
    n0 = n0 - n0(end); 
    n0(n0 < 0) = 0; % Salvaguarda numérica
    
    % Fator K de normalização da integral numérica (Garante que a área total seja = Pop_Ini)
    n0 = n0 * (Pop_Ini / (sum(n0)*da)); 
end

% -------------------------------------------------------------------------
% UTILITÁRIO GRÁFICO GLOBAL
% -------------------------------------------------------------------------
function set_ax(xlbl, ylbl, x_lim)
    xlabel(xlbl); 
    ylabel(ylbl); 
    grid on;
    if nargin == 3 && ~isempty(x_lim)
        xlim(x_lim); 
    end
end