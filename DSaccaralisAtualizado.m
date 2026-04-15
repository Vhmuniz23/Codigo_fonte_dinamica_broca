% =========================================================================
% MODELO MCKENDRICK-VON FOERSTER DE 3 ESTÁGIOS (OVO-LARVA-ADULTO)               
% Espécie: Diatraea saccharalis (Broca-da-cana)                                 
% Referência: Baseado em parâmetros de Carbognin et al. (2023)                  
%                                                                               
% DESCRIÇÃO MATEMÁTICA:                                                         
% Resolve a Equação Diferencial Parcial (EDP) de transporte com estrutura       
% etária usando o "Método das Linhas" (Method of Lines).                        
% dn/dt + d(g*n)/da = -mu*n                                                     
% Onde:                                                                         
% n(a,t): Densidade populacional na idade 'a' e tempo 't'                     
% g(T,N): Taxa de desenvolvimento (depende de Temp e Densidade)               
% mu(T) : Taxa de mortalidade (depende da Temp)                               
% =========================================================================
clear all;  
clc;        
close all;  

% =========================================================================
% 1. CONFIGURAÇÃO DO AMBIENTE E FORÇANTE CLIMÁTICA
% =========================================================================
p.T_media_C = 25.0;       % Temperatura média base (°C)
p.amplitude_C = 3.0;      % Amplitude da variação térmica (ex: 20°C a 30°C)
p.periodo = 365;          % Período da oscilação (1 ano)

% Função anônima para temperatura diária (Forçante senoidal)
% Simula a sazonalidade ao longo do ano
p.T_func = @(t) p.T_media_C + p.amplitude_C * sin(2 * pi * t / p.periodo);

% =========================================================================
% 2. PARÂMETROS BIOLÓGICOS
% =========================================================================

% 2.1. Desenvolvimento (Modelo Brière-1)
% O modelo Brière captura a não-linearidade: crescimento lento em baixas temps,
% pico no ótimo, e queda abrupta no limite letal superior.
p.briere_egg_a = 1.85e-4;     % Constante de inclinação (Ovos)
p.briere_egg_Tmin = 12.8;     % Limiar térmico inferior (Ovos)
p.briere_egg_Tmax = 36.8;     % Limiar térmico superior (Ovos)

p.briere_larva_a = 3.60e-5;   % Constante de inclinação (Larvas)
p.briere_larva_Tmin = 14.3;   % Limiar térmico inferior (Larvas)
p.briere_larva_Tmax = 36.3;   % Limiar térmico superior (Larvas)

% 2.2. Fecundidade (Modelo Gaussiano)
% Define quantos ovos uma fêmea põe por dia dependendo da temperatura atual.
p.b_max = 93.63;   % Taxa máxima de oviposição (ovos/dia no pico)
p.T_opt_b = 24.27; % Temperatura ótima para postura
p.sigma_b = 4.46;  % Largura da curva (tolerância térmica para reprodução)

% 2.3. Mortalidade (Polinomial: aT^2 + bT + c)
% A mortalidade aumenta nos extremos de temperatura (Curva em "U").
% Imaturos (Ovo + Larva)
p.mu_imm_a = 0.0002;
p.mu_imm_b = -0.0104;
p.mu_imm_c = 0.1362;
p.mu_lethal_T = 36.0;  % Temperatura letal crítica
p.mu_val_lethal = 1.0; % Mortalidade massiva acima da temp crítica (100%/dia)

% Adultos
p.mu_ad_a = 0.011;
p.mu_ad_b = -0.572;
p.mu_ad_c = 7.586;

% 2.4. Competição Intraespecífica (Densidade-Dependência)
% Parâmetro crucial para estabilidade do modelo.
% Afeta o 'g' (desenvolvimento): Com muitos indivíduos, o crescimento desacelera.
p.q_crowding = 5e-5; % Coeficiente de sensibilidade ao crowding

% 2.5. Configuração da Grade de Idade (Discretização)
T_ref = 25.0; % Temperatura de referência para definir o "tamanho" da grade

% Função auxiliar local para calcular Brière
calc_briere = @(T, a, Tmin, Tmax) max(0, (T>Tmin & T<Tmax) .* (a*T*(T-Tmin)*sqrt(abs(Tmax-T))));

% Calcula taxas na temperatura de referência para normalizar o tempo
rate_egg_ref = calc_briere(T_ref, p.briere_egg_a, p.briere_egg_Tmin, p.briere_egg_Tmax);
rate_larva_ref = calc_briere(T_ref, p.briere_larva_a, p.briere_larva_Tmin, p.briere_larva_Tmax);

% Armazena taxas de referência (usadas para calcular o fator 'g' na ODE)
p.rate_egg_ref = rate_egg_ref;
p.rate_larva_ref = rate_larva_ref;

% Durações médias (em dias) na temperatura de referência
dur_ovo_ref = 1 / rate_egg_ref;
dur_larva_ref = 1 / rate_larva_ref;
dur_adulto_ref = 10.0; % Duração fixa arbitrária para adultos
% dur_ovo_ref = 5.0;
% dur_larva_ref = 39.0
% dur_adulto_ref = 6.0;

% Define os pontos de corte na grade de idade (Idade Cronológica Fixa)
idade_eclosao = dur_ovo_ref;
idade_emergencia = dur_ovo_ref + dur_larva_ref;
idade_max = idade_emergencia + dur_adulto_ref;

% Criação do vetor de idades (eixo 'a' da EDP)
n_classes_total = 600;                             % Resolução espacial (número de compartimentos)
idades = linspace(0, idade_max, n_classes_total)'; % Vetor coluna de idades
da = idades(2) - idades(1);                        % Passo da discretização (delta a)
p.da = da;

% Identificação dos índices que separam os estágios no vetor
idx_eclosao = find(idades >= idade_eclosao, 1);
idx_emergencia = find(idades >= idade_emergencia, 1);

% Separação dos vetores de idade para plotagem
idades_O = idades(1:idx_eclosao-1);
idades_L = idades(idx_eclosao:idx_emergencia-1);
idades_A = idades(idx_emergencia:end);

% Armazena o tamanho (número de células) de cada estágio
p.n_O = length(idades_O);
p.n_L = length(idades_L);
p.n_A = length(idades_A);

% =========================================================================
% 3. CONDIÇÕES INICIAIS
% =========================================================================
pop_ini_O = 5000; pop_ini_L = 2000; pop_ini_A = 500; 

% Distribuição Gaussiana inicial: Evita choques numéricos de um pulso quadrado.
% Distribui a população inicial ao redor da média de idade de cada estágio.
n0_O = pop_ini_O * exp(-((idades_O - mean(idades_O)).^2)/(2*2^2));
n0_L = pop_ini_L * exp(-((idades_L - mean(idades_L)).^2)/(2*10^2));
n0_A = pop_ini_A * exp(-((idades_A - mean(idades_A)).^2)/(2*2^2));

% Monta o vetor de estado único (y0) para o solver
y0 = [n0_O; n0_L; n0_A]; 

% =========================================================================
% 4. SOLVER ODE45 (INTEGRAÇÃO TEMPORAL)
% =========================================================================
t_span = [0 3650]; % Simulação de 10 anos (3650 dias)
% t_span = [0 5475]; % Simulação de 15 anos (5475 dias)

disp('Iniciando simulação (Competição no Desenvolvimento)...');
% Tolerâncias ajustadas para garantir precisão na conservação de massa
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4); 

tic;
% Chama o integrador numérico. A mágica acontece em 'odefun_3stage'
[t, y_sol] = ode45(@(t,y) odefun_3stage(t, y, p), t_span, y0, opts);
tempo_exec = toc;
fprintf('Simulação concluída em %.2f segundos.\n', tempo_exec);

% =========================================================================
% 5. PÓS-PROCESSAMENTO E PLOTS
% =========================================================================

% Recupera índices para fatiar a matriz solução 'y_sol'
idx1 = p.n_O; 
idx2 = p.n_O + p.n_L;

% Fatiamento: y_sol tem dimensões [tempo x (nO+nL+nA)]
n_sol_O = y_sol(:, 1:idx1);
n_sol_L = y_sol(:, idx1+1:idx2);
n_sol_A = y_sol(:, idx2+1:end);

% Integração numérica (Soma de Riemann) para obter População Total vs Tempo
total_O = sum(n_sol_O, 2) * da;
total_L = sum(n_sol_L, 2) * da;
total_A = sum(n_sol_A, 2) * da;

% Definição de Cores
green = [0, 1, 0.498];          
yellow = [1, 0.871, 0.0129];    
red = [1, 0.173, 0.173];        
colors = [green; yellow; red];  
color_taxa = [0 0.4470 0.7410]; 
color_temp = [0.8500 0.3250 0.0980]; 

% 5.1 Plot das Condições Iniciais
figure('Name', 'Condições Iniciais – Ovos'); 
plot(idades_O, n0_O, 'Color', colors(1,:), 'LineWidth', 2); 
title('Distribuição Inicial: Ovos'); xlabel('Idade (dias)'); ylabel('Densidade'); grid on; drawnow;

figure('Name', 'Condições Iniciais – Larvas'); 
plot(idades_L, n0_L, 'Color', colors(2,:), 'LineWidth', 2); 
title('Distribuição Inicial: Larvas'); xlabel('Idade (dias)'); ylabel('Densidade'); grid on; drawnow;
 
figure('Name', 'Condições Iniciais – Adultos'); 
plot(idades_A, n0_A, 'Color', colors(3,:), 'LineWidth', 2); 
title('Distribuição Inicial: Adultos'); xlabel('Idade (dias)'); ylabel('Densidade'); grid on; drawnow;

% 5.2 Gráficos 3D de Evolução (Superfície)
% Permite visualizar ondas populacionais viajando no tempo e idade
figure('Name', 'Evolução 3D – Ovos'); 
surf(idades_O, t, n_sol_O, 'EdgeColor', 'none'); 
shading interp; xlabel('Idade (dias)'); ylabel('Tempo (dias)'); zlabel('Densidade');
title('Dinâmica Espaço-Temporal: Ovos'); colormap('jet'); view(3); grid on; drawnow;

figure('Name', 'Evolução 3D – Larvas'); 
surf(idades_L, t, n_sol_L, 'EdgeColor', 'none');
shading interp; xlabel('Idade (dias)'); ylabel('Tempo (dias)'); zlabel('Densidade');
title('Dinâmica Espaço-Temporal: Larvas'); colormap('jet'); view(3); grid on; drawnow;

figure('Name', 'Evolução 3D – Adultos'); 
surf(idades_A, t, n_sol_A, 'EdgeColor', 'none');
shading interp; xlabel('Idade (dias)'); ylabel('Tempo (dias)'); zlabel('Densidade');
title('Dinâmica Espaço-Temporal: Adultos'); colormap('jet'); view(3); grid on; drawnow;

% 5.3 Dinâmica Total por Estágio
figure('Name', 'Dinâmica Populacional Total por Estágio');
subplot(3,1,1); 
plot(t, total_O, 'Color', colors(1,:), 'LineWidth', 2); 
title('População Total: Ovos'); ylabel('Indivíduos'); grid on;

subplot(3,1,2); 
plot(t, total_L, 'Color', colors(2,:), 'LineWidth', 2); 
title('População Total: Larvas'); ylabel('Indivíduos'); grid on;

subplot(3,1,3); 
plot(t, total_A, 'Color', colors(3,:), 'LineWidth', 2); 
title('População Total: Adultos'); ylabel('Indivíduos'); xlabel('Tempo (dias)'); grid on;
drawnow;

% 5.4 Temperatura Sazonal
temp_t_C = arrayfun(p.T_func, t); 
figure('Name', 'Variação Sazonal da Temperatura');
plot(t, temp_t_C, 'LineWidth', 2, 'Color', color_temp); 
xlabel('Tempo (dias)'); ylabel('Temperatura (°C)'); title('Temperatura Sazonal'); grid on;
drawnow;

% =========================================================================
% CÁLCULO VETORIZADO PARA ANÁLISE PÓS-SIMULAÇÃO
% =========================================================================
% Recalculamos as taxas vitais fora do loop ODE para poder plotá-las
disp('Recalculando taxas vitais e fluxos para plotagem...');
briere_vec = @(T, a, Tmin, Tmax) max(0, (T>Tmin & T<Tmax) .* (a.*T.*(T-Tmin).*sqrt(abs(Tmax-T))));

% 1. Desenvolvimento e Competição
r_O_vec = briere_vec(temp_t_C, p.briere_egg_a, p.briere_egg_Tmin, p.briere_egg_Tmax);
r_L_vec_pure = briere_vec(temp_t_C, p.briere_larva_a, p.briere_larva_Tmin, p.briere_larva_Tmax);

v_taxa_g_O = r_O_vec / p.rate_egg_ref;

% AQUI ESTÁ A CHAVE DA ANÁLISE:
% v_taxa_g_L_pure: Desenvolvimento teórico baseado apenas na temperatura
% v_taxa_g_L_eff:  Desenvolvimento real após penalidade por superpopulação
v_taxa_g_L_pure = r_L_vec_pure / p.rate_larva_ref;
v_taxa_g_L_eff = v_taxa_g_L_pure ./ (1 + p.q_crowding * total_L); 

% 2. Mortalidade
mu_im_vec = p.mu_imm_a * temp_t_C.^2 + p.mu_imm_b * temp_t_C + p.mu_imm_c;
mu_im_vec(temp_t_C > p.mu_lethal_T) = p.mu_val_lethal;

mu_ad_vec = p.mu_ad_a * temp_t_C.^2 + p.mu_ad_b * temp_t_C + p.mu_ad_c;
mu_ad_vec(temp_t_C > p.mu_lethal_T) = p.mu_val_lethal;

% 3. Fecundidade
fec = p.b_max * exp(-0.5 * ((temp_t_C - p.T_opt_b)/p.sigma_b).^2);
v_taxa_b = fec; 

% 4. Fluxos entre estágios (Boundary Conditions Outputs)
% Fluxo = Velocidade de desenvolvimento * Densidade na última idade
v_fluxo_eclosao = v_taxa_g_O .* n_sol_O(:, end); 
v_fluxo_emergencia = v_taxa_g_L_eff .* n_sol_L(:, end);

% 5.5 Plots de Taxas Vitais e Efeito da Competição
figure('Name', 'Taxas Vitais e Efeito da Competição');
subplot(2,1,1); 
yyaxis left; plot(t, v_taxa_b, 'LineWidth', 2, 'Color', color_taxa); ylabel('Ovos/Fêmea'); title('Fecundidade Potencial (Influência da Temperatura)');
yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp); ylabel('Temperatura (°C)'); ax=gca; ax.YColor = color_temp; grid on;

subplot(2,1,2); 
% Visualização do atraso no desenvolvimento causado pela densidade
plot(t, v_taxa_g_L_pure, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
plot(t, v_taxa_g_L_eff, 'LineWidth', 2, 'Color', color_taxa); 
legend('g Potencial (Só Temperatura)', 'g Real (Temperatura + Competição)');
title('Aceleração do Desenvolvimento Larval (g)'); ylabel('Fator g'); xlabel('Tempo (dias)'); grid on;
drawnow;

% =========================================================================
% 5.5b GRÁFICOS DE MORTALIDADE COM SOBREPOSIÇÃO DE TEMPERATURA
% =========================================================================
% Visualiza como o risco de morte varia ao longo do ano devido à temperatura.
% O modelo usa um polinômio quadrático (Parábola), o que significa que a 
% mortalidade aumenta tanto no frio extremo quanto no calor extremo.

figure('Name', 'Mortalidade (Modelo Quadrático) e Temperatura');

% Subplot 1: Imaturos (Ovos e Larvas) 
subplot(2,1,1); 
yyaxis left; % Eixo para a Taxa de Mortalidade
plot(t, mu_im_vec, 'Color', colors(2,:), 'LineWidth', 2); 
ylabel('Taxa Instantânea (dia^{-1})');
title('Mortalidade Natural: Imaturos (Ovos/Larvas)');
grid on;

yyaxis right; % Eixo para a Temperatura
plot(t, temp_t_C, '--', 'Color', color_temp);
ylabel('Temperatura (°C)');
ax = gca; ax.YColor = color_temp; % Ajusta a cor do eixo para combinar com a linha

% Subplot 2: Adultos 
subplot(2,1,2); 
yyaxis left; % Eixo para a Taxa de Mortalidade
plot(t, mu_ad_vec, 'Color', colors(3,:), 'LineWidth', 2); 
ylabel('Taxa Instantânea (dia^{-1})');
title('Mortalidade Natural: Adultos');
grid on;

yyaxis right; % Eixo para a Temperatura
plot(t, temp_t_C, '--', 'Color', color_temp);
ylabel('Temperatura (°C)');
xlabel('Tempo (dias)');
ax = gca; ax.YColor = color_temp;

% Comentário de interpretação:
% Note que nos vales e picos da linha tracejada (temperatura), 
% a linha sólida (mortalidade) sobe, confirmando o comportamento 
% parabólico do modelo mu(T).
drawnow;

% 5.6 Fluxos 
figure('Name', 'Fluxos Populacionais');
subplot(2,1,1);
yyaxis left; plot(t, v_fluxo_eclosao, 'LineWidth', 2, 'Color', color_taxa); ylabel('Indivíduos/dia'); title('Eclosão (Ovo -> Larva)');
yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp); ylabel('Temperatura (°C)'); ax=gca; ax.YColor = color_temp; grid on;

subplot(2,1,2);
yyaxis left; plot(t, v_fluxo_emergencia, 'LineWidth', 2, 'Color', color_taxa); ylabel('Indivíduos/dia'); title('Maturação (Larva -> Adulto)'); xlabel('Tempo (dias)');
yyaxis right; plot(t, temp_t_C, '--', 'Color', color_temp); ylabel('Temperatura (°C)'); ax=gca; ax.YColor = color_temp; grid on;
drawnow;

disp('Plots gerados com sucesso.');

% =========================================================================
% FUNÇÃO PRINCIPAL DA EDO (CÁLCULO DAS DERIVADAS)
% Esta função é chamada a cada passo de tempo pelo ode45
% =========================================================================
function dydt = odefun_3stage(t, y, p)
    % 1. Atualiza temperatura para o tempo atual 't'
    T = p.T_func(t); 
    
    % 2. Recupera índices para separar o vetorzão 'y'
    idx1 = p.n_O; 
    idx2 = p.n_O + p.n_L;
    
    % 3. Desempacota o estado atual
    n_O = y(1:idx1);      
    n_L = y(idx1+1:idx2); 
    n_A = y(idx2+1:end);  
    
    % 4. Calcula Populações Totais (Integração no espaço 'a')
    % Necessário para a densidade-dependência e nascimentos
    N_L_tot = sum(n_L)*p.da; 
    N_A_tot = sum(n_A)*p.da; 
    
    % CÁLCULO DAS TAXAS VITAIS INSTANTÂNEAS
    
    % Função local Brière
    briere = @(temp, a, Tmin, Tmax) max(0, (temp>Tmin & temp<Tmax) .* (a*temp*(temp-Tmin)*sqrt(abs(Tmax-temp))));
    
    % Taxas brutas de desenvolvimento (1/dias)
    r_egg = briere(T, p.briere_egg_a, p.briere_egg_Tmin, p.briere_egg_Tmax);
    r_larva = briere(T, p.briere_larva_a, p.briere_larva_Tmin, p.briere_larva_Tmax);
    
    % Fator de Aceleração 'g' (Adimensional ou normalizado)
    % Se g=1, o inseto envelhece na mesma velocidade do tempo cronológico.
    % Se g=2, ele envelhece o dobro (temperatura quente).
    g_O = r_egg / p.rate_egg_ref; 
    g_L_temp = r_larva / p.rate_larva_ref; 
    
    % APLICAÇÃO DA COMPETIÇÃO NO DESENVOLVIMENTO
    % Modelo: g_efetivo = g_temperatura / (1 + q * População)
    % Efeito: Aumenta o tempo que a larva fica no estágio (atraso fenológico)
    g_L = g_L_temp / (1 + p.q_crowding * N_L_tot);
    
    g_A = 1.0; % Adultos envelhecem cronologicamente (sem dependência de T neste modelo)
    
    % Mortalidade (Polinomial) 
    mu_imm = p.mu_imm_a * T^2 + p.mu_imm_b * T + p.mu_imm_c;
    mu_imm = max(0, mu_imm); if T > p.mu_lethal_T, mu_imm = p.mu_val_lethal; end
    
    mu_ad = p.mu_ad_a * T^2 + p.mu_ad_b * T + p.mu_ad_c;
    mu_ad = max(0, mu_ad); if T > p.mu_lethal_T, mu_ad = p.mu_val_lethal; end

    % Cria vetores de mortalidade (pode ser expandido para idade-dependente no futuro)
    mu_vec_O = mu_imm * ones(size(n_O));
    mu_vec_L = mu_imm * ones(size(n_L)); 
    mu_vec_A = mu_ad * ones(size(n_A));
    
    % Reprodução (Condição de Contorno em a=0)
    fec = p.b_max * exp(-0.5 * ((T - p.T_opt_b)/p.sigma_b)^2);
    birth_rate = fec * N_A_tot; % Total de ovos produzidos pela população adulta
    
    % =====================================================================
    % DISCRETIZAÇÃO DA EDP (MÉTODO UPWIND / DIFERENÇAS FINITAS) 
    % EDP: dn/dt = -g * dn/da - mu * n
    % Discretizado: dn[i]/dt = - (g/da)*(n[i] - n[i-1]) - mu[i]*n[i]
    % =====================================================================
    
    da = p.da;
    dydt = zeros(size(y));
    
    % ESTÁGIO 1: OVOS 
    % Célula 1 (Condição de Contorno Esquerda): Entrada de nascimentos
    dydt(1) = -(g_O/da)*n_O(1) + (g_O/da)*birth_rate - mu_vec_O(1)*n_O(1);
    
    % Células 2 até fim (Transporte interno)
    dydt(2:idx1) = -(g_O/da)*(n_O(2:end) - n_O(1:end-1)) - mu_vec_O(2:end).*n_O(2:end);
    
    % Fluxo de saída (Ovos -> Larvas): O que sai da última célula de ovo
    fluxo_O_L = (g_O/da) * n_O(end) * da; 
    
    % ESTÁGIO 2: LARVAS
    % Célula 1: Recebe o fluxo de eclosão dos ovos
    dydt(idx1+1) = -(g_L/da)*n_L(1) + (fluxo_O_L/da) - mu_vec_L(1)*n_L(1);
    
    % Células internas (Transporte com g_L afetado por densidade)
    dydt(idx1+2:idx2) = -(g_L/da)*(n_L(2:end) - n_L(1:end-1)) - mu_vec_L(2:end).*n_L(2:end);
    
    % Fluxo de saída (Larvas -> Adultos)
    fluxo_L_A = (g_L/da) * n_L(end) * da; 
    
    % ESTÁGIO 3: ADULTOS
    % Célula 1: Recebe o fluxo de emergência das larvas
    dydt(idx2+1) = -(g_A/da)*n_A(1) + (fluxo_L_A/da) - mu_vec_A(1)*n_A(1);
    
    % Células internas
    dydt(idx2+2:end) = -(g_A/da)*(n_A(2:end) - n_A(1:end-1)) - mu_vec_A(2:end).*n_A(2:end);
end