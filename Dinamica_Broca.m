% =========================================================================
% MODELO POPULACIONAL ESTRUTURADO EM IDADE (EDP) - DIATRAEA SACCHARALIS
% =========================================================================
% TIPO DE MODELO: 
% Sistema de Equacoes Diferenciais Parciais (EDP) de Transporte 
% (Baseado na equacao de McKendrick-von Foerster).
% 
% EQUACAO GOVERNANTE GERAL (Para cada estagio):
% dn/dt + d(g(T)*n)/da = -mu(T,a)*n
%
% ONDE:
%    n(t,a): Densidade de individuos (indivíduos por dia de idade)
%         t: Tempo Cronologico (dias)
%         a: Idade fisiologica (adimensional ou normalizada)
%      g(T): Taxa de adveccao (velocidade de envelhecimento fisiologico 
%            dependente da temperatura T)
%   mu(T,a): Taxa de mortalidade instantanea (depende da temperatura T e 
%            da idade a, representando estresse termico e senescencia)
%
% METODO DE RESOLUCAO NUMERICA:
% Metodo das Linhas (Method of Lines - MOL). A dimensao espacial da idade 
% ‘a’ e discretizada em compartimentos (caixas numericas) usando diferencas 
% finitas (esquema Upwind). Isso transforma a EDP num grande sistema acoplado 
% de Equacoes Diferenciais Ordinarias (EDOs), que e resolvido ao longo do 
% tempo ’t’ utilizando o integrador ‘ode15s’ (adequado para sistemas rigidos).
% =========================================================================

% [MAT] Limpeza completa do ambiente para evitar conflitos de memoria
clearvars; 
clc; 
close all;

% — PADRONIZACAO GRAFICA GLOBAL —
set(groot, ‘defaultAxesFontSize’, 12);                   % Eixos e numeros (Tick labels) tamanho 12
set(groot, ‘defaultTextFontSize’, 12);                   % Textos gerais e legendas tamanho 12
set(groot, ‘defaultAxesTitleFontSizeMultiplier’, 14/12); % Forca os titulos para tamanho 14
