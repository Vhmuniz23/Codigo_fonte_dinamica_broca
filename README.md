# Modelagem Populacional da *Diatraea saccharalis*

Este repositório contém uma série de scripts desenvolvidos em **MATLAB** para a modelagem matemática e simulação da dinâmica populacional da broca-da-cana (*Diatraea saccharalis*). 

O modelo é estruturado em idade fisiológica e utiliza Equações Diferenciais Parciais (EDPs) de transporte e métodos numéricos para análise de sensibilidade térmica e sazonalidade. Além disso, captura os "trade-offs" ecológicos e janelas de surtos populacionais da praga, fornecendo subsídios quantitativos para estratégias de Manejo Integrado de Pragas (MIP).

# Principais Funcionalidades do Modelo

* **Estruturação por Estágios:** Divide a população em três estágios fisiológicos (ovo, larva/pupa e adulto).
* **Dependência Térmica:** Simula o desenvolvimento da praga sob regimes de temperatura constante e sob forçantes climáticas sazonais (variáveis no tempo).
* **Regulação Densidade-Dependente:** Incorpora uma função de recrutamento do tipo Beverton-Holt para limitar o crescimento populacional.
* **Cálculo do Número Reprodutivo Básico (R_0):** Avalia os limites de tolerância letal e a viabilidade da praga no ecossistema.
* **Método das Características (MC):** Resolve analítica e numericamente a conversão das EDPs em equações integrais ao longo do tempo.

## Estrutura do Repositório

O projeto é composto por três módulos principais:

1. **`Ajuste_Curva_Desenvolvimento.m`**: Realiza o ajuste não-linear dos dados experimentais ao modelo de Brière-1, determinando parâmetros térmicos críticos ($T_{min}$, $T_{max}$, $a$) para o desenvolvimento dos diferentes estágios da broca.
2. **`Ajuste_Curva_Mortalidade.m`**: Processa as taxas de mortalidade através de regressão polinomial e estima o componente de senescência para a fase adulta, correlacionando temperatura e idade fisiológica.
3. **`Dinamica_Broca_final.m`**: O núcleo do projeto. Implementa o sistema de EDPs de McKendrick-von Foerster via Método das Linhas (MOL). Este script simula a dinâmica populacional ao longo do tempo, considerando forçantes climáticas (sazonalidade) e competição intraespecífica.

# Ajustes de Parâmetros Termobiológicos: *Diatraea saccharalis*

Os scripts foram desenvolvidos para calibrar as funções contínuas de taxa de desenvolvimento e mortalidade da broca-da-cana, *Diatraea saccharalis*. 

Os parâmetros extraídos nestes ajustes são fundamentais para alimentar simulações numéricas de dinâmica populacional (modelos estruturados por idade) e delimitar o nicho térmico fundamental da espécie.

## 🧮 Modelos Matemáticos Implementados nos Ajustes de Parâmetros 

### 1. Taxa de Desenvolvimento (Estágios Imaturos)
A taxa de desenvolvimento absoluta, definida como $R(T) = 1/\tau$, foi ajustada utilizando o modelo termobiológico não-linear de **Brière-1**. Este modelo captura a assimetria biológica do desenvolvimento de ectotérmicos:

$$R_i(T) = \gamma_i T (T - T_{min,i}) \sqrt{T_{max,i} - T}$$

Para temperaturas fora do intervalo biologicamente viável ($T \leq T_{min,i}$ ou $T \geq T_{max,i}$), assume-se $R_i(T) = 0$, para i \in \{1,2\}.

### 2. Taxa de Mortalidade (Estágios Imaturos)
A mortalidade dependente da temperatura para os Estágios 1 (ovos a L3) e 2 (L4 a pupa) foi descrita por uma equação polinomial de segundo grau (parábola com concavidade voltada para cima, refletindo o estresse térmico em extremos):

$$d_i(T) = a_{i} T^2 + b_{i} T + c_{i}, i \in \{1,2\}.$$

### 3. Taxa de Mortalidade com Senescência (Adultos)
A mortalidade na fase adulta (Estágio 3) requer a integração do fator térmico com o desgaste fisiológico natural (idade, $a$). A função implementada garante que a mortalidade tenda ao infinito quando o inseto se aproxima da sua longevidade máxima assumida ($A_3$):

$$\mu_3(a, T) = (a_{3} T^2 + b_{3} T + c_{3}) \cdot d_{3} \frac{a}{(A_{3} - a)^2},$$

onde $d_3$ é o parâmetro escalar de senescência ajustado iterativamente via rotinas de otimização.

## 📌 Base de Dados

Os ajustes empíricos foram realizados tomando como base os dados experimentais de tempo de desenvolvimento ($\tau$) e sobrevivência reportados na literatura especializada, com destaque primário para:
* **Carbognin et al. (2023)**: *Unraveling the effect of temperature and humidity on the life cycle of Diatraea saccharalis...* Environmental Entomology.

## Pré-requisitos

* **MATLAB** (Recomendado versão R2020a ou superior).
* **Optimization Toolbox** (Necessário para as funções `lsqcurvefit`).
* **Parallel Computing Toolbox** (Opcional, utilizado para acelerar as simulações em lote/varredura).

## Como Utilizar

1.  Clone este repositório ou baixe os arquivos `.m`.
2.  Abra o MATLAB e defina o diretório atual para a pasta onde os arquivos foram salvos.
3.  **Para análise de dados**: Execute `Ajuste_Curva_Desenvolvimento.m` e `Ajuste_Curva_Mortalidade.m` para visualizar os ajustes estatísticos e obter os coeficientes necessários.
4.  **Para simulação populacional**: Execute `Dinamica_Broca_final.m`. O script irá:
    * Configurar a malha de idade fisiológica.
    * Executar simulações em lote para diferentes cenários térmicos.
    * Gerar gráficos de abundância populacional (mínimos e máximos) em função da temperatura média anual e amplitude sazonal.

## Metodologia

O modelo utiliza o **Método das Linhas (Method of Lines - MOL)** para discretizar a idade em compartimentos (esquema *upwind*) e o integrador `ode15s` para resolver a evolução temporal do sistema rígido. As curvas de desenvolvimento seguem o modelo de Brière-1, e a mortalidade é tratada como uma função dependente de temperatura e senescência etária.

## Contribuições

Sinta-se à vontade para abrir *Issues* ou enviar *Pull Requests* caso identifique melhorias nos modelos ou na implementação numérica.

---
*Desenvolvido como parte das experimentações computacionais para uma pesquisa de Pós-Graduação, com o objetivo de monitorar a dinâmica populacional da praga e otimizar o controle biológico em culturas de cana-de-açúcar.*