# Modelo de Dinâmica Populacional da Broca-da-cana

Este repositório contém os códigos em MATLAB para a simulação de um modelo matemático que descreve a dinâmica populacional da broca-da-cana (*Diatraea saccharalis*). 

O modelo é estruturado em idade fisiológica e utiliza Equações Diferenciais Parciais (EDPs) para capturar os "trade-offs" ecológicos e janelas de surtos populacionais da praga, fornecendo subsídios quantitativos para estratégias de Manejo Integrado de Pragas (MIP).

# Principais Funcionalidades do Modelo

* **Estruturação por Estágios:** Divide a população em três estágios fisiológicos (ovo, larva/pupa e adulto).
* **Dependência Térmica:** Simula o desenvolvimento da praga sob regimes de temperatura constante e sob forçantes climáticas sazonais (variáveis no tempo).
* **Regulação Densidade-Dependente:** Incorpora uma função de recrutamento do tipo Beverton-Holt para limitar o crescimento populacional.
* **Cálculo do Número Reprodutivo Básico (R_0):** Avalia os limites de tolerância letal e a viabilidade da praga no ecossistema.
* **Método das Características (MC):** Resolve analítica e numericamente a conversão das EDPs em equações integrais ao longo do tempo.

# Tecnologias Utilizadas

* MATLAB

# Aplicação

Este código foi desenvolvido como parte das experimentações computacionais para uma pesquisa de Pós-Graduação, com o objetivo de otimizar o controle biológico em culturas de cana-de-açúcar.

# Ajuste de Parâmetros Termobiológicos: *Diatraea saccharalis*

Esta parte do repositório contém os scripts computacionais em MATLAB desenvolvidos para calibrar as funções contínuas de taxa de desenvolvimento e mortalidade da broca-da-cana, *Diatraea saccharalis* (Lepidoptera: Crambidae). 

Os parâmetros extraídos nestes ajustes são fundamentais para alimentar simulações numéricas de dinâmica populacional (modelos estruturados por idade) e delimitar o nicho térmico fundamental da espécie.

## 📌 Base de Dados

Os ajustes empíricos foram realizados tomando como base os dados experimentais de tempo de desenvolvimento ($\tau$) e sobrevivência reportados na literatura especializada, com destaque primário para:
* **Carbognin et al. (2023)**: *Unraveling the effect of temperature and humidity on the life cycle of Diatraea saccharalis...* Environmental Entomology.

## 🧮 Modelos Matemáticos Implementados

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