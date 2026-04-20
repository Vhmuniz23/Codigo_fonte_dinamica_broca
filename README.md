% Modelo de Dinâmica Populacional da Broca-da-cana

Este repositório contém os códigos em MATLAB para a simulação de um modelo matemático que descreve a dinâmica populacional da broca-da-cana (*Diatraea saccharalis*). 

O modelo é estruturado em idade fisiológica e utiliza Equações Diferenciais Parciais (EDPs) para capturar os "trade-offs" ecológicos e janelas de surtos populacionais da praga, fornecendo subsídios quantitativos para estratégias de Manejo Integrado de Pragas (MIP).

% Principais Funcionalidades do Modelo

* **Estruturação por Estágios:** Divide a população em três estágios fisiológicos (ovo, larva/pupa e adulto).
* **Dependência Térmica:** Simula o desenvolvimento da praga sob regimes de temperatura constante e sob forçantes climáticas sazonais (variáveis no tempo).
* **Regulação Densidade-Dependente:** Incorpora uma função de recrutamento do tipo Beverton-Holt para limitar o crescimento populacional.
* **Cálculo do Número Reprodutivo Básico (R_0):** Avalia os limites de tolerância letal e a viabilidade da praga no ecossistema.
* **Método das Características (MC):** Resolve analítica e numericamente a conversão das EDPs em equações integrais ao longo do tempo.

% Tecnologias Utilizadas
* MATLAB

% Aplicação
Este código foi desenvolvido como parte das experimentações computacionais para uma pesquisa de Pós-Graduação, com o objetivo de otimizar o controle biológico em culturas de cana-de-açúcar.
