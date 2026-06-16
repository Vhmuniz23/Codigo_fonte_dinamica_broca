# Ajuste de Parâmetros Termobiológicos: *Diatraea saccharalis*

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Este repositório contém os scripts computacionais em MATLAB desenvolvidos para calibrar as funções contínuas de taxa de desenvolvimento e mortalidade da broca-da-cana, *Diatraea saccharalis* (Lepidoptera: Crambidae). 

Os parâmetros extraídos nestes ajustes são fundamentais para alimentar simulações numéricas de dinâmica populacional (modelos estruturados por idade) e delimitar o nicho térmico fundamental da espécie.

## 📌 Base de Dados

Os ajustes empíricos foram realizados tomando como base os dados experimentais de tempo de desenvolvimento ($\tau$) e sobrevivência reportados na literatura especializada, com destaque primário para:
* **Carbognin et al. (2023)**: *Unraveling the effect of temperature and humidity on the life cycle of Diatraea saccharalis...* Environmental Entomology.

## 🧮 Modelos Matemáticos Implementados

### 1. Taxa de Desenvolvimento (Estágios Imaturos)
A taxa de desenvolvimento absoluta, definida como $R(T) = 1/\tau$, foi ajustada utilizando o modelo termobiológico não-linear de **Brière-1**. Este modelo captura a assimetria biológica do desenvolvimento de ectotérmicos:

$$R(T) = \gamma T (T - T_{min}) \sqrt{T_{max} - T}$$

Para temperaturas fora do intervalo biologicamente viável ($T \leq T_{min}$ ou $T \geq T_{max}$), assume-se $R(T) = 0$.

### 2. Taxa de Mortalidade (Estágios Imaturos)
A mortalidade dependente da temperatura para os Estágios 1 (ovos a L3) e 2 (L4 a pupa) foi descrita por uma equação polinomial de segundo grau (parábola com concavidade voltada para cima, refletindo o estresse térmico em extremos):

$$d(T) = a_{i} T^2 + b_{i} T + c_{i}$$

### 3. Taxa de Mortalidade com Senescência (Adultos)
A mortalidade na fase adulta (Estágio 3) requer a integração do fator térmico com o desgaste fisiológico natural (idade, $a$). A função implementada garante que a mortalidade tenda ao infinito quando o inseto se aproxima da sua longevidade máxima assumida ($A_3$):

$$\mu_3(a, T) = (a_{3} T^2 + b_{3} T + c_{3}) \cdot d_{3} \frac{a}{(A_{3} - a)^2}$$

Onde $d_3$ é o parâmetro escalar de senescência ajustado iterativamente via rotinas de otimização.

---

## 📂 Estrutura do Repositório

* `Ajuste_Curva_Desenvolvimento.m`: Script que carrega os dados de tempo de desenvolvimento, calcula $R(T)$ e utiliza a função `lsqcurvefit` para extrair os parâmetros $a$, $T_{min}$ e $T_{max}$ de cada estágio imaturo. Gera gráficos comparativos (Dados vs. Ajuste).
* `Ajuste_Curva_Mortalidade.m`: Script para obtenção das taxas vitais de mortalidade. Usa `polyfit` para as componentes térmicas e `lsqcurvefit` para estimar o parâmetro de senescência $d_3$ do estágio adulto. Gera subplots de validação e uma superfície 3D de mortalidade em função de $a$ e $T$.

## 🚀 Como Executar

1. Certifique-se de ter o **MATLAB** instalado em sua máquina.
2. O pacote **Optimization Toolbox** é estritamente necessário, pois os códigos dependem de `lsqcurvefit` para a resolução dos Mínimos Quadrados Não-Lineares.
3. Clone este repositório:
   ```bash
   git clone [https://github.com/SEU-USUARIO/NOME-DO-REPOSITORIO.git](https://github.com/SEU-USUARIO/NOME-DO-REPOSITORIO.git)
