import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Funções matemáticas
def f(V, a, b, N, P, T, k):
    termo1 = P + a * (N/V)**2
    termo2 = V - N * b
    return termo1 * termo2 - k * N * T

def df(V, a, b, N, P):
    termo1 = P + a * (N/V)**2
    termo2 = -2 * a * N**2 / V**3
    return termo1 + (V - N * b) * termo2

# Método da Bisseção
def bissecao(a_intervalo, b_intervalo, tolerancia, a, b, N, P, T, k):
    if f(a_intervalo, a, b, N, P, T, k) * f(b_intervalo, a, b, N, P, T, k) >= 0:
        raise ValueError("f(a) e f(b) devem ter sinais opostos")
    
    Va = a_intervalo
    Vb = b_intervalo
    iteracao_contador = 0
    
    while (Vb - Va) > tolerancia and iteracao_contador < 1000:
        Vm = (Va + Vb) / 2
        if f(Vm, a, b, N, P, T, k) == 0:
            return Vm, iteracao_contador
        elif f(Va, a, b, N, P, T, k) * f(Vm, a, b, N, P, T, k) < 0:
            Vb = Vm
        else:
            Va = Vm
        iteracao_contador += 1
    
    return (Va + Vb) / 2, iteracao_contador

# Método da Falsa Posição
def falsa_posicao(a_intervalo, b_intervalo, tolerancia, a, b, N, P, T, k):
    if f(a_intervalo, a, b, N, P, T, k) * f(b_intervalo, a, b, N, P, T, k) >= 0:
        raise ValueError("f(a) e f(b) devem ter sinais opostos")

    Va = a_intervalo
    Vb = b_intervalo
    iteracao_contador = 0

    while abs(Vb - Va) > tolerancia and iteracao_contador < 1000:
        Vm = (Va * f(Vb, a, b, N, P, T, k) - Vb * f(Va, a, b, N, P, T, k)) / (f(Vb, a, b, N, P, T, k) - f(Va, a, b, N, P, T, k))
        if f(Vm, a, b, N, P, T, k) == 0:
            return Vm, iteracao_contador
        elif f(Va, a, b, N, P, T, k) * f(Vm, a, b, N, P, T, k) < 0:
            Vb = Vm
        else:
            Va = Vm
        iteracao_contador += 1
    
    return (Va + Vb) / 2, iteracao_contador

# Método de Newton-Raphson
def newton_raphson(V0, tolerancia, a, b, N, P, T, k):
    V = V0
    iteracao_contador = 0

    while iteracao_contador < 1000:
        f_V = f(V, a, b, N, P, T, k)
        df_V = df(V, a, b, N, P)
        
        if abs(df_V) < 1e-10:
            raise ValueError("Derivada muito próxima de zero")
            
        V_novo = V - f_V / df_V

        if abs(V_novo - V) < tolerancia:
            return V_novo, iteracao_contador
        V = V_novo
        iteracao_contador += 1
    
    raise ValueError("Método não convergiu após 1000 iterações")

# Interface Streamlit
st.set_page_config(page_title="Cálculo do Volume do CO₂", layout="wide")
st.title("✨ Cálculo do Volume Ocupado por Moléculas de CO₂ ✨")
st.markdown("""
Este aplicativo resolve a equação de estado do dióxido de carbono (CO₂) para determinar o volume ocupado por 1000 moléculas, utilizando métodos numéricos clássicos: **Bisseção**, **Falsa Posição** e **Newton-Raphson**. Os parâmetros são fixos conforme o enunciado do Exercício 8, e os resultados são apresentados com gráficos e análises detalhadas.
""")

# Sidebar com parâmetros fixos
st.sidebar.header("🔧 Parâmetros Fixos do CO₂")
a = st.sidebar.number_input('Coeficiente a (Pa m³)', value=0.401, format='%f', disabled=True)
b = st.sidebar.number_input('Coeficiente b (m³)', value=42.7e-6, format='%e', disabled=True)
N = st.sidebar.number_input('Número de moléculas', value=1000, min_value=1, disabled=True)
P = st.sidebar.number_input('Pressão (Pa)', value=3.5e7, format='%e', disabled=True)
T = st.sidebar.number_input('Temperatura (K)', value=300.0, disabled=True)
k = st.sidebar.number_input('Constante de Boltzmann (J/K)', value=1.3806503e-23, format='%e', disabled=True)
tolerancia = st.sidebar.number_input('Tolerância', value=1e-12, format='%e', disabled=True)

# Intervalos iniciais
Va = N * b
Vb = Va * 1.001
V0 = (Va + Vb) / 2

# Explicação detalhada
st.markdown("### Passo a Passo da Resolucão do Exercicio")
st.write("""
Para resolver o Exercicio 8 seguimos a equação de estado do CO2 dada por P + a N/V^2 vezes V - N vezes b = k vezes N vezes T. Nosso objetivo e encontrar o volume V ocupado por 1000 moleculas de CO2 com os valores fixos fornecidos: a = 0.401 Pa m^3, b = 42.7 vezes 10^-6 m^3, N = 1000, P = 3.5 vezes 10^7 Pa, T = 300 K, k = 1.3806503 vezes 10^-23 J/K e tolerancia de 10^-12. Vamos resolver isso passo a passo com os tres metodos pedidos.

1. Reorganizacao da Equacao
Primeiro reorganizamos a equacao para a forma f(V) = 0: f(V) = P + a N/V^2 vezes V - N vezes b - k vezes N vezes T = 0. Essa funcao sera usada para encontrar a raiz V que e o volume procurado.

2. Metodo da Bissecao
O metodo da bissecao requer um intervalo inicial Va e Vb onde f(Va) e f(Vb) possuem sinais opostos. Escolhemos Va = N vezes b = 1000 vezes 42.7 vezes 10^-6 = 4.27 vezes 10^-2 m^3 como ponto proximo do limite fisico onde V - N vezes b = 0 e Vb = 1.001 vezes Va para garantir um intervalo pequeno mas suficiente. Iteramos dividindo o intervalo ao meio ate que a diferenca seja menor que 10^-12.

3. Metodo da Falsa Posicao
Similar a bissecao usamos o mesmo intervalo inicial. Porem em vez de dividir o intervalo ao meio calculamos um ponto Vm pela formula Vm = Va vezes f(Vb) - Vb vezes f(Va) dividido por f(Vb) - f(Va). Atualizamos Va ou Vb com base no sinal de f(Vm) ate atingir a tolerancia.

4. Metodo de Newton-Raphson
Este metodo requer um chute inicial V0 = Va + Vb dividido por 2 e a derivada f'(V): f'(V) = P + a N/V^2 + V - N vezes b vezes -2 vezes a vezes N^2/V^3. Iteramos com Vnovo = V - f(V) dividido por f'(V) ate que a diferenca entre iteracoes seja menor que 10^-12.

5. Analise e Comparacao
Cada metodo converge para um volume proximo mas com diferencas sutis devido as suas abordagens. A bissecao e robusta mas lenta, a falsa posicao e mais rapida em intervalos bem definidos e Newton-Raphson converge rapidamente com um bom chute inicial. Os resultados sao validados pelo grafico de f(V) e pela proximidade dos valores encontrados.
""")

# Botão para calcular
if st.button('🔍 Calcular Resultados'):
    try:
        with st.spinner('Processando os cálculos...'):
            V_bissecao, iteracao_bissecao = bissecao(Va, Vb, tolerancia, a, b, N, P, T, k)
            V_falsa, iteracao_falsa = falsa_posicao(Va, Vb, tolerancia, a, b, N, P, T, k)
            V_newton, iteracao_newton = newton_raphson(V0, tolerancia, a, b, N, P, T, k)

            st.markdown("### 🎯 Resultados Obtidos")
            resultados = {
                'Método': ['Bisseção', 'Falsa Posição', 'Newton-Raphson'],
                'Volume (m³)': [f"{V_bissecao:.15e}", f"{V_falsa:.15e}", f"{V_newton:.15e}"],
                'Iterações': [iteracao_bissecao, iteracao_falsa, iteracao_newton],
                'f(V)': [f"{f(V_bissecao, a, b, N, P, T, k):.15e}", 
                         f"{f(V_falsa, a, b, N, P, T, k):.15e}", 
                         f"{f(V_newton, a, b, N, P, T, k):.15e}"]
            }
            st.dataframe(resultados)

            # Gráfico
            V_media = (V_bissecao + V_falsa + V_newton) / 3
            delta_V = abs(V_media) * 0.001
            V_range_zoom = np.linspace(V_media - delta_V, V_media + delta_V, 1000)
            f_values_zoom = [f(V, a, b, N, P, T, k) for V in V_range_zoom]

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(V_range_zoom, f_values_zoom, 'b-', label='f(V)')
            ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
            ax.plot(V_bissecao, f(V_bissecao, a, b, N, P, T, k), 'ro', label='Bisseção')
            ax.plot(V_falsa, f(V_falsa, a, b, N, P, T, k), 'go', label='Falsa Posição')
            ax.plot(V_newton, f(V_newton, a, b, N, P, T, k), 'yo', label='Newton-Raphson')
            ax.grid(True)
            ax.set_xlabel('Volume (m³)')
            ax.set_ylabel('f(V)')
            ax.set_title('Comportamento da Equação de Estado do CO₂')
            ax.legend()
            plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
            st.pyplot(fig)

            st.markdown("### 🔬 Análise das Diferenças")
            st.write(f"Diferença Bisseção - Falsa Posição: {(V_bissecao - V_falsa):.15e} m³")
            st.write(f"Diferença Bisseção - Newton-Raphson: {(V_bissecao - V_newton):.15e} m³")
            st.write(f"Diferença Falsa Posição - Newton-Raphson: {(V_falsa - V_newton):.15e} m³")

    except Exception as e:
        st.error(f"Erro ao calcular: {str(e)}")

# Seção de referências
st.markdown("### 📚 Referências Citadas")
st.write("""
- Stewart, J. (2016). *Calculus: Early Transcendentals*. 8ª ed. Cengage Learning;
- Burden, R. L., & Faires, J. D. (2011). *Numerical Analysis*. 9ª ed. Brooks/Cole.
""")