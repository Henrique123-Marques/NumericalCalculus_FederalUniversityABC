import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def f(V, a, b, N, P, T, k):
    termo1 = P + a * (N/V)**2
    termo2 = V - N * b
    return termo1 * termo2 - k * N * T

#Primeira derivada da funcao
def df(V, a, b, N, P):
    termo1 = P + a * (N/V)**2
    termo2 = -2 * a * N**2 / V**3
    return termo1 + (V - N * b) * termo2

#Metodo da bissecao
def bissecao(a_intervalo, b_intervalo, tolerancia, a, b, N, P, T, k):
    # Debug inicial
    st.write(f"Iniciando bisseção com intervalo: [{a_intervalo}, {b_intervalo}]")
    st.write(f"f(a_intervalo) = {f(a_intervalo, a, b, N, P, T, k)}")
    st.write(f"f(b_intervalo) = {f(b_intervalo, a, b, N, P, T, k)}")
    
    if f(a_intervalo, a, b, N, P, T, k) * f(b_intervalo, a, b, N, P, T, k) >= 0:
        raise ValueError("f(a) e f(b) devem ter sinais opostos")
    
    Va = a_intervalo
    Vb = b_intervalo
    iteracao_contador = 0
    
    while (Vb-Va) > tolerancia and iteracao_contador < 1000:  # Adicionado limite de iterações
        Vm = (Va + Vb) / 2
        if f(Vm, a, b, N, P, T, k) == 0:
            return Vm, iteracao_contador
        elif f(Va, a, b, N, P, T, k) * f(Vm, a, b, N, P, T, k) < 0:
            Vb = Vm
        else:
            Va = Vm
        iteracao_contador += 1
    
    return (Va + Vb) / 2, iteracao_contador

#Metodo da falsa-posicao
def falsa_posicao(a_intervalo, b_intervalo, tolerancia, a, b, N, P, T, k):
    if f(a_intervalo, a, b, N, P, T, k) * f(b_intervalo, a, b, N, P, T, k) >= 0:
        raise ValueError("f(a) e f(b) devem ter sinais opostos")

    Va = a_intervalo
    Vb = b_intervalo
    iteracao_contador = 0

    while abs(Vb-Va) > tolerancia and iteracao_contador < 1000:  # Adicionado limite de iterações
        Vm = (Va * f(Vb, a, b, N, P, T, k) - Vb * f(Va, a, b, N, P, T, k)) / (f(Vb, a, b, N, P, T, k) - f(Va, a, b, N, P, T, k))
        if f(Vm, a, b, N, P, T, k) == 0:
            return Vm, iteracao_contador
        elif f(Va, a, b, N, P, T, k) * f(Vm, a, b, N, P, T, k) < 0:
            Vb = Vm
        else:
            Va = Vm
        iteracao_contador += 1
    
    return (Va + Vb) / 2, iteracao_contador

#Metodo de Newton-Raphson
def newton_raphson(V0, tolerancia, a, b, N, P, T, k):
    V = V0
    iteracao_contador = 0

    while iteracao_contador < 1000:  # Adicionado limite de iterações
        f_V = f(V, a, b, N, P, T, k)
        df_V = df(V, a, b, N, P)
        
        if abs(df_V) < 1e-10:  # Evitar divisão por zero
            raise ValueError("Derivada muito próxima de zero")
            
        V_novo = V - f_V / df_V

        if abs(V_novo - V) < tolerancia:
            return V_novo, iteracao_contador
        V = V_novo
        iteracao_contador += 1
    
    raise ValueError("Método não convergiu após 1000 iterações")

# Interface Streamlit
st.title('Equação de Estado do CO₂')
st.write("""Enunciado do Exercicio 8 - Para o dióxido de carbono (CO2), a equação de estado de um gás tem os
seguintes coeficientes: a = 0, 401Pa m^3 e b = 42,7.10^−6 m^3. Determine o volume ocupado por 1000 moléculas deste gás, 
à temperatura de 300K e pressão de 3, 5.10^7Pa, pelo método de bisecção, com uma tolerância de 10^−12. Constante de Boltzmann é k = 1, 3806503*10^−23 Joule K^−1.
[p + a*(N/V )^2]*(V − N*b) = k*N*T. Faça também outros algoritmos para resolução usando os métodos de falsa-
posição e de Newton-Raphson.

Este aplicativo calcula o volume ocupado por moléculas de CO₂ usando diferentes métodos numéricos:
- Método da Bisseção
- Método da Falsa Posição
- Método de Newton-Raphson""")

# Sidebar para parâmetros
st.sidebar.header('Parâmetros')

# Parâmetros principais com valores padrão
a = st.sidebar.number_input('Coeficiente a (Pa m³)', value=0.401, format='%f')
b = st.sidebar.number_input('Coeficiente b (m³)', value=42.7e-6, format='%e')
N = st.sidebar.number_input('Número de moléculas', value=1000, min_value=1)
P = st.sidebar.number_input('Pressão (Pa)', value=3.5e7, format='%e')
T = st.sidebar.number_input('Temperatura (K)', value=300.0)
k = st.sidebar.number_input('Constante de Boltzmann (J/K)', value=1.3806503e-23, format='%e')
tolerancia = st.sidebar.number_input('Tolerância', value=1e-12, format='%e')

# Intervalos para os métodos
Va = N * b  # próximo de N*b
Vb = Va * 1.001  # Aumentado o intervalo inicial
V0 = (Va + Vb) / 2  # chute inicial para Newton-Raphson

# Debug dos valores iniciais
st.write("### Valores iniciais:")
st.write(f"Intervalo inicial: Va = {Va:.10e}, Vb = {Vb:.10e}")
st.write(f"Chute inicial para Newton-Raphson: V0 = {V0:.10e}")

if st.button('Calcular'):
    try:
        with st.spinner('Calculando...'):
            # Cálculos
            V_bissecao, iteracao_bissecao = bissecao(Va, Vb, tolerancia, a, b, N, P, T, k)
            V_falsa, iteracao_falsa = falsa_posicao(Va, Vb, tolerancia, a, b, N, P, T, k)
            V_newton, iteracao_newton = newton_raphson(V0, tolerancia, a, b, N, P, T, k)

            # Exibindo resultados em uma tabela
            resultados = {
                'Método': ['Bisseção', 'Falsa Posição', 'Newton-Raphson'],
                'Volume (m³)': [V_bissecao, V_falsa, V_newton],
                'Iterações': [iteracao_bissecao, iteracao_falsa, iteracao_newton],
                'f(V)': [f(V_bissecao, a, b, N, P, T, k), 
                         f(V_falsa, a, b, N, P, T, k), 
                         f(V_newton, a, b, N, P, T, k)]
            }
            
            st.write("### Resultados")
            st.dataframe(resultados)

            # Gráfico
            V_media = (V_bissecao + V_falsa + V_newton) / 3
            delta_V = abs(V_media) * 0.0001
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
            ax.set_title('Equação de Estado do CO₂ - Zoom nos Resultados')
            ax.legend()
            plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
            st.pyplot(fig)

            # Diferenças entre os métodos
            st.write("### Diferenças entre os métodos")
            st.write(f"Bisseção - Falsa Posição: {(V_bissecao - V_falsa):.15e} m³")
            st.write(f"Bisseção - Newton-Raphson: {(V_bissecao - V_newton):.15e} m³")
            st.write(f"Falsa Posição - Newton-Raphson: {(V_falsa - V_newton):.15e} m³")

    except Exception as e:
        st.error(f"Erro nos cálculos: {str(e)}")
        st.write("### Debug de valores:")
        st.write(f"f(Va) = {f(Va, a, b, N, P, T, k)}")
        st.write(f"f(Vb) = {f(Vb, a, b, N, P, T, k)}")