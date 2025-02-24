import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

# Configuração inicial para multipáginas
st.set_page_config(page_title="Cálculo Numérico - Lista 1", layout="wide")

# Funções do Exercício 8 (CO₂)
def f(V, a, b, N, P, T, k):
    termo1 = P + a * (N/V)**2
    termo2 = V - N * b
    return termo1 * termo2 - k * N * T

def df(V, a, b, N, P):
    termo1 = P + a * (N/V)**2
    termo2 = -2 * a * N**2 / V**3
    return termo1 + (V - N * b) * termo2

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

# Função da Sucessão (Exercício 3 - original)
def calcular_sucessao(tol=1e-1, n_max=16):
    I = [(1 / math.e) * (math.e - 1)]  # I0
    for n in range(n_max):
        I_next = 1 - (n + 1) * I[-1]
        I.append(I_next)
        if abs(I_next) < tol:  # Critério de convergência
            break
    return I

# Estrutura multipágina
pages = {
    "Exercício 1": "ex1",
    "Exercício 2": "ex2",
    "Exercício 3": "ex3",
    "Exercício 4": "ex4",
    "Exercício 5": "ex5",
    "Exercício 6": "ex6",
    "Exercício 7": "ex7",
    "Exercício 8": "ex8",
    "Exercício 9": "ex9"
}

page = st.sidebar.selectbox("Escolha um exercício", list(pages.keys()))

# Exercício 1 - IEEE 754
if page == "Exercício 1":
    st.title("🌟 Explorando o Padrão IEEE 754 - Exercício 1")
    st.markdown("""
    Este aplicativo investiga as propriedades do padrão IEEE 754 para números de ponto flutuante em 64 bits, calculando o maior e menor número representáveis, o epsilon da máquina e analisando a expressão (1 + x - 1) / x para diferentes valores de x. Vamos mergulhar no fascinante mundo da precisão numérica!
    """)

    st.markdown("### 📘 Explicação do Problema")
    st.write("""
    O padrão IEEE 754 define como números de ponto flutuante sao representados em computadores. Em 64 bits, usamos 1 bit para o sinal, 11 bits para o expoente e 52 bits para a mantissa. Com isso, podemos determinar:

    - Maior numero representavel: (2 - epsilon) vezes 2^1023
    - Menor numero normalizado: 1.0 vezes 2^-1022
    - Epsilon da maquina: o menor valor u tal que 1 + u e diferente de 1

    Alem disso, analisamos a expressao (1 + x - 1) / x para x = 10^-15 e x = 10^15, comparando o valor aproximado com o exato (1) e calculando erros absoluto e relativo. Por fim, apresentamos uma solucao alternativa simplificando a expressao para 1.
    """)

    if st.button("🔍 Executar Cálculos"):
        with st.spinner("Calculando os valores..."):
            Valormax = (2 - np.finfo(float).eps) * 2**1023
            st.markdown("### 🎯 Resultados")
            st.write(f"Maior número representável (Valormax): {Valormax:e}")

            Valormin = 1.0 * 2**-1022
            st.write(f"Menor número representável normalizado (Valormin): {Valormin:e}")

            u = 1.0
            while 1.0 + u != 1.0:
                u *= 0.5
            epsilon = 2 * u
            st.write(f"Epsilon da máquina: {epsilon:e}")

            Valores_X = [1e-15, 1e+15]
            Valor_exato = 1
            st.markdown("### 🔬 Análise da Expressão (1 + x - 1) / x")
            for x in Valores_X:
                f_expre = (1 + x - 1) / x
                erro_absoluto = abs(Valor_exato - f_expre)
                erro_relativo = erro_absoluto / abs(Valor_exato)
                st.write(f"\nPara x = {x:e}:")
                st.write(f"Valor aproximado: {f_expre:e}")
                st.write(f"Erro absoluto: {erro_absoluto:e}")
                st.write(f"Erro relativo: {erro_relativo:e}")

            st.markdown("### 🌟 Solução Alternativa")
            for x in Valores_X:
                f_expre = 1
                erro_absoluto = abs(Valor_exato - f_expre)
                erro_relativo = erro_absoluto / abs(Valor_exato)
                st.write(f"\nPara x = {x:e}:")
                st.write(f"Valor aproximado: {f_expre:e}")
                st.write(f"Erro absoluto: {erro_absoluto:e}")
                st.write(f"Erro relativo: {erro_relativo:e}")

# Exercício 2 - Raiz da função f(x) com bisseção e gráfico
elif page == "Exercício 2":
    st.title("✨ Análise da Função Polinomial - Exercício 2")
    st.markdown("""
    Este aplicativo avalia a função \( f(x) = x^7 - 7x^6 + 21x^5 - 35x^4 + 35x^3 - 21x^2 + 7x - 1 \) em um pequeno intervalo ao redor de \( x = 1 \) usando o método da bisseção para encontrar a raiz e exibe seu comportamento graficamente.
    """)

    st.markdown("### 📘 Explicação do Problema")
    st.write("""
    A função dada é um polinômio de grau 7: \( f(x) = x^7 - 7x^6 + 21x^5 - 35x^4 + 35x^3 - 21x^2 + 7x - 1 \). O objetivo é encontrar a raiz no intervalo \( [1 - 2 \times 10^{-8}, 1 + 2 \times 10^{-8}] \) utilizando o método da bisseção em duas implementações diferentes. O método da bisseção divide repetidamente o intervalo ao meio, selecionando o subintervalo onde ocorre uma mudança de sinal, até que o tamanho do intervalo seja menor que a tolerância (\( 10^{-10} \)) ou \( f(x) \) seja suficientemente pequeno. Além disso, o comportamento da função é visualizado em um gráfico.
    """)

    # Função f(x)
    def f(x):
        return x**7 - 7*x**6 + 21*x**5 - 35*x**4 + 35*x**3 - 21*x**2 + 7*x - 1

    # Método da Bisseção - Código 1
    def bissecao_1(f, a, b, tol=1e-10, max_iter=100):
        if f(a) * f(b) >= 0:
            st.error("ERRO! f(a) e f(b) precisam ter sinais opostos. O intervalo não contém uma raiz.")
            return None
        for i in range(max_iter):
            m = (a + b) / 2
            if (b - a) / 2 < tol or abs(f(m)) < tol:
                return m
            if f(m) * f(a) < 0:
                b = m
            else:
                a = m
        return (a + b) / 2

    # Método da Bisseção - Código 2
    def bissecao_2(f, a, b, tol=1e-10, max_iter=100):
        if f(a) * f(b) >= 0:
            st.error("ERRO! f(a) e f(b) precisam ter sinais opostos. O intervalo não contém uma raiz.")
            return None
        iter = 0
        while (b - a) / 2 > tol and iter < max_iter:
            m = (a + b) / 2
            if f(m) == 0:
                return m
            if f(a) * f(m) < 0:
                b = m
            else:
                a = m
            iter += 1
        return (a + b) / 2

    if st.button("🔍 Calcular Raiz e Gerar Gráfico"):
        with st.spinner("Calculando a raiz e gerando o gráfico..."):
            # Intervalo
            a = 1 - 2 * 10**-8
            b = 1 + 2 * 10**-8

            # Aplicação dos métodos
            raiz_b1 = bissecao_1(f, a, b)
            raiz_b2 = bissecao_2(f, a, b)

            st.markdown("### 🎯 Resultados do Método da Bisseção")
            if raiz_b1 is not None:
                st.write(f"Raiz encontrada pelo método da bisseção (Código 1): {raiz_b1:.15e}")
                st.write(f"Valor de f(raiz) (Código 1): {f(raiz_b1):.15e}")
            if raiz_b2 is not None:
                st.write(f"Raiz encontrada pelo método da bisseção (Código 2): {raiz_b2:.15e}")
                st.write(f"Valor de f(raiz) (Código 2): {f(raiz_b2):.15e}")

            # Gráfico
            valores_x = np.linspace(a, b, 401)
            valores_y = f(valores_x)

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(valores_x, valores_y, label="f(x)")
            ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
            if raiz_b1 is not None:
                ax.plot(raiz_b1, f(raiz_b1), 'ro', label='Raiz')
            if raiz_b2 is not None:
                ax.plot(raiz_b2, f(raiz_b2), 'go')
            ax.set_title("Comportamento da função f(x)", fontsize=20)
            ax.set_xlabel("x", fontsize=15)
            ax.set_ylabel("f(x)", fontsize=15)
            ax.grid(True)
            ax.legend()
            st.pyplot(fig)

# Exercício 3 - Sucessão
elif page == "Exercício 3":
    st.title("🌟 Análise da Sucessão Recursiva - Exercício 3")
    st.markdown("""
    Este aplicativo calcula e analisa uma sucessão definida pela fórmula recursiva \( I_{n+1} = 1 - (n + 1) I_n \), com \( I_0 = (1/e) (e - 1) \). Vamos explorar os valores da sequência e verificar sua convergência com um critério de tolerância!
    """)

    st.markdown("### 📘 Explicação do Problema")
    st.write("""
    A sucessao e definida por I0 = 1/e vezes (e - 1) e para n maior ou igual a 0, In+1 = 1 - (n + 1) vezes In. O objetivo e calcular os termos da sucessao ate um maximo de 16 iteracoes ou ate que o valor absoluto de In seja menor que a tolerancia 10^-1. Depois, verificamos se a sucessao converge para 0 com precisao de 10^-10. Vamos calcular e analisar os resultados passo a passo!
    """)

    if st.button("🔍 Calcular Sucessão"):
        with st.spinner("Calculando os termos da sucessão..."):
            sucessao = calcular_sucessao(tol=1e-1, n_max=16)
            
            st.markdown("### 🎯 Resultados da Sucessão")
            st.write("Valores da sucessao:")
            for n, val in enumerate(sucessao):
                st.write(f"I{n} = {val:e}")

            st.markdown("### 🔬 Verificação de Convergência")
            if abs(sucessao[-1]) < 1e-10:
                st.write(f"A sucessão converge para 0 após {len(sucessao) - 1} iterações.")
            else:
                st.write("A sucessão não atingiu a precisão desejada (10^-10).")
                st.write(f"Último valor (I{len(sucessao)-1}): {sucessao[-1]:e}")

    st.markdown("### 📚 Referências Citadas")
    st.write("""
    - Burden, R. L., & Faires, J. D. (2011). *Numerical Analysis*. 9ª ed. Brooks/Cole.
    """)

# Exercício 4 - Estimativa de Pi
elif page == "Exercício 4":
    st.title("✨ Estimativa de π - Exercício 4")
    st.markdown("""Este aplicativo resolve duas partes do Exercício 4: estima o valor de π usando o método de Monte Carlo!""")

    st.markdown("### 📘 Explicação do Problema")
    st.write("""Estimativa de pi: Usamos o metodo de Monte Carlo gerando pontos aleatorios em um quadrado de lado 1 e verificando quantos caem dentro de um quarto de circulo. A formula e pi estimado = 4 vezes m dividido por n, onde m e o numero de pontos dentro do circulo e n e o total de pontos. Calculamos o erro para diferentes valores de n.""")

    if st.button("🔍 Calcular Resultados"):
        with st.spinner("Processando os cálculos..."):
            
            def estimar_pi(n):
                x = np.random.uniform(0, 1, n)
                y = np.random.uniform(0, 1, n)
                dentro_do_circulo = (x**2 + y**2) <= 1
                m = np.sum(dentro_do_circulo)
                pi_estimado = 4 * m / n
                return pi_estimado

            def calcular_erro(n_valores):
                erros = []
                valores_pi = []
                for n in n_valores:
                    pi_estimado = estimar_pi(n)
                    valores_pi.append(pi_estimado)
                    erro = abs(np.pi - pi_estimado)
                    erros.append(erro)
                return valores_pi, erros

            n_valores = np.logspace(2, 6, num=20, dtype=int)
            valores_pi, erros = calcular_erro(n_valores)

            st.markdown("### 🎯 Resultados da Estimativa de π")
            st.write("Valores estimados de pi e erros para diferentes n:")
            for n, pi_est, erro in zip(n_valores, valores_pi, erros):
                st.write(f"n = {n}: π estimado = {pi_est:.6f}, Erro = {erro:.6f}")

            # Gráfico do erro
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.plot(n_valores, erros, marker='o', linestyle='-', label='Erro absoluto')
            ax.axhline(y=0, color='r', linestyle='--', label='Erro zero')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel('Número de pontos (n)')
            ax.set_ylabel('Erro absoluto')
            ax.set_title('Evolução do erro na estimativa de π')
            ax.legend()
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)
            st.pyplot(fig)

# Exercício 5 - Vazio
elif page == "Exercício 5":
    st.title("✨ Exercício 5")
    st.markdown("""
    Esta página está reservada para o Exercício 5. No momento, não há conteúdo disponível!
    """)

# Exercício 6 - Comprimento Máximo da Barra
elif page == "Exercício 6":
    st.title("✨ Exercício 6 - Comprimento Máximo da Barra")
    st.markdown("""Esta página calcula o comprimento máximo de uma barra usando os métodos de Newton-Raphson e Bisseção!""")
    st.markdown("""📘 Explicação do Problema""")
    st.write("""Comprimento máximo da barra: Resolvemos a equação f(alpha) = l2 vezes cos(pi - gamma - alpha) dividido por sen(pi - gamma - alpha)^2 - l1 vezes cos(alpha) dividido por sen(alpha)^2 = 0, com l1 = 8, l2 = 10 e gamma = 3 vezes pi dividido por 5. Usamos Newton-Raphson e Bisseção para encontrar alpha e depois calculamos L = l2 dividido por sen(pi - gamma - alpha) + l1 dividido por sen(alpha).
   """)

    if st.button("Calcular Resultados"):
        with st.spinner('Processando os cálculos...'):
            l2 = 10
            l1 = 8
            gamma = 3 * np.pi / 5

            def func(alpha, l1=l1, l2=l2, gamma=gamma):
                return (l2 * np.cos(np.pi - gamma - alpha) / np.sin(np.pi - gamma - alpha)**2 -
                        l1 * np.cos(alpha) / np.sin(alpha)**2)

            def dfunc(alpha, l1=l1, l2=l2, gamma=gamma):
                return (-l2 * (np.sin(np.pi - gamma - alpha) + 2 * np.cos(np.pi - gamma - alpha) / np.sin(np.pi - gamma - alpha)**3) -
                        l1 * (-np.sin(alpha) + 2 * np.cos(alpha) / np.sin(alpha)**3))

            # Newton-Raphson
            alpha0 = 0.5
            tolerancia = 1e-6
            alpha_newton, iter_newton = newton_raphson(alpha0, tolerancia, l1, l2, 0, 0, 0, 0)  # Parâmetros extras não usados
            L_newton = l2 / np.sin(np.pi - gamma - alpha_newton) + l1 / np.sin(alpha_newton)

            # Bisseção
            alpha_bisect, iter_bissect = bissecao(0.1, np.pi/2, tolerancia, l1, l2, 0, 0, 0, 0)  # Parâmetros extras não usados
            L_bisect = l2 / np.sin(np.pi - gamma - alpha_bisect) + l1 / np.sin(alpha_bisect)

            st.markdown("### 🎯 Resultados do Comprimento Máximo da Barra")
            st.write(f"Alpha encontrado por Newton-Raphson: {alpha_newton:.6f}")
            st.write(f"Comprimento máximo L (Newton-Raphson): {L_newton:.6f}")
            st.write(f"Alpha encontrado por Bisseção: {alpha_bisect:.6f}")
            st.write(f"Comprimento máximo L (Bisseção): {L_bisect:.6f}")

# Exercício 7 - Vazio
elif page == "Exercício 7":
    st.title("✨ Exercício 7")
    st.markdown("""
    Esta página está reservada para o Exercício 7. No momento, não há conteúdo disponível!
    """)

# Exercício 8 - CO₂
elif page == "Exercício 8":
    st.title("✨ Cálculo do Volume Ocupado por Moléculas de CO₂ - Exercício 8")
    st.markdown("""
    Este aplicativo resolve a equação de estado do dióxido de carbono (CO₂) para determinar o volume ocupado por 1000 moléculas, utilizando métodos numéricos clássicos: **Bisseção**, **Falsa Posição** e **Newton-Raphson**. Os parâmetros são fixos conforme o enunciado do Exercício 8, e os resultados são apresentados com gráficos e análises detalhadas.
    """)

    st.sidebar.header("🔧 Parâmetros Fixos do CO₂")
    a = st.sidebar.number_input('Coeficiente a (Pa m³)', value=0.401, format='%f', disabled=True)
    b = st.sidebar.number_input('Coeficiente b (m³)', value=42.7e-6, format='%e', disabled=True)
    N = st.sidebar.number_input('Número de moléculas', value=1000, min_value=1, disabled=True)
    P = st.sidebar.number_input('Pressão (Pa)', value=3.5e7, format='%e', disabled=True)
    T = st.sidebar.number_input('Temperatura (K)', value=300.0, disabled=True)
    k = st.sidebar.number_input('Constante de Boltzmann (J/K)', value=1.3806503e-23, format='%e', disabled=True)
    tolerancia = st.sidebar.number_input('Tolerância', value=1e-12, format='%e', disabled=True)

    Va = N * b
    Vb = Va * 1.001
    V0 = (Va + Vb) / 2

    st.markdown("### 📘 Passo a Passo da Resolucão do Exercicio")
    st.write("""
    Para resolver o Exercicio 8 seguimos a equação de estado do CO2 dada por P + a N/V^2 * V - N * b = k * N * T. Nosso objetivo e encontrar o volume V ocupado por 1000 moleculas de CO2 com os valores fixos fornecidos: a = 0.401 Pa m^3, b = 42.7 * 10^-6 m^3, N = 1000, P = 3.5 * 10^7 Pa, T = 300 K, k = 1.3806503 * 10^-23 J/K e tolerancia de 10^-12. Vamos resolver isso passo a passo com os tres metodos pedidos.

    1. Reorganização da Equação:
    Primeiro reorganizamos a equacao para a forma f(V) = 0: f(V) = P + a N/V^2 * V - N * b - k * N * T = 0. Essa funcao sera usada para encontrar a raiz V que e o volume procurado.

    2. Metodo da Bisseção
    O metodo da bissecao requer um intervalo inicial Va e Vb onde f(Va) e f(Vb) possuem sinais opostos. Escolhemos Va = N * b = 1000 * 42.7 * 10^-6 = 4.27 * 10^-2 m^3 como ponto proximo do limite fisico onde V - N * b = 0 e Vb = 1.001 * Va para garantir um intervalo pequeno mas suficiente. Iteramos dividindo o intervalo ao meio ate que a diferenca seja menor que 10^-12.

    3. Metodo da Falsa Posição
    Similar a bissecao usamos o mesmo intervalo inicial. Porem em vez de dividir o intervalo ao meio calculamos um ponto Vm pela formula Vm = Va * f(Vb) - Vb * f(Va) dividido por f(Vb) - f(Va). Atualizamos Va ou Vb com base no sinal de f(Vm) ate atingir a tolerancia.

    4. Metodo de Newton-Raphson
    Este metodo requer um chute inicial V0 = Va + Vb dividido por 2 e a derivada f'(V): f'(V) = P + a N/V^2 + V - N * b * -2 * a * N^2/V^3. Iteramos com Vnovo = V - f(V) dividido por f'(V) ate que a diferenca entre iteracoes seja menor que 10^-12.

    5. Analise e Comparação
    Cada metodo converge para um volume proximo mas com diferencas sutis devido as suas abordagens. A bissecao e robusta mas lenta, a falsa posicao e mais rapida em intervalos bem definidos e Newton-Raphson converge rapidamente com um bom chute inicial. Os resultados sao validados pelo grafico de f(V) e pela proximidade dos valores encontrados.
    """)

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

                V_media = (V_bissecao + V_falsa + V_newton) / 3
                delta_V = abs(V_media) * 0.001
                V_range_zoom = np.linspace(V_media - delta_V, V_media + delta_V, 1000)
                f_values_zoom = [f(V, a, b, N, P, T, k) for V in V_range_zoom]

                fig, ax = plt.subplots(figsize=(12, 6))
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

    st.markdown("### 📚 Referências Citadas")
    st.write("""
    - Stewart, J. (2016). *Calculus: Early Transcendentals*. 8ª ed. Cengage Learning.
    - Burden, R. L., & Faires, J. D. (2011). *Numerical Analysis*. 9ª ed. Brooks/Cole.
    """)

# Exercício 9 - Resolução de f(x) = -1/x³ - 1/x² + 2
elif page == "Exercício 9":
    st.title("✨ Resolução de Equação Não Linear - Exercício 9")
    st.markdown("""
    Este aplicativo resolve a equação \( f(x) = -\frac{1}{x^3} - \frac{1}{x^2} + 2 = 0 \) utilizando os métodos numéricos **Bisseção**, **Falsa Posição** e **Newton-Raphson**, com intervalo inicial \([-2.0, -0.5]\) e um critério de parada baseado no epsilon da máquina.
    """)

    st.markdown("### 📘 Explicação do Problema")
    st.write("""
    A função \( f(x) = -\frac{1}{x^3} - \frac{1}{x^2} + 2 \) deve ser resolvida para encontrar uma raiz no intervalo \([-2.0, -0.5]\). Os métodos numéricos utilizados são:

    - **Bisseção**: Divide o intervalo ao meio iterativamente até que o erro relativo seja menor que o epsilon da máquina.
    - **Falsa Posição**: Usa uma interpolação linear para estimar a raiz, ajustando o intervalo com base nos sinais.
    - **Newton-Raphson**: Utiliza a derivada \( f'(x) = \frac{3}{x^4} + \frac{2}{x^3} \) e um chute inicial \( x_0 = -0.75 \) para convergir rapidamente à raiz.

    O critério de parada é baseado no erro relativo: \( |x_{novo} - x_{velho}| \leq \epsilon \cdot \max(1, |x_{novo}|) \), onde \(\epsilon\) é o epsilon da máquina.
    """)

    # Funções do Exercício 9
    eps = np.finfo(float).eps

    def f(x):
        return -1/(x**3) - 1/(x**2) + 2

    def df(x):
        return 3/(x**4) + 2/(x**3)

    def erro_relativo(x_novo, x_velho):
        return abs(x_novo - x_velho) <= eps * max(1, abs(x_novo))

    def bissecao_9(a, b):
        if f(a) * f(b) >= 0:
            raise ValueError("f(a) e f(b) devem ter sinais opostos")
        x_velho = a
        iter_count = 0
        while True:
            x_novo = (a + b) / 2
            if erro_relativo(x_novo, x_velho):
                return x_novo, iter_count
            elif f(a) * f(x_novo) < 0:
                b = x_novo
            else:
                a = x_novo
            x_velho = x_novo
            iter_count += 1

    def falsa_posicao_9(a, b):
        if f(a) * f(b) >= 0:
            raise ValueError("f(a) e f(b) devem ter sinais opostos")
        x_velho = a
        iter_count = 0
        while True:
            x_novo = (a * f(b) - b * f(a)) / (f(b) - f(a))
            if erro_relativo(x_novo, x_velho):
                return x_novo, iter_count
            elif f(a) * f(x_novo) < 0:
                b = x_novo
            else:
                a = x_novo
            x_velho = x_novo
            iter_count += 1

    def newton_raphson_9(x0):
        x_velho = x0
        iter_count = 0
        while True:
            x_novo = x_velho - f(x_velho) / df(x_velho)
            if erro_relativo(x_novo, x_velho):
                return x_novo, iter_count
            x_velho = x_novo
            iter_count += 1

    if st.button("🔍 Calcular Resultados"):
        with st.spinner("Calculando as raízes..."):
            # Parâmetros iniciais
            a = -2.0
            b = -0.5
            x0 = -0.75

            st.write(f"Epsilon da máquina: {eps:.2e}")
            st.write(f"f({a}) = {f(a):.2f}, f({b}) = {f(b):.2f}")

            # Executando os métodos
            try:
                x_biss, iter_biss = bissecao_9(a, b)
                st.markdown("### 🎯 Resultados")
                st.write(f"Bisseção: x = {x_biss:.15f}, Iterações = {iter_biss}")
            except ValueError as e:
                st.error(f"Erro na bisseção: {e}")
                x_biss = None

            try:
                x_falsa, iter_falsa = falsa_posicao_9(a, b)
                st.write(f"Falsa Posição: x = {x_falsa:.15f}, Iterações = {iter_falsa}")
            except ValueError as e:
                st.error(f"Erro na falsa posição: {e}")
                x_falsa = None

            try:
                x_newton, iter_newton = newton_raphson_9(x0)
                st.write(f"Newton-Raphson: x = {x_newton:.15f}, Iterações = {iter_newton}")
            except Exception as e:
                st.error(f"Erro no Newton-Raphson: {e}")
                x_newton = None

            # Verificação
            st.markdown("### 🔬 Verificação")
            if x_biss is not None:
                st.write(f"f(x_biss) = {f(x_biss):.2e}")
            if x_falsa is not None:
                st.write(f"f(x_falsa) = {f(x_falsa):.2e}")
            if x_newton is not None:
                st.write(f"f(x_newton) = {f(x_newton):.2e}")