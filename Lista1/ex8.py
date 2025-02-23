'''Exercicio 8
Para o dióxido de carbono (CO2), a equação de estado de um gás tem os
seguintes coeficientes: a = 0, 401Pa m^3 e b = 42,7.10^−6 m^3
. Determine o volume ocupado por 1000 moléculas deste gás, à temperatura de 300K e pressão de 3, 5.10^7Pa,
 pelo método de bisecção, com uma tolerância de 10^−12. Constante de Boltzmann é k = 1, 3806503*10^−23 Joule K^−1.
[p + a*(N/V )^2]*(V − N*b) = k*N*T .
Faça também outros algoritmos para resolução usando os métodos de falsa-
posição e de Newton-Raphson.
'''
import numpy as np
import matplotlib.pyplot as plt

a = 0.401
b = 42.7e-6
N = 1000 #Numero de moleculas
P = 3.5e7 #pressao
T = 300 #temperatura em Kelvin
k = 1.3806503e-23 #Cte de Boltzman
tolerancia = 1e-12

#Funcao f(V)=0 baseada na equacao de estado
def f(V):
  termo1 = P + a * (N/V)**2
  termo2 = V - N * b
  return termo1 * termo2 - k * N * T

#Derivada de f(V) para Newton-Raphson
def df(V):
  termo1 = P + a * (N/V)**2
  termo2 = -2 * a * N**2 / V**3
  return termo1 + (V - N * b) * termo2

#Metodo da Bisseção
def bissecao(a_intervalo, b_intervalo, tolerancia):
  if f(a_intervalo) * f(b_intervalo) >= 0:
    raise ValueError("f(a) e f(b) devem ter sinais opostos")
  
  Va = a_intervalo
  Vb = b_intervalo
  iteracao_contador = 0

  while (Vb-Va) > tolerancia:
    Vm = (Va + Vb) / 2
    if f(Vm) == 0:
      return Vm, iteracao_contador
    elif f(Va) * f(Vm) < 0:
      Vb = Vm
    else:
      Va = Vm
    iteracao_contador += 1
  
  return (Va + Vb) / 2, iteracao_contador

#Metodo da Falsa Posicao
def falsa_posicao(a_intervalo, b_intervalo, tolerancia):
  if f(a_intervalo) * f(b_intervalo) >= 0:
    raise ValueError("f(a) e f(b) devem ter sinais opostos")

  Va = a_intervalo
  Vb = b_intervalo
  iteracao_contador = 0

  while abs(Vb-Va) > tolerancia:
    Vm = (Va * f(Vb) - Vb * f(Va)) / (f(Vb) - f(Va))
    if f(Vm) == 0:
      return Vm, iteracao_contador
    elif(Va) * f(Vm) < 0:
      Vb = Vm
    else:
      Va = Vm
    iteracao_contador += 1
  
  return (Va + Vb) / 2, iteracao_contador

#Metodo Newton-Raphson
def newton_raphson(V0, tolerancia):
  V = V0
  iteracao_contador = 0

  while True:
    f_V = f(V)
    df_V = df(V)
    V_novo = V - f_V / df_V

    if abs(V_novo - V) < tolerancia:
      return V_novo, iteracao_contador
    V = V_novo
    iteracao_contador += 1

#Intervalos definidos
Va = 0.0427 #proximo de N*b
Vb = 0.0428 #um pouco maior
V0 = 0.04275 # chute inicial para o metodo Newton-Raphson

#Executando os metodos
try:
  V_bissecao, iteracao_bissecao = bissecao(Va, Vb, tolerancia)
  print(f'Bissecao: V = {V_bissecao:.15f} m³, Iteracoes = {iteracao_bissecao}')
except ValueError as e:
  print(f'Erro na bissecao: {e}')

try:
  V_falsa, iteracao_falsa = falsa_posicao(Va, Vb, tolerancia)
  print(f'Falsa Posicao: V = {V_falsa:.15f} m³, Iteracoes = {iteracao_falsa}')
except ValueError as e:
  print(f'Erro na falsa posicao: {e}')

try:
  V_newton, iteracao_newton = newton_raphson(V0, tolerancia)
  print(f'Newton-Raphson: V = {V_newton:.15f} m³, Iteracoes = {iteracao_newton}')
except Exception as e:
  print(f'Erro no Newton-Raphson: {e}')

#Verificando os resultados
print("\nVerificação dos valores encontrados:")
print(f"f(V_bissecao) = {f(V_bissecao):.2e}")
print(f"f(V_falsa) = {f(V_falsa):.2e}")
print(f"f(V_newton) = {f(V_newton):.2e}")

# Novo intervalo para o gráfico com zoom
# Usaremos uma janela de ±0.000005 ao redor da média das soluções
V_media = (V_bissecao + V_falsa + V_newton) / 3
delta_V = 0.00000000001
V_range_zoom = np.linspace(V_media - delta_V, V_media + delta_V, 1000)
f_values_zoom = [f(V) for V in V_range_zoom]

# Criando o gráfico com zoom
plt.figure(figsize=(10, 6))
plt.plot(V_range_zoom, f_values_zoom, 'b-', label='f(V)')
plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
plt.plot(V_bissecao, f(V_bissecao), 'ro', label=f'Bisseção (V={V_bissecao:.10f})')
plt.plot(V_falsa, f(V_falsa), 'go', label=f'Falsa Posição (V={V_falsa:.10f})')
plt.plot(V_newton, f(V_newton), 'yo', label=f'Newton-Raphson (V={V_newton:.10f})')
plt.grid(True)
plt.xlabel('Volume (m³)')
plt.ylabel('f(V)')
plt.title('Equação de Estado do CO₂ - Zoom nos Resultados')
plt.legend()
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))  # Formato científico para eixo x
plt.show()

# Exibindo a diferença entre os métodos
print("\nDiferenças entre os métodos:")
print(f"Bisseção - Falsa Posição: {(V_bissecao - V_falsa):.15e} m³")
print(f"Bisseção - Newton-Raphson: {(V_bissecao - V_newton):.15e} m³")
print(f"Falsa Posição - Newton-Raphson: {(V_falsa - V_newton):.15e} m³")