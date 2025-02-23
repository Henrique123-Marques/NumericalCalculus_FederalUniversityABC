'''
'''

import numpy as np

# Epsilon da máquina
eps = np.finfo(float).eps
print(f"Epsilon da máquina: {eps}")

# Função corrigida f(x) = -1/x³ - 1/x² + 2 (para E = -2)
def f(x):
    return -1/(x**3) - 1/(x**2) + 2

# Derivada de f(x): f'(x) = 3/x⁴ + 2/x³
def df(x):
    return 3/(x**4) + 2/(x**3)

# Critério de parada
def erro_relativo(x_novo, x_velho):
    return abs(x_novo - x_velho) <= eps * max(1, abs(x_novo))

# Método de Bisseção
def bissecao(a, b):
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

# Método da Falsa Posição
def falsa_posicao(a, b):
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

# Método de Newton-Raphson
def newton_raphson(x0):
    x_velho = x0
    iter_count = 0
    
    while True:
        x_novo = x_velho - f(x_velho) / df(x_velho)
        if erro_relativo(x_novo, x_velho):
            return x_novo, iter_count
        x_velho = x_novo
        iter_count += 1

# Intervalo inicial corrigido
a = -2.0
b = -0.5
print(f"f({a}) = {f(a)}, f({b}) = {f(b)}")  # Verifica sinais

# Chute inicial para Newton-Raphson
x0 = -0.75

# Inicializando variáveis como None
x_biss = None
x_falsa = None
x_newton = None

# Executando os métodos
try:
    x_biss, iter_biss = bissecao(a, b)
    print(f"Bisseção: x = {x_biss:.15f}, Iterações = {iter_biss}")
except ValueError as e:
    print(f"Erro na bisseção: {e}")

try:
    x_falsa, iter_falsa = falsa_posicao(a, b)
    print(f"Falsa Posição: x = {x_falsa:.15f}, Iterações = {iter_falsa}")
except ValueError as e:
    print(f"Erro na falsa posição: {e}")

try:
    x_newton, iter_newton = newton_raphson(x0)
    print(f"Newton-Raphson: x = {x_newton:.15f}, Iterações = {iter_newton}")
except Exception as e:
    print(f"Erro no Newton-Raphson: {e}")

# Verificação apenas para valores definidos
print("\nVerificação:")
if x_biss is not None:
    print(f"f(x_biss) = {f(x_biss):.2e}")
if x_falsa is not None:
    print(f"f(x_falsa) = {f(x_falsa):.2e}")
if x_newton is not None:
    print(f"f(x_newton) = {f(x_newton):.2e}")