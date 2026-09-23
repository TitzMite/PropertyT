import grp_CSA2 as A2
import sympy as sp

radius = 1
identity = A2.A2Word.identity()

# Include the identity in the Gram factors.
factors = A2.geometric_ball(radius)
products = A2.geometric_ball(2 * radius)

generators = A2.geometric_gens
num_generators = len(generators)
product_index = {g: i for i, g in enumerate(products)}

# Gram multiplication table, indexed into products.
mult = sp.zeros(len(factors), len(factors))

for i, g in enumerate(factors):
    for j, h in enumerate(factors):
        mult[i, j] = product_index[g.inverse() * h]

# delta, represented in the products basis.
delta = sp.zeros(len(products), 1)
delta[0] = num_generators

for g in generators:
    delta[product_index[g]] -= 1

minus_delta = -delta

# delta^2
delta_support = [identity] + generators
delta_sq = sp.zeros(len(products), 1)

for g in delta_support:
    for h in delta_support:
        k = product_index[g * h]
        delta_sq[k] += (
            delta[product_index[g]]
            * delta[product_index[h]]
        )

d = 4 + 6*sp.sqrt(2)
c = -6 - 2*sp.sqrt(2)

# Diagonal contribution
P1 = sp.Matrix([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, d, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, d, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, d, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, d, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d],
])

# Common contribution of 1
P2 = sp.Matrix([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
])

# Incidence contribution
P3 = sp.Matrix([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, c, c, 0, c, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, c, 0, 0, 0, c, c],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c, c, 0, c],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, c, c, 0, c, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, c, 0, c, 0, 0, 0, c],
    [0, 0, 0, 0, 0, 0, 0, 0, c, 0, 0, 0, c, c, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c, c, 0, c, 0],
    [0, c, 0, 0, 0, c, c, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, c, c, 0, c, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, c, c, 0, c, 0, 0, 0, 0, 0, 0, 0],
    [0, c, 0, c, 0, 0, 0, c, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, c, c, 0, c, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, c, 0, 0, 0, c, c, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, c, c, 0, c, 0, 0, 0, 0, 0, 0, 0, 0, 0],
])

P = P1 + P2 + P3
epsilon = 10 - 6*sp.sqrt(2)


def check_Gram_matrix(M, P, epsilon):
    target = delta_sq + epsilon * minus_delta
    sos = sp.zeros(len(products), 1)

    for i in range(len(factors)):
        for j in range(len(factors)):
            k = int(M[i, j])
            sos[k] += P[i, j]

    return (sos - target).applyfunc(sp.simplify)


residual = check_Gram_matrix(mult, P, epsilon)
identity_holds = residual == sp.zeros(len(products), 1)

# Verify positivity of the explicitly symmetric matrix P.
a = 8 + 12*sp.sqrt(2)
b = 22 + 12*sp.sqrt(2)
I = sp.eye(15)

R = (P * (P - a*I) * (P - b*I)).applyfunc(sp.expand)
is_psd = R == sp.zeros(15, 15)

print("Exact solution:", identity_holds)
print("Positive semidefinite:", is_psd)