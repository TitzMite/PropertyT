import numpy as np
import grp_CSA2 as A2
import Ozawa_SDP as OZ

radius = 1
identity = A2.A2Word.identity()

generators = A2.geometric_gens

#exclude the
factors = A2.geometric_ball(radius)

products = A2.geometric_ball(2 * radius)
product_index = {g: i for i, g in enumerate(products)}

num_generators = len(generators)

# Gram multiplication table, indexed directly into products.
mult = np.zeros((len(factors), len(factors)), dtype=int)

for i, g in enumerate(factors):
    for j, h in enumerate(factors):
        mult[i, j] = product_index[g.inverse() * h]

# Coefficients of delta, including its identity term.
delta = np.zeros(len(products), dtype=int)
delta[product_index[identity]] = num_generators

for g in generators:
    delta[product_index[g]] -= 1

minus_delta = -delta

# Compute delta^2 using the support of delta,
# independently of the Gram factors.
delta_support = [identity] + generators

delta_sq = np.zeros(len(products), dtype=int)

for g in delta_support:
    for h in delta_support:
        k = product_index[g * h]
        delta_sq[k] += (
            delta[product_index[g]]
            * delta[product_index[h]]
        )

def check_T():
    return OZ.prove_property_T(
        mult,
        delta_sq,
        minus_delta,
        radius,
        num_generators,
    )

result = check_T()