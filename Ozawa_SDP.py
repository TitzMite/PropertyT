import numpy as np
import cvxpy as cp
from fractions import Fraction
import math

def print_fraction(label, x):
    print(label, x)
    print("   ≈", float(x))

def setup_property_T_sdp(mult, delta_sq, minus_delta):
    m = mult.shape[0]
    n = minus_delta.shape[0]
    epsilon = cp.Variable(nonneg=True)
    target = delta_sq + epsilon * minus_delta
    P = cp.Variable((m, m), symmetric=True)
    positions = [[] for _ in range(n)]
    for i in range(m):
        for j in range(m):
            positions[mult[i, j]].append((i, j))
    constraints = [P >> 0]
    for k in range(n):
        constraints.append(cp.sum([P[i, j] for i, j in positions[k]]) == target[k])
    problem = cp.Problem(cp.Maximize(epsilon), constraints)
    return problem, P, epsilon, constraints

def setup_fixed_property_T_sdp(mult, delta_sq, minus_delta):
    m = mult.shape[0]
    n = minus_delta.shape[0]
    epsilon_param = cp.Parameter(nonneg=True)
    target = delta_sq + epsilon_param * minus_delta
    P = cp.Variable((m, m), symmetric=True)
    positions = [[] for _ in range(n)]
    for i in range(m):
        for j in range(m):
            positions[mult[i, j]].append((i, j))
    constraints = [P >> 0]
    for k in range(n):
        constraints.append(cp.sum([P[i, j] for i, j in positions[k]]) == target[k])
    problem = cp.Problem(cp.Minimize(0), constraints)
    return problem, P, epsilon_param

def maximize_epsilon(problem, epsilon):
    problem.solve(solver=cp.CLARABEL, verbose=True)
    if problem.status not in ["optimal", "optimal_inaccurate"]:
        print("Maximum epsilon could not be found.")
        return None
    if epsilon.value is None:
        print("No epsilon was returned.")
        return None
    epsilon_max = float(epsilon.value)
    print("\nNumerical maximal epsilon =", epsilon_max)
    return epsilon_max

def factor_gram_matrix(P_value, tol=1e-8):
    if P_value is None:
        print("No Gram matrix was returned.")
        return None
    P_num = np.array(P_value, dtype=float, copy=True)
    P_num = 0.5 * (P_num + P_num.T)
    eigenvalues, V = np.linalg.eigh(P_num)
    smallest = eigenvalues[0]
    print("Smallest eigenvalue:", smallest)
    if smallest < -tol:
        print("P is not positive semidefinite within the chosen tolerance.")
        return None
    small = np.abs(eigenvalues) <= tol
    if np.any(small):
        print(f"Eigenvalues adjusted: {np.sum(small)} eigenvalue(s) set to zero.")
        eigenvalues[small] = 0.0
    else:
        print("Eigenvalues were not adjusted.")
    eigenvalues = np.maximum(eigenvalues, 0.0)
    return np.diag(np.sqrt(eigenvalues)) @ V.T

def find_sos_factor(problem, P, tol=1e-8):
    problem.solve(solver=cp.CLARABEL, verbose=True)
    if problem.status not in ["optimal", "optimal_inaccurate"]:
        print("No solution found.")
        return None
    if P.value is None:
        print("No solution matrix P was returned.")
        return None
    print("Solution found.")
    return factor_gram_matrix(P.value, tol=tol)

def integer_augmentation_roots(Q, scale=10**12, correction_col=0):
    Q_int = np.rint(scale * Q).astype(np.int64)
    row_sums = np.sum(Q_int, axis=1)
    print("Rows requiring augmentation correction:", np.count_nonzero(row_sums))
    print("Largest augmentation correction:", np.max(np.abs(row_sums)))
    Q_int[:, correction_col] -= row_sums
    if not np.all(np.sum(Q_int, axis=1) == 0):
        raise RuntimeError("Augmentation correction failed.")
    print("All rows are now in the augmentation ideal.")
    return Q_int, scale

def sos_error_1norm(Q_int, mult, target, scale=10**12):
    n = len(target)
    scale_sq = scale**2
    Q_int = Q_int.astype(object)
    P_int = Q_int.T @ Q_int
    sos_int = np.zeros(n, dtype=object)
    m = Q_int.shape[1]
    for i in range(m):
        for j in range(m):
            sos_int[mult[i, j]] += P_int[i, j]
    # Compare exact rational coefficients—no rounding.
    error = Fraction(0)
    for k in range(n):
        target_coefficient = Fraction(target[k])
        sos_coefficient = Fraction(int(sos_int[k]), scale_sq)
        error += abs(target_coefficient - sos_coefficient)
    print_fraction("l1-error:", error)
    return error

def check_property_T(epsilon, error, radius, num_generators):
    epsilon = Fraction(epsilon)
    error = Fraction(error)
    C = 2 ** (2 * math.ceil(math.log2(radius)))
    omega = C * error
    epsilon_cert = epsilon - omega
    print("C =", C)
    print_fraction("l1-error =", error)
    print_fraction("correction C*error =", omega)
    print_fraction("certified epsilon =", epsilon_cert)
    if epsilon_cert <= 0:
        print("Property (T) is not certified by this computation.")
        return False, epsilon_cert, None
    print("Property (T) is certified.")
    kappa_sq = Fraction(2 * epsilon_cert, num_generators)
    print_fraction("Kazhdan constant lower bound squared =", kappa_sq)
    return True, epsilon_cert, kappa_sq

def certify_factor(Q, epsilon, mult, delta_sq, minus_delta, radius, num_generators, scale):
    Q_int, scale_used = integer_augmentation_roots(Q, scale=scale)
    target_exact = np.array(delta_sq, dtype=object) + epsilon * np.array(minus_delta, dtype=object)
    error = sos_error_1norm(Q_int, mult, target_exact, scale=scale_used)
    has_T, epsilon_cert, kappa_sq = check_property_T(epsilon, error, radius, num_generators)
    if not has_T:
        return None
    kappa = math.sqrt(float(kappa_sq))
    print("Kazhdan constant lower bound =", kappa)
    return {
        "epsilon": epsilon,
        "error": error,
        "epsilon_cert": epsilon_cert,
        "kappa_sq": kappa_sq,
        "kappa": kappa,
        "Q_int": Q_int,
        "scale": scale_used
    }

def prove_property_T(mult, delta_sq, minus_delta, radius, num_generators, scale=10**12, tol=1e-8):
    problem, P, epsilon_var, constraints = setup_property_T_sdp(mult, delta_sq, minus_delta)
    epsilon_max = maximize_epsilon(problem, epsilon_var)
    if epsilon_max is None:
        return None
    print("\nTrying to certify the Gram matrix from the maximizing SDP directly.")
    Q = factor_gram_matrix(P.value, tol=tol)
    if Q is not None:
        epsilon_first = Fraction.from_float(epsilon_max)
        print_fraction("trying maximizing epsilon =", epsilon_first)
        result = certify_factor(Q, epsilon_first, mult, delta_sq, minus_delta, radius, num_generators, scale)
        if result is not None:
            print("Certification succeeded using the initial maximizing Gram matrix.")
            result["epsilon_max"] = epsilon_max
            result["used_initial_gram"] = True
            return result
    print("\nDirect certification failed. Starting fixed-epsilon descent.")
    step = Fraction(1, 1000)
    epsilon = Fraction(math.floor(1000 * epsilon_max), 1000)
    if float(epsilon) >= epsilon_max:
        epsilon -= step
    fixed_problem, fixed_P, epsilon_param = setup_fixed_property_T_sdp(mult, delta_sq, minus_delta)
    while epsilon > 0:
        print()
        print_fraction("trying certification epsilon =", epsilon)
        epsilon_param.value = float(epsilon)
        Q = find_sos_factor(fixed_problem, fixed_P, tol=tol)
        if Q is None:
            epsilon -= step
            continue
        result = certify_factor(Q, epsilon, mult, delta_sq, minus_delta, radius, num_generators, scale)
        if result is not None:
            result["epsilon_max"] = epsilon_max
            result["used_initial_gram"] = False
            return result
        epsilon -= step
    print("Property (T) was not certified.")
    return None