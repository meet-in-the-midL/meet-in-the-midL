import math
from sage.all import GF, is_prime, next_prime, prime_range

from custom_types import *

# ======================================================================== #
#
# MATHS UTILS
#
# ======================================================================== #

def check_pL_input(p: Prime, L: list[Prime]) -> bool:
    return is_prime(p) and all([is_prime(ell) for ell in L]) #and (p % 12 == 1)

def p_mod12_range(pmin, pmax) -> list[Prime]:
    primes = prime_range(pmin, pmax)
    return [p for p in primes if p % 12 == 1]

def eig_threshold(L: list[Prime]) -> float:
    d = sum(L) + len(L)
    threshold = 2 * math.sqrt(d - 1)
    return threshold

def gen_chain(L1 : list[Prime], L2 : list[Prime]) -> list[list[Prime]]:
    extras = [x for x in L2 if x not in L1]
    current = L1
    chain: list[list[Prime]] = []
    chain.append(current.copy())
    for x in extras:
        current = current + [x]
        chain.append(list(current))
    return chain

def pearson(x: list[float], y: list[float]) -> float | None:
    n = len(x)
    if n == 0:
        raise ValueError("len(x) == 0")
    mx = sum(x) / n
    my = sum(y) / n
    vx = sum((xi - mx) ** 2 for xi in x)
    vy = sum((yi - my) ** 2 for yi in y)
    if vx == 0 or vy == 0:
        return None
    cov = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    return float(cov / math.sqrt(vx * vy))

def next_good_prime(pmin):
    p = next_prime(pmin)
    while p % 12 != 1:
        p = next_prime(p)
    return p

def get_conj_map(p) -> Callable[[Vertex], Vertex]:
    F = GF(p ** 2) # heavily assuming this is deterministic
    def inv(v : Vertex) -> Vertex:
        return str(F(v).conjugate())
    return inv

def conj_edge(
        edge : Edge,
        conjugate : Callable[[Vertex], Vertex],
    ) -> Edge:
    return (conjugate(edge[0]), conjugate(edge[1]), edge[2])

def conj_path(
        path : Path,
        conjugate : Callable[[Vertex], Vertex],
    ) -> Path:
    """
    Given a path from E -> E', returns the conjugate path E^(p) -> E'^(p)
    """
    conj_edge_unary = lambda edge : conj_edge(edge, conjugate)
    path_c = list(map(conj_edge_unary, path))
    return path_c

def rep_path(path : Path) -> list[Vertex]:
    """
    Returns the sequence of vertices lying on a path
    """
    if path == []:
        return []
    return [e[0] for e in path] + [path[-1][1]]

def start(edge : Edge) -> Vertex:
    return edge[0]

def end(edge : Edge) -> Vertex:
    return edge[1]

def reverse_edge(edge : Edge) -> Edge:
    """
    Excessive as copied over from other project
    """
    j0, j1, label = edge
    return j1, j0, label

def reverse_path(path: Path) -> Path:
    return [reverse_edge(edge) for edge in path[::-1]]

def edge_degree(edge: Edge) -> int:
    ell, _ = edge[2].split(",")
    return int(ell)

def edge_weight(edge: Edge) -> float:
    # TODO: work out what this should be
    return float(edge_degree(edge))

def L_id(L) -> str:
    return "-".join(f"{x}" for x in L)