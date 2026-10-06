from sage.all import EllipticCurve_from_j
from sage.graphs.digraph import DiGraph
from sage.modular.ssmod.ssmod import supersingular_j

from data import *
from ssl_pari import ssl_graph_edges

def get_L_graph(
        p : Prime,
        L : list[Prime],
        module : str = 'pari'
    ) -> DiGraph:
    edges = []
    for ell in L:
        try:
            data = fetch_data(p, ell, DataSets.GRAPH)
        except FileNotFoundError:
            print(f"computing G_{p}_{ell}")
            data = gen_ell_graph_edges(p, ell, module)
            store_data(p, ell, DataSets.GRAPH, data)
        edges = edges + data
    G = DiGraph(edges, loops = True, multiedges = True)
    return G

def gen_ell_graph_edges(
        p : Prime,
        ell : Prime,
        module : str
    ) -> list[Edge]:
    """
    Computes the L-isogeny graph over Fp2 for the given list of L's.
    - Each edge is of the form (start_vertex, end_vertex, "{degree},{i}) where i is a counter that is useful to distinguish multi-edges
    """
    if module == "sage":
        Fp2 = GF(p ** 2)
        j = supersingular_j(Fp2)
        E = EllipticCurve_from_j(j)
        G_ell = E.isogeny_ell_graph(ell, directed = True, label_by_j = True)
        tmp_edges = G_ell.edges()
    elif module == "pari":
        tmp_edges = ssl_graph_edges(p, ell)
    else:
        raise NotImplementedError
    seen : list[Edge] = []
    edges = []
    for edge in tmp_edges:
        count = seen.count(edge)
        edges.append((edge[0], edge[1], f"{ell},{count}"))
        seen.append(edge)
    return edges