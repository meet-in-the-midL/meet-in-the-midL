from sage.all import pari, ZZ, PolynomialRing, GF
from sage.interfaces.gp import gp

gp.read("Isogeny/isogeny.gp")
from sage.libs.pari.convert_sage import gen_to_sage

def ssl_graph_edges(p : int, l : int) -> list[tuple[str, str, str]]:
    gpgraph = gp.ssl_graph(p, l)
    jinv = pari(gpgraph[1])
    edges : list[list[int]] = pari(gpgraph[2])
    # First, we sort out the labels
    labels : list[str] = []
    Fp2 = GF(p ** 2)
    # It's annoying to convert finite fields to Sage.
    # Make a polynomial ring in the integers
    F = PolynomialRing(ZZ, name = "y")
    y = F.gen(0)
    # Relies on ffgen always outputting the same generator. We also use y, the same variable as ssl_graph.
    gpgen = pari(p ** 2).ffgen(y)
    # Let's find the min poly for gpgen=y, which is the variable used in jinv.
    mpol = y ** 2 - eval(gen_to_sage((gpgen ** 2).Str()))
    # The polynomial
    mpol_Fp2 = mpol.change_ring(Fp2)
    # We have written it in Fp2=GF(p^2)
    gpgen_inFp2 = mpol_Fp2.any_root()
    # The j-invariants, as polynomials. The integral ones are OK, the others we must convert.
    sage_j = eval(gen_to_sage(jinv.Str()))
    notintclass = type(mpol)
    # Let's go through and convert them!
    for j in sage_j:
        # Integer, i.e. F_p j-invariant.
        if type(j) != notintclass:
            labels.append(str(j))
            continue
        # Not an integer, so must do more work.
        coefs = j.list()
        j_inv = coefs[0] + coefs[1] * gpgen_inFp2
        j_label = str(j_inv)
        labels.append(j_label)

    # Initialize empty graph, then loop over the edges, adding them in
    out_edges = []
    for i, vtx in enumerate(edges):
        for j, edge in enumerate(vtx):
            out = labels[i], labels[edge - 1], f'{l}'
            out_edges.append(out)
    return out_edges

