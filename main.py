from random import choice as rchoice
from sage.all import Graphics, floor, ceil, scatter_plot

from mitm import *

def cumulative(B, vals) -> int:
    foo = [x for x in vals if x <= B]
    return len(foo)

def cdf_approx(
        vals : list[float],
        Bmin : float,
        Bmax : float,
        color : str
    ) -> Graphics:
    Bs = range(floor(Bmin), ceil(Bmax) + 1)
    pairs = [(B, cumulative(B, vals)) for B in Bs]
    plot = scatter_plot(
        pairs, markersize = 3, facecolor = color, edgecolor = "none"
    )
    return plot

NUM_SAMPLES = 1000

mitm_strategy = MitMStrategy(
    parts_iterator = PartsIterator.non_inc,
    stopping_condition = StoppingCondition.ISOG_CONJ | StoppingCondition.CYCLE,
    sides = Side.left
)

p = next_prime(10 ** 4)
vs = get_L_graph(p, [2]).vertices()
ells = [2, 3, 5]
samples = [rchoice(vs) for _ in range(NUM_SAMPLES)]
dicts : list[dict[vPair, float]] = []

for i in range(1, 4):
    for L in combinations(ells, r = i):
        L_list = list(L)
        cost_dict = get_MitM_cost_dict(p, L_list, mitm_strategy, samples)
        cost_function = lambda ell : get_isogeny_cost(ell, "modular")

        pair_cost_map : Callable[[vPair], float] = (
            lambda pair : cost_pair(cost_dict[pair], cost_function)
        )

        pair_cost_dict = {pair : pair_cost_map(pair) for pair in cost_dict}
        dicts.append(pair_cost_dict)

vals = [list(d.values()) for d in dicts]
maxB = int(max([max(vals) for vals in vals]))

plot = cdf_approx(vals[0], 0, int(maxB), "red")
save(plot, "plot.png")