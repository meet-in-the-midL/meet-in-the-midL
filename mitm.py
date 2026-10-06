from dataclasses import dataclass
from itertools import combinations
import statistics

from graphs import *
from strategies import *

class Side(Enum):
    both = auto()
    left = auto()
    right = auto()

class StoppingCondition(IntFlag):
    none = 0
    MEET = auto()
    CYCLE = auto()
    MEET_CONJ = auto()
    ISOG_CONJ = auto()

@dataclass
class MitMStrategy:
    r"""
    3 arguments may be passed that determine how a MitM style search is performed.

    - `parts_iterator` : B-smooth degree isogenies decompose in to isogenies with degrees that form a sequence of primes. These sequences form a tree structure, and this argument determines how this tree is explored
        **Options:**
        - 'non_increasing' (default)
        - 'non_decreasing'
    - `stopping_condition` : Standard MitM forms trees expanding from the starting vertex(es) until a collision with the other tree is found. In the context of isogenies, there are other phenomena of interest that can be found.
        **Options:**
        - 'none' (default)
        - 'loop' : Halts if a two different paths to the same vertex are found
        - 'conjugacy' : Halts if a collision is found up to Galois conjugacy
        - 'conjugacy_loop' : Both of the above two conditions
    - ``sides` : Search is typically conducted by expanding from both sides, but optionally just one can be chosen
        **Options:**
        - 'both' (default)
        - 'left'
        - 'right'
    """
    parts_iterator : PartsIterator = PartsIterator.non_inc
    stopping_condition : StoppingCondition = StoppingCondition.MEET
    sides : Side = Side.both

    def __repr__(self) -> str:
        return ",".join([
            str(self.parts_iterator.value),
            str(self.stopping_condition.value),
            str(self.sides.value)
        ])

def weighted_Mitm(
        G : DiGraph,
        s : Vertex,
        t : Optional[Vertex],
        strategy : MitMStrategy,
        conjugate : Optional[Callable[[Vertex], Vertex]] = None,
        path : bool = False
    ) -> tuple[Optional[Path], dict[Prime, int]]:

    if s == t:
        if path:
            return [], {}
        else:
            return None, {}

    # ------------------------------------------------------------------ #
    # Helpers 
    # ------------------------------------------------------------------ #  

    def outgoing_edges_iterator_l(
            v : Vertex,
            l : Prime,
        ) -> Iterable[Edge]:

        edge_gen: Iterable[Edge] = G.outgoing_edge_iterator(v)
        for e in edge_gen:
            if edge_degree(e) == l:
                yield e

    prev_f : dict[Vertex, list[Edge]] = {s : []}
    reached_f : set[Vertex] = set([s])
    pq_f : list[tuple[PrimePartition, Vertex]] = [
        (PrimePartition.root(), s)
    ]
    dbl_reached : set[Vertex] = set()

    prev_b : dict[Vertex, list[Edge]] = {}
    reached_b : set[Vertex] = set()
    pq_b : list[tuple[PrimePartition, Vertex]] = []

    if t is not None:
        prev_b : dict[Vertex, list[Edge]] = {t : []}
        reached_b : set[Vertex] = set([t])
        dbl_reached : set[Vertex] = set()
        pq_b : list[tuple[PrimePartition, Vertex]] = [
            (PrimePartition.root(), t)
        ]

    edge_comp : dict[Prime, int] = {}
    stopping = strategy.stopping_condition

    # -------------------------------------------------------------------- #
    # Main iteration
    # -------------------------------------------------------------------- #

    def expand_one(
            g : PrimePartition,
            u : Vertex,
            pq : list[tuple[PrimePartition, Vertex]],
            prev : dict[Vertex, list[Edge]],
            reached : set[Vertex],
            other_reached : set[Vertex],
            dbl_reached : set[Vertex],
            stopping : StoppingCondition,
        ) -> StoppingCondition:

        ell = g[-1]
        prevs = prev[u]
        cond = StoppingCondition.none

        for edge in outgoing_edges_iterator_l(u, ell):
            rev = reverse_edge(edge)
            if rev in prevs:
                continue

            inc = edge_comp.get(ell, 0) + 1
            edge_comp[ell] = inc
            v = end(edge)

            if v in reached:
                # indicates a small cycle
                cond |= StoppingCondition.CYCLE
                prev[v].insert(0, edge)
                dbl_reached.add(v)

            reached.add(v)
            prev[v] = prev.get(v, []) + [edge]

            if v in other_reached:
                cond |= StoppingCondition.MEET

            if not conjugate is None:
                vc = conjugate(v)
                if vc in reached:
                    # either because its an Fp curve we've just added, or otherwise
                    # we have found a separable path from v to vc via s or t
                    cond |= StoppingCondition.ISOG_CONJ
                if vc in other_reached:
                    cond |= StoppingCondition.MEET_CONJ

            if cond & stopping:
                break
            
            heapq.heappush(pq, (g, v))
        return cond & stopping

    # -------------------------------------------------------------------- #
    # Main loop
    # -------------------------------------------------------------------- #
    
    start_b = 0 if strategy.sides != Side.right else 1
    b_bound = 2 if strategy.sides != Side.left else 1
    gen = PrimePartitionIterator(strategy.parts_iterator)

    """
    The loop invariant(s) depend on the strategy chosen. The loop exits when the stopping condition of interest is met.
    """

    b = start_b
    collision = 0
    while not collision:
        b = start_b
        g = next(gen)
        sub_pq_f = [e for e in pq_f if e[0] == g.parent()]
        sub_pq_b = [e for e in pq_b if e[0] == g.parent()]
        while not collision and b < b_bound:
            try:
                if b == 0:
                    _, u = sub_pq_f.pop()
                    collision = expand_one(
                        g, u, pq_f, prev_f, reached_f, reached_b, dbl_reached, stopping
                    )
                else:
                    _, u = sub_pq_b.pop()
                    collision = expand_one(
                        g, u, pq_b, prev_b, reached_b, reached_f, dbl_reached, stopping
                    )
            except IndexError:
                b += 1
                continue

    # -------------------------------------------------------------------- #
    # Reconstruct path
    # -------------------------------------------------------------------- #

    if not path:
        return None, edge_comp

    prev1 : Callable[[Vertex], Optional[Edge]]
    prev2 : Callable[[Vertex], Optional[Edge]]
    edge1 : Optional[Edge]
    edge2 : Optional[Edge]
    loop : Optional[Edge]
    path_f : Path = []
    path_b : Path = []

    if collision & StoppingCondition.MEET:
        meeting = reached_f.intersection(reached_b).pop()
        v1 = v2 = meeting
        prev1 = lambda v : next(iter(prev_f.get(v, [])), None)
        prev2 = lambda v : next(iter(prev_b.get(v, [])), None)

        edge1 = prev1(v1)
        edge2 = prev2(v2)
        path_f = []
        path_b = []
    elif collision & StoppingCondition.CYCLE:
        meeting = dbl_reached.pop()
        if b == 0:
            prev1 = lambda v : next(iter(prev_f.get(v, [])), None)
            prev2 = lambda v : next(iter(prev_f.get(v, [])), None)
            es = iter(prev_f.get(meeting, []))
        else:
            prev1 = lambda v : next(iter(prev_b.get(v, [])), None)
            prev2 = lambda v : next(iter(prev_b.get(v, [])), None)
            es = iter(prev_b.get(meeting, []))

        if collision == 2:
            edge1 = next(es, None)
            edge2 = next(es, None)
        else:
            loop = next(es, None)
            edge1 = edge2 = next(es, None)
            path_f = [loop] if not loop is None else []
    elif collision & StoppingCondition.ISOG_CONJ:
        assert not conjugate is None
        if b == 0:
            conj_set = map(conjugate, reached_f)
            meeting = reached_f.intersection(conj_set).pop()
            prev1 = lambda v : next(iter(prev_f.get(v, [])), None)
            prev2 = lambda v : next(iter(prev_f.get(v, [])), None)
        else:
            conj_set = map(conjugate, reached_b)
            meeting = reached_b.intersection(conj_set).pop()
            prev1 = lambda v : next(iter(prev_b.get(v, [])), None)
            prev2 = lambda v : next(iter(prev_b.get(v, [])), None)
        meeting_c = conjugate(meeting)
        assert not meeting_c is None
        edge1 = prev1(meeting)
        edge2 = prev2(meeting_c)
        path_f = []
        path_b = []
    elif collision & StoppingCondition.MEET_CONJ:
        assert not conjugate is None
        conj_set = map(conjugate, reached_b)
        meeting = reached_f.intersection(conj_set).pop()
        meeting_c = conjugate(meeting)
        assert not meeting_c is None
        prev1 = lambda v : next(iter(prev_f.get(v, [])), None)
        prev2 = lambda v : next(iter(prev_b.get(v, [])), None)
        edge1 = prev1(meeting)
        edge2 = prev2(meeting_c)
        frob : Edge = (meeting, meeting_c, "p,1")
        path_f = [frob]
        path_b = []
    else:
        raise ValueError
        
    while edge1 is not None:
        meeting_s = start(edge1)
        path_f = [edge1] + path_f
        edge1 = prev1(meeting_s)

    while edge2 is not None:
        meeting_t = start(edge2)
        path_b = [edge2] + path_b
        edge2 = prev2(meeting_t)

    if collision & StoppingCondition.ISOG_CONJ:
        assert not conjugate is None
        out = path_f + reverse_path(conj_path(path_b, conjugate))
    else:
        out = path_f + reverse_path(path_b)
    return out, edge_comp

def gen_MitM_cost_dict(
        p : Prime,
        L : list[Prime],
        strategy : MitMStrategy,
        samples : Optional[list[Vertex]] = None
    ) -> dict[vPair, dict[Prime, int]]:
    G = get_L_graph(p, L)
    vs : list[Vertex] = G.vertices()
    pairs = combinations(vs, 2)
    conj = get_conj_map(p)
    cost_dict = lambda pair : weighted_Mitm(
        G, pair[0], pair[1], strategy, conjugate = conj, path = False
    )[1]

    out : dict[vPair, dict[Prime, int]] = {}
    if strategy.sides == Side.both:
        n_pairs = binomial((p / 12), 2)
        for i, pair in enumerate(pairs):
            print(f"progress : {round(100 * i / n_pairs, 2)}%", end="\r")
            out[pair] = cost_dict(pair)
    elif strategy.sides == Side.left:
        n_vs = len(vs)
        dummy = vs[0]
        if not samples is None:
            for i, v in enumerate(samples):
                pair = v, dummy
                print(f"progress : {round(100 * i / len(samples), 2)}%", end="\r")
                out[pair] = cost_dict(pair)
        else:
            for i, v in enumerate(vs):
                pair = v, dummy
                print(f"progress : {round(100 * i / n_vs, 2)}%", end="\r")
                out[pair] = cost_dict(pair)
    elif strategy.sides == Side.right:
        n_vs = len(vs)
        dummy = vs[0]
        for i, v in enumerate(vs):
            pair = dummy, v
            print(f"progress : {round(100 * i / n_vs, 2)}%", end="\r")
            out[pair] = cost_dict(pair)
    else:
        raise ValueError
    return out

def get_MitM_cost_dict(
        p : Prime,
        L : list[Prime],
        strategy : MitMStrategy,
        samples : Optional[list[Vertex]]
    ) -> dict[vPair, dict[Prime, int]]:
    if not samples is None:
        data = gen_MitM_cost_dict(p, L, strategy, samples)
    else:
        try:
            data = fetch_data(p, L, DataSets.MITM, strategy)
        except FileNotFoundError:
            data = gen_MitM_cost_dict(p, L, strategy)
            store_data(p, L, DataSets.MITM, data, strategy)
    return data

def get_isogeny_cost(
        ell : Prime,
        formula : str = "modular",
    ) -> float:
    match formula:
        case "modular":
            return float(ell)
        case "sqrt":
            return float(ell) ** (1 / 2)
        case _:
            return float(ell) ** (3 / 2)

def cost_pair(
        steps_dict : dict[Prime, int],
        cost : Optional[Callable[[Prime], float]],
    ) -> float:
    if cost is None:
        cost = lambda ell : get_isogeny_cost(ell)
    total_cost = sum(
        [n * cost(ell) for ell, n in steps_dict.items()]
    )
    return total_cost

def MitM_cost(
        p : Prime, L : list[Prime],
        strategy : MitMStrategy,
        samples : Optional[list[Vertex]],
        cost : Optional[Callable[[Prime], float]] = None,
    ) -> dict[vPair, float]:
    """
    Returns the total computation involved in a MitM attack for the fixed pair, with the associated edge cost function
    """
    cost_dict = get_MitM_cost_dict(p, L, strategy, samples)
    if cost is None:
        cost = lambda ell : get_isogeny_cost(ell)

    total_cost : Callable[[vPair], float] = (
        lambda pair : cost_pair(cost_dict[pair], cost)
    )
    out = {pair : total_cost(pair) for pair in cost_dict}
    return out

def MitM_cost_summary(
        p : Prime,
        L : list[Prime],
        strategy : MitMStrategy,
        samples : Optional[list[Vertex]],
        cost_f : Optional[Callable[[Prime], float]] = None,
    ) -> tuple[float, float]:
    cost = MitM_cost(p, L, strategy, samples, cost_f).values()
    avg_comp = float(statistics.mean(cost))
    stdv = float(statistics.stdev(cost))
    data = avg_comp, stdv
    return data