from enum import Enum, auto
import heapq
from itertools import count
from sage.all import next_prime

from custom_types import *

class PartsIterator(Enum):
    non_inc = auto()
    non_dec = auto()

class PrimePartition(list[int]):
    """

    """

    @classmethod
    def root(cls) -> "PrimePartition":
        """The root of the tree: the empty partition."""
        return cls()

    def is_root(self) -> bool:
        return len(self) == 0

    def parent(self) -> Optional["PrimePartition"]:
        """The partition this one extends, or None if this is the root."""
        if self.is_root():
            return None
        return PrimePartition(self[:-1])

    def first_child(self) -> "PrimePartition":
        return PrimePartition([*self, 2])

    def next_sibling(self) -> Optional["PrimePartition"]:
        """
        """
        if self.is_root():
            return None
        next_p = next_prime(self[-1])
        return PrimePartition([*self[:-1], next_p])

    def children(self):
        """
        Generator yielding all children of this node, in sort order.
        """
        child = self.first_child()
        while child is not None:
            yield child
            child = child.next_sibling()

    def prime_increase(self) -> bool:
        parent = self.parent()
        if parent is None:
            return True
        if parent.is_root():
            return True
        return self[-1] >= parent[-1]

    def prime_decrease(self) -> bool:
        parent = self.parent()
        if parent is None:
            return True
        if parent.is_root():
            return True
        return self[-1] <= parent[-1]

    def __repr__(self) -> str:
        return f"{list(self)}"
    
    def sort_key(self):
        return sum(self), - len(self), [x for x in self]

class PrimePartitionIterator:
    """
    Iterates over all PrimePartitions (except the root) in the order:
    sum ascending; ties broken by fewest parts; remaining ties broken by
    descending lexicographic order.

        [2] -> [3] -> [2,2] -> [5] -> [3,2] -> [3,3] -> [2,2,2] -> ...
    """

    def __init__(self, strategy : PartsIterator):
        self._heap : list[tuple[Any, Any, Any, Any, PrimePartition]] = []
        self._counter = count()
        self._push(PrimePartition.root().first_child())
        if strategy == PartsIterator.non_inc:
            self.nextf = self.next2
        elif strategy == PartsIterator.non_dec:
            self.nextf = self.next1
        else:
            raise ValueError


    def __iter__(self) -> "PrimePartitionIterator":
        return self

    def next1(self) -> PrimePartition:
        *_, node = heapq.heappop(self._heap)
        children = node.children()
        child = next(children)
        while not child.prime_increase():
            child = next(children)
        self._push(child)
        sibling = node.next_sibling()
        if sibling is not None:
            self._push(sibling)
        return node

    def next2(self) -> PrimePartition:
        *_, node = heapq.heappop(self._heap)
        child = node.first_child()
        self._push(child)
        sibling = node.next_sibling()
        if sibling is not None:
            if sibling.prime_decrease():
                self._push(sibling)
        return node
    
    def __next__(self) -> PrimePartition:
        return self.nextf()

    def _push(self, node: PrimePartition) -> None:
        heapq.heappush(
            self._heap, (*node.sort_key(), next(self._counter), node)
        )
