from collections import Counter


class PMF():
    """Probability Mass Function
    Having bins with the length of n
    each bin is a counter dictionary, representing the number of
    occurence of a symbol at a position

    EX:
        >>> p = PMF(3)
        >>> p.increment(0, 'A')
        >>> p.get_distribution(0)
        {'A': 1.0}
        >>> p.increment(0, 'B')
        >>> p.get_distribution(0)
        {'A': 0.5, 'B': 0.5}

    """
    def __init__(self, length):
        self.bins = [Counter() for i in range(length)]

    def get_distribution(self, pos):
        return {k: v*(1.0/sum(self.bins[pos].values())) for k, v in self.bins[pos].items()}

    def increment(self, pos, symbol):
        self.bins[pos][symbol] += 1
