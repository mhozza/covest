import math


class IntervalTree:
    LEFT = 0
    RIGHT = 1

    def __init__(self, start=0, end=0, default=0, resize_to_pow_2=True):
        self.start = start
        self.original_end = end
        self.end = IntervalTree._next_pow_2(self.original_end)\
            if resize_to_pow_2 else self.original_end
        self.children = [None, None]
        self.value = default
        self._add = 0

    @staticmethod
    def _next_pow_2(n):
        return int(pow(2, math.ceil(math.log(n, 2))))

    @property
    def left(self):
        return self.children[IntervalTree.LEFT]

    @property
    def right(self):
        return self.children[IntervalTree.RIGHT]

    def _get_or_create(self, child):
        if self.children[child] is None:
            mid = (self.start + self.end) / 2
            self.children[IntervalTree.LEFT] = IntervalTree(
                self.start, mid, self.value, resize_to_pow_2=False)
            self.children[IntervalTree.RIGHT] = IntervalTree(
                mid, self.end, self.value, resize_to_pow_2=False)
        return self.children[child]

    def add_to_interval(self, start, end, value=1):
        if start < self.start or end > self.end:
            return
        if self.start == start and self.end == end:
            if self.value is not None:
                self.value += value + self._add
                self._add = 0
            else:
                self._add += value
        else:
            if start < self._get_or_create(IntervalTree.LEFT).end:
                self.left.add_to_interval(start, self.left.end, value)
            if end >= self._get_or_create(IntervalTree.RIGHT).start:
                self.right.add_to_interval(self.right.start, end, value)
            self.value = None

    def get_one(self, pos):
        if self.value is not None:
            self.value += self._add
            self._add = 0
            return self.value
        else:
            if self._add > 0:
                self.left._add += self._add
                self.right._add += self._add
                self._add = 0
            if pos < (self.start + self.end) / 2:
                return self.left.get_one(pos)
            else:
                return self.right.get_one(pos)
