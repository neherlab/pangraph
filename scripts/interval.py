import intervals

class Partition(object):
    def __init__(self, L = None):
        self.intervals = []
        self.mapped = intervals.empty()
        if L is not None:
            self.intervals.append(intervals.closedopen(0, L).to_atomic())

    def __repr__(self):
        return f"Intervals: {self.intervals}"

    def __str__(self):
        return f"{self.intervals}"

    def __len__(self):
        return len(self.intervals)

    def __eq__(self, other):
        return len(self.intervals) == len(other.intervals) and \
               all(Is == Io for (Is, Io) in zip(self.intervals, other.intervals))

    def add_interval(self, low, high):
        I = intervals.closedopen(low, high).to_atomic()
        if self.mapped is None:
            self.mapped = self.mapped | I

        new_intervals = []
        for I0 in self.intervals:
            if I.is_empty() or not I0.overlaps(I):
                new_intervals.append(I0)
                continue

            Id = I0 - I
            if isinstance(Id, intervals.AtomicInterval):
                new_intervals.append(Id)
            else:
                new_intervals.extend([i for i in list(Id) if not i.is_empty()])
            new_intervals.append(I0 & I)
            I = I - I0

            if not isinstance(I, intervals.AtomicInterval):
                assert len(I) <= 1
                I = I.to_atomic()

        self.intervals = sorted(new_intervals)

def unittest():
    P1 = Partition(100)
    P1.add_interval(0, 10)
    P1.add_interval(90, 100)

    P2 = Partition(100)
    P2.add_interval(90, 100)
    P2.add_interval(0, 10)

    assert P1 == P2

    print("Passed test 1")

    P1.add_interval(8, 27)
    P1.add_interval(24, 63)

    P3 = Partition(100)
    P3.add_interval(0, 8)
    P3.add_interval(8, 10)
    P3.add_interval(10, 24)
    P3.add_interval(24, 27)
    P3.add_interval(27, 63)
    P3.add_interval(63, 90)
    P3.add_interval(90, 100)

    assert P1 == P3

    print("Passed test 2")


if __name__ == "__main__":
    unittest()
