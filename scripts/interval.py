import intervals

class Partition(object):
    def __init__(self, L = None):
        self.intervals = []
        if L is not None:
            self.intervals.append(intervals.closedopen(0, L).to_atomic())

    def __repr__(self):
        return f"Intervals: {self.intervals}"

    def __str__(self):
        return f"{self.intervals}"

    def __len__(self):
        return len(self.intervals)

    def add_interval(self, low, high):
        I = intervals.closedopen(low, high).to_atomic()

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

