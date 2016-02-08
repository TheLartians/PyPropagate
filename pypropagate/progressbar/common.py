from time import time
from math import sqrt


class ProgressBarBase(object):

    def __init__(self, iterable_or_max, title=None, key=None, autohide=False, quiet=False):
        try:
            # if any iterable given - use it
            self.iterable = iter(iterable_or_max)
            try:
                # if has len - it's the max
                self.max = len(iterable_or_max)
            except TypeError:
                # iterable without len - indeterminate
                self.max = None
        except TypeError:
            # not iterable - assume max is given
            self.iterable = xrange(iterable_or_max)
            self.max = iterable_or_max

        self.title = title
        self.key = key
        self.autohide = autohide
        self.quiet = quiet
        
        self.current = None
        self.square_time = 0
        
    def start(self):
        self.start_time = time()
        self.last_time = time()
        self.current = 0

    def advance(self):
        current_time = time()
        diff = (current_time - self.last_time)
        self.square_time += diff * diff
        self.last_time = time()
        self.current += 1

    def finish(self):
        if self.autohide:
            self.hide()

    def __iter__(self):
        self.start()
        for elem in self.iterable:
            yield elem
            self.advance()
        self.finish()

    def set_extra_text(self, text):
        self.extra_text = text

    def hide(self):
        pass

    @property
    def frac(self):
        return self.current / float(self.max)

    @property
    def percent(self):
        return int(self.frac * 100)

    @property
    def elapsed(self):
        return time() - self.start_time

    @property
    def avg_time(self):
        if self.current == 0: return 0
        return (self.last_time - self.start_time) / self.current
    
    @property
    def std_time(self):
        if self.current == 0: return 0
        return sqrt(max(0,self.square_time/self.current - self.avg_time**2 ))
        
    @property
    def eta_avg(self):
        return (self.max - self.current) * self.avg_time

    @property
    def eta_stddev(self):
        return (self.max - self.current) * self.std_time

    def __getitem__(self, key):
        # for % formatting
        try:
            return getattr(self, key)
        except:
            return float('NaN')
