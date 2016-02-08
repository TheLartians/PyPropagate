#encoding: utf-8
from .common import ProgressBarBase
import sys


class ProgressBarTerminal(ProgressBarBase):

    def __init__(self,
                 iterable_or_max,
                 title='Progress', key=None, autohide=False, quiet=False,
                 format_str='%(title)s %(currentp1)d/%(maxp1)d: %(percent)3d%% [%(bar)s] [%(elapsed).1f s] [eta %(eta_avg).0fs +- %(eta_stddev).0fs]',
                 width=90,framerate = 1):
        super(ProgressBarTerminal, self).__init__(iterable_or_max, title, key, autohide, quiet)
        self.format_strs = format_str.split('%(bar)s')
        self.width = width
        self.phases = (u' ', u'▏', u'▎', u'▍', u'▌', u'▋', u'▊', u'▉', u'█')
        self.last_percent = 0
        self.framerate = framerate
        #self.phases = (' ','=')
        
    def print_output(self):
        self.currentp1 = self.current+1
        self.maxp1 = self.max+1
        parts = [format % self for format in self.format_strs]
        parts[1:1] = self.bar(self.width - sum(map(len, parts)))
        print '\r' + ''.join(parts),
        sys.stdout.flush()

    def start(self):
        super(ProgressBarTerminal, self).start()
        print
        self.last_print_time = self.start_time
        self.print_output()

    def advance(self):
        super(ProgressBarTerminal, self).advance()
        if self.percent != self.last_percent or (self.last_time - self.last_print_time) > 1./self.framerate:
            self.last_percent = self.percent
            self.print_output()
            sys.stdout.flush()
            self.last_print_time = self.last_time

    def finish(self):
        super(ProgressBarTerminal, self).finish()
        if not self.autohide:
            print

    def hide(self):
        super(ProgressBarTerminal, self).hide()
        print (' ' * self.width) + '\r',
        sys.stdout.flush()

    def bar(self, bar_width):
        completely_filled = self.current * bar_width / self.max
        phase = (self.current * bar_width * len(self.phases) / self.max) % len(self.phases)

        return (self.phases[-1] * completely_filled +
                (self.phases[phase] if completely_filled < bar_width else '') +
                self.phases[0] * (bar_width - completely_filled))
