#encoding: utf-8
from __future__ import print_function
from .common import ProgressBarBase
import sys

class ProgressBarTerminal(ProgressBarBase):

    def __init__(self,
                 iterable_or_max,
                 desc='iterating', key=None, autohide=False, quiet=False,
                 format_str='%(title)s:|%(bar)s| %(currentp1)d/%(maxp1)d [%(elapsed).1fs < %(eta_avg).0f(%(eta_stddev).0f)s]',
                 width=60,framerate = 2):
        super(ProgressBarTerminal, self).__init__(iterable_or_max, desc, key, autohide, quiet)
        self.format_strs = format_str.split('%(bar)s')
        self.width = width
        self.phases = (u' ', u'▏', u'▎', u'▍', u'▌', u'▋', u'▊', u'▉', u'█')
        #self.phases = (u' ', u'█')
        self.last_percent = 0
        self.framerate = framerate

    def print_output(self):
        self.currentp1 = self.current
        self.maxp1 = self.max
        parts = [format % self for format in self.format_strs]
        parts[1:1] = self.bar(self.width - sum(map(len, parts)))
        print(''.join(parts),end='\r')
        sys.stdout.flush()

    def start(self):
        super(ProgressBarTerminal, self).start()
        self.last_print_time = self.start_time
        self.print_output()

    def advance(self):
        super(ProgressBarTerminal, self).advance()
        if self.percent != self.last_percent or (self.last_time - self.last_print_time) > 1./self.framerate:
            self.last_percent = self.percent
            self.print_output()
            self.last_print_time = self.last_time

    def finish(self):
        super(ProgressBarTerminal, self).finish()
        if not self.autohide:  
            self.print_output()
            print()
    
    def hide(self): 
        super(ProgressBarTerminal, self).hide()
        print((' ' * self.width)),
        sys.stdout.flush()
        
    def bar(self, bar_width):
        completely_filled = self.current * bar_width / self.max
        phase = (self.current * bar_width * len(self.phases) / self.max) % len(self.phases)

        return (self.phases[-1] * completely_filled +
                (self.phases[phase] if completely_filled < bar_width else '') +
                self.phases[0] * (bar_width - completely_filled))

