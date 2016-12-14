from .propagator import Propagator
from ..hankel import *
import numpy as np


class FresnelPropagatorRS(Propagator):
    ndim = 1
    dtype = np.complex128

    def __init__(self, settings, thread_count=1):
        import expresso.pycas as pc
        import numpy as np

        pe = settings.partial_differential_equation
        x, y, z = settings.partial_differential_equation.coordinates

        xi = getattr(settings.simulation_box,x.name + 'i')
        self.__ximin = float(settings.get_unitless(xi.subs(x.symbol, 0)))
        self.__ximax = float(settings.get_unitless(xi.subs(x.symbol, x.max)))
        self.__sx = self.__ximax - self.__ximin

        super(FresnelPropagatorRS, self).__init__(settings)
        self._thread_count = thread_count

        self._set_initial_field(settings)

        R, = self._get_evaluators([pc.exp(pe.F * z.step)], settings, return_type=pc.Types.Complex)

        D = pc.numpyfy(self._evaluate(pc.exp(-pe.A * z.step * (pc.Symbol('kx') ** 2)), settings))

        if self._F_is_constant_in_z:
            self.__R_step = R(*self._get_indices())
        else:
            self.__R = R

        kx = hankel_freq(self._nx) * ((self._nx - 1) * 1. / self.__sx)
        self.__D_step = D(kx=kx, **self._get_indices_dict())

        self.__hankel_resample_matrix = hankel_resample_matrix(self._nx, (np.arange(self._nx) - self.__ximin) * (self._nx * 1. / self.__sx), cache_key=(self._nx, self.__ximin, self.__sx), xmax=self._nx)

    def _create_x_indices(self):
        return self.__ximin + hankel_samples(self._nx, xmax=self.__sx)

    def _step(self):
        if self._F_is_constant:
            self.__freq_data *= self.__D_step * self.__R_step
        elif self._F_is_constant_in_z:
            self.__freq_data = hankel(self.__data, xmax=self._nx)
            self.__freq_data *= self.__D_step
            self.__data = hankel(self.__freq_data, kmax=self._nx)
            self.__data *= self.__R_step
        else:
            self.__freq_data = hankel(self.__data, xmax=self._nx)
            self.__freq_data *= self.__D_step
            self.__data = hankel(self.__freq_data, kmax=self._nx)
            self.__data *= self.__R(*self._get_indices())

    def _get_field(self):
        import numpy as np

        if self._F_is_constant:
            self.__data = hankel(self.__freq_data, kmax=self._nx)
            return np.dot(self.__hankel_resample_matrix, self.__data)
        else:
            return np.dot(self.__hankel_resample_matrix, self.__data)

    def _set_field(self, field):
        import numpy as np
        self.__data = field.astype(np.complex128)
        self.__freq_data = hankel(self.__data, xmax=self._nx)
