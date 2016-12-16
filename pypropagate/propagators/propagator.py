from ..solver import Solver

cached_evaluators = {}


class Propagator(Solver):
    def __init__(self, settings):
        super(Propagator, self).__init__(settings)

        coordinates = settings.partial_differential_equation.coordinates

        self._x, self._t = coordinates[0].symbol, coordinates[2].symbol
        if self.ndim > 1: self._y = coordinates[1].symbol

        self._nx, self._nt = self._get_as((coordinates[0].steps, coordinates[2].steps), int, settings)
        if self.ndim > 1: self._ny = self._get_as(coordinates[1].steps, int, settings)

        self._xmin, self._xmax, self._tmin, self._tmax = settings.get_numeric(
            (coordinates[0].min, coordinates[0].max, coordinates[2].min, coordinates[2].max))
        if self.ndim > 1: self._ymin, self._ymax = settings.get_numeric((coordinates[1].min, coordinates[1].max))

        import expresso.pycas as pc
        pe = settings.partial_differential_equation

        self._F_is_zero = settings.get_unitless(pe.F) == pc.Zero
        self._F_is_constant_in_z = settings.get_numeric(pc.derivative(pe.F, self._t)) == pc.Zero
        self._F_is_constant = self._F_is_constant_in_z and settings.get_numeric(pc.derivative(pe.F, self._x)) == pc.Zero
        if self.ndim > 1:
            self._F_is_constant &= settings.get_numeric(pc.derivative(pe.F, self._y)) == pc.Zero

        if self.ndim > 1:
            self.__coordinate_names = [c.index.name for c in coordinates]
        else:
            self.__coordinate_names = [c.index.name for c in (coordinates[0], coordinates[2])]

    def _set_initial_field(self, settings):
        sb = settings.simulation_box
        z = sb.coordinates[2]
        u0 = self._get_evaluators(settings.partial_differential_equation.u0.subs(z.symbol, z.min), settings)
        initial = u0(*self._get_indices())
        self.__initial = initial
        self.set_field(initial)

    def _reset(self):
        self.set_field(self.__initial)

    def _evaluate(self, expr, settings):
        import expresso.pycas as pc
        expr = pc.S(expr)

        if self.ndim == 1:
            sb = settings.simulation_box
            pde = settings.partial_differential_equation
            x, y, z = pde.coordinates
            y0 = getattr(pde, y.name + '0')
            expr = settings.get_optimized(expr.subs(y.symbol, y0))
            try:
                y0i = settings.get_as(y.step.subs(y.symbol, y0), int)
                expr = settings.get_optimized(expr.subs(y.index, y0i))
            except:
                pass
            return expr
        else:
            return settings.get_optimized(expr)

    def _get_as(self, expr, type, settings):
        return settings.get_as(self._evaluate(expr, settings), type)

    def _create_x_indices(self):
        import numpy as np
        return np.arange(self._nx, dtype=np.uint)

    def _get_x_indices(self):
        try:
            return self.__x_indices
        except AttributeError:
            self.__x_indices = self._create_x_indices()
            return self._get_x_indices()

    def _create_y_indices(self):
        import numpy as np
        return np.arange(self._ny, dtype=np.uint)

    def _get_y_indices(self):
        try:
            return self.__y_indices
        except AttributeError:
            self.__y_indices = self._create_y_indices()
            return self._get_y_indices()

    def __get_xy_indices(self):
        import numpy as np
        npy, npx = np.meshgrid(self._create_y_indices(), self._create_x_indices())
        return npx, npy

    def _get_indices(self):
        try:
            self.__z_indices.fill(self._i)
            return self.__indices + [self.__z_indices]
        except AttributeError:
            import numpy as np
            self.__indices = [self._create_x_indices()] if self.ndim == 1 else list(self.__get_xy_indices())
            self.__z_indices = np.zeros(self.__indices[0].shape, dtype=np.uint)
            return self._get_indices()

    def _get_transposed_indices(self):
        try:
            self.__transposed_z_indices.fill(self._i)
            return self.__transposed_indices + [self.__transposed_z_indices]
        except AttributeError:
            import numpy as np
            self.__transposed_indices = [self._create_x_indices()] if self.ndim == 1 else [np.ascontiguousarray(i.transpose()) for i in self.__get_xy_indices()]
            self.__transposed_z_indices = np.zeros(self.__transposed_indices[0].shape, dtype=np.uint)
            return self._get_transposed_indices()

    def _get_indices_dict(self):
        return {n: v for n, v in zip(self.__coordinate_names, self._get_indices())}

    def _get_initial(self):
        return self.__initial

    def _get_evaluators(self, expressions, settings, compile_to_c=None, args=None, **kwargs):
        import expresso.pycas as pc
        if not isinstance(expressions, (list, tuple)):
            return_single = True
            expressions = [expressions]
        else:
            return_single = False

        x, y, z = settings.simulation_box.coordinates

        if args is None:
            args = (x.index, y.index, z.index) if self.ndim == 2 else (x.index, z.index)
        expressions = [self._evaluate(expr, settings) for expr in expressions]

        if 'return_type' not in kwargs:
            kwargs['return_type'] = pc.Types.Complex

        def is_constant(expr):
            try:
                expr.N()
                return True
            except:
                return False

        definitions = [pc.FunctionDefinition('f%s' % i, args, expr, **kwargs)
                       for i, expr in enumerate(expressions) if not is_constant(expr)]
        if compile_to_c == None:
            compile_to_c = self.ndim > 1

        if not compile_to_c:
            lib = pc.ncompile(*definitions)
        else:
            key = tuple(expressions)
            if key in cached_evaluators:
                lib = cached_evaluators[key]
            else:
                lib = pc.ccompile(*definitions)
                if len(cached_evaluators) > 100:
                    cached_evaluators.clear()
                cached_evaluators[key] = lib

        def get_constant_expression(expr):
            try:
                c = float(expr.N(20))
            except:
                c = complex(expr.N(20))
            import numpy as np
            def constant_expression(*args, **kwargs):
                res = kwargs.pop('res', None)
                if res is not None:
                    res.fill(c)
                    return res
                return np.full(np.array(args[0]).shape, c, dtype=self.dtype)

            return constant_expression

        res = [getattr(lib, 'f%s' % i) if hasattr(lib, 'f%s' % i) else get_constant_expression(expressions[i]) for i in
               range(len(expressions))]

        if return_single:
            return res[0]
        return res

    def get_boundary_indices(self):
        if self.ndim == 1:
            return [0, self._nx-1]
        else:
            return [range(self._nx) * 2 + [0] * self._nx + [self._nx-1] * self._nx,
                    [0] * self._ny + [self._ny-1] * self._ny + range(self._ny) * 2]

