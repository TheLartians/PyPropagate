
from ..solver import Solver

class Propagator(Solver):
    
    def __init__(self,settings):
        super(Propagator,self).__init__(settings)

        sb = settings.simulation_box

        self._x = sb.x
        if self.ndim > 1: self._y = sb.y
        self._t = sb.z

        self._nx,self._nt = settings.get_as((sb.Nx,sb.Nz),int)
        if self.ndim > 1: self._ny = settings.get_as(sb.Ny,int)
        self._xmin,self._xmax = settings.get_numeric((sb.xmin,sb.xmax))
        if self.ndim > 1: self._ymin,self._ymax = settings.get_numeric((sb.ymin,sb.ymax))
        self._tmin,self._tmax = settings.get_numeric((sb.zmin,sb.zmax))

        import expresso.pycas as pc
        pe = settings.partial_differential_equation

        self._F_is_zero = settings.get_unitless( pe.F ) == pc.Zero
        self._F_is_constant_in_z = settings.get_numeric(pc.derivative(pe.F, sb.z)) == pc.Zero
        self._F_is_constant = self._F_is_constant_in_z and settings.get_numeric(pc.derivative(pe.F, sb.x)) == pc.Zero
	if self.ndim > 1: self._F_is_constant &= settings.get_numeric(pc.derivative(pe.F, sb.y)) == pc.Zero

    def _set_initial_field(self,settings):
        sb = settings.simulation_box
        u0 = self._get_evaluators(settings.partial_differential_equation.u0.subs(sb.z,sb.zmin),settings)
        initial = u0(*self._get_coordinates())
        self.__initial = initial
        self.set_field(initial)

    def _reset(self):
        self.set_field(self.__initial)

    def __get_x_coordinates(self):
        import numpy as np
        return np.arange(self._nx,dtype=np.uint)

    def _get_x_coordinates(self):
        try:
            return self.__x_coordinates
        except AttributeError:
            self.__x_coordinates = self.__get_x_coordinates()
            return self._get_x_coordinates()

    def __get_y_coordinates(self):
        import numpy as np
        return np.arange(self._ny,dtype=np.uint)

    def _get_y_coordinates(self):
        try:
            return self.__y_coordinates
        except AttributeError:
            self.__y_coordinates = self.__get_y_coordinates()
            return self._get_y_coordinates()

    def __get_xy_coordinates(self):
        import numpy as np
        npy,npx = np.meshgrid(self.__get_y_coordinates(),self.__get_x_coordinates())
        return npx,npy

    def _get_coordinates(self):
        try:
            self.__z_coordinates.fill(self._i)
            return self.__coordinates + [self.__z_coordinates]
        except AttributeError:
            import numpy as np
            self.__coordinates = [self.__get_x_coordinates()] if self.ndim == 1 else list(self.__get_xy_coordinates())
            self.__z_coordinates = np.zeros(self.__coordinates[0].shape,dtype=np.uint)
            return self._get_coordinates()

    def _get_initial(self):
        return self.__initial

    def _get_evaluators(self,expressions,settings,compile_to_c = None,**kwargs):
        import expresso.pycas as pc

        if not isinstance(expressions,(list,tuple)):
            return_single = True
            expressions = [expressions]
        else:
            return_single = False

        sb = settings.simulation_box
        args = (sb.xi,sb.yi,sb.zi) if self.ndim == 2 else (sb.xi,sb.zi)

        if self.ndim == 1:
            expressions = [settings.get_optimized(expr.subs(sb.y,sb.fy)) for expr in expressions]
        else:
            expressions = [settings.get_optimized(expr) for expr in expressions]

        if 'return_type' not in kwargs:
            kwargs['return_type'] = pc.Types.Complex

        def is_constant(expr):
	    try:
            	expr.N()
		return True
	    except:
 		return False   

        definitions = [pc.FunctionDefinition('f%s' % i,args,expr,**kwargs)
                       for i,expr in enumerate(expressions) if not is_constant(expr)]

        if compile_to_c == None:
 	    compile_to_c = self.ndim > 1       
 
        if not compile_to_c:
            lib = pc.ncompile(*definitions)
        else:
            lib = pc.ccompile(*definitions)
        
	def get_constant_expression(expr):
	    try:
	        c = float(expr.N(20))
            except:
                c = complex(expr.N(20))
            import numpy as np
	    def constant_expression(*args,**kwargs):
                res = kwargs.pop('res',None)
		if res is not None:
		    res.fill(c)
		    return res
	        return np.full(args[0].shape,c)      
            return constant_expression

        res = [ getattr(lib,'f%s' % i) if hasattr(lib,'f%s' % i) else get_constant_expression(expressions[i]) for i in range(len(expressions)) ]

        if return_single:
            return res[0]
        return res
















