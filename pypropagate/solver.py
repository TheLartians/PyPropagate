
from .coordinate_ndarray import CoordinateNDArray 
from numpy import complex128

def rebin(ndarray, new_shape):

    compression_pairs = [(d, c//d) for d,c in zip(new_shape, ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    
    ndarray = ndarray.reshape(flattened)
    
    for i in range(len(new_shape)):
        ndarray = ndarray.mean(-1*(i+1))
        
    return ndarray

class Solver(object):
    """
    Solver Base Class for simulations. Time axis is labeled t, space ax
    
    Necessary Attributes for 1D Solver
    ----------------------------------
    ndim:                       1 
                                Dimensions
    _x,_t:                      Symbols
                                Axis symbols
    _nx,_nt:                    int
                                Voxel counts
    _xmin,_xmax,_tmin,_tmax:    Symbols
                                Simulation Boundaries
    _transform:                 callable object [optional, defaults to lambda x:x]
                                Function that converts symbolic expression to numeric values
    dtype:                      numpy dtype
                                internal data type    
    _step():                    function
                                run one simulation step
    _get_field():               function
                                return the current field as a numpy array
    _set_field(field):          function
                                set the current field
    _reset():                   function
                                resets the simulation
    
    Notes
    -----
    For 2 and 3 dimensional solvers, define values for y and z as well. For higher dimensions implement the functions _get_boundary, _get_axis_symbol and _get_box_size which return the appropriate value for an integer axis index. Index convention: 0: t, 1: x, 2: y, 3: z.  
    
    """
    
    def __init__(self,settings):
        settings.initialize()
        #if len(settings.undefined_keys()) > 0:
        #    raise ValueError("there are undefined symbols: %s" % ", ".join([str(s) for s in settings.undefined_keys()]))
        self._updaters = settings.updaters.copy()
        self._i = 0
        self._transform = settings.get_numeric_transform()

    def _step(self):
        raise NotImplementedError('step function not implemented')
    
    def _reset(self):
        raise NotImplementedError('reset function not implemented')
    
    def _get_field(self):
        raise NotImplementedError('get field function not implemented')
        
    def _set_field(self,field):
        raise NotImplementedError('set field function not implemented')
    
    def _get_transform(self):
        try:
            return self._transform
        except AttributeError:
            return lambda x:x

    def _get_boundary(self,axis):
        if axis == 0:
            return (self._tmin,self._tmax)
        if axis == 1:
            return (self._xmin,self._xmax)
        if axis == 2:
            return (self._ymin,self._ymax)
        if axis == 3:
            return (self._zmin,self._zmax)
        raise IndexError('axis out of range, define custom _get_boundary')
    
    def _get_axis_symbol(self,axis):
        if axis == 0:
            return self._t
        if axis == 1:
            return self._x
        if axis == 2:
            return self._y
        if axis == 3:
            return self._z
        raise IndexError('axis out of range, define custom _get_axis_symbol')
        
    def _get_box_size(self,axis):
        if axis == 0:
            return self._nt
        if axis == 1:
            return self._nx
        if axis == 2:
            return self._ny
        if axis == 3:
            return self._nz
        raise IndexError('axis out of range, define custom _get_box_size')
        
    def _get_nd_box_size(self,n = None,**kwargs):
        if n == None:
            n = self.ndim
        res = []
        for i in range(n+1):
            res.append(self._get_box_size(i,**kwargs))
        return res
    
    def _get_nd_boundary(self,n = None):
        if n == None:
            n = self.ndim
        res = []
        for i in range(n+1):
            res.append(self._get_boundary(i))
        return res
        
    def _get_nd_axis_symbols(self,n = None):
        if n == None:
            n = self.ndim
        res = []
        for i in range(n+1):
            res.append(self._get_axis_symbol(i))
        return res
    
    def get_field(self):
        res = self._get_field()
        return CoordinateNDArray(res,self._get_nd_boundary()[1:],self._get_nd_axis_symbols()[1:],self._get_transform())

    def set_field(self,field):
        
        if isinstance(field,CoordinateNDArray):
            field = field.data
        
        for x,xi,xj in zip(self._get_nd_axis_symbols()[1:],self._get_nd_box_size()[1:],field.shape):
            if xi != xj:
                raise ValueError('Field size in %s direction (%s) doesn\'t match size defined in settings: %s' % (x,xj,xi))
        self._set_field(field)
            
    def reset(self):
        self._reset()
        self._i = 0

    def _call_updaters(self):
        for updater in self._updaters.values():
            updater(self)

    def step(self,callback=None):
        self._i += 1
        self._call_updaters()
        self._step()
        if callback is not None:
            callback(self)

    def run_slice(self,**kwargs):

        class RunSliceAgent:
            def __init__(self,parent,shape,kwargs):
                self.parent = parent
                self.kwargs = kwargs
                self.shape = shape
            def __getitem__(self,sliced):
                return self.parent._run_slice(sliced,**kwargs)

        agent_axis = self._get_nd_axis_symbols()
        agent_axis = agent_axis[1:] + [agent_axis[0]]
        agent_bounds = self._get_nd_boundary()
        agent_bounds= agent_bounds[1:] + [agent_bounds[0]]
        agent_shape = self._get_nd_box_size()
        agent_shape = agent_shape[1:] + [agent_shape[0]]

        slice_agent = CoordinateNDArray(RunSliceAgent(self,agent_shape,kwargs),agent_bounds,agent_axis,self._get_transform())
        return slice_agent

    def run(self,callback = None,display_progress=True, autohide_progress=False):
        from .progressbar import ProgressBar

        run_steps = range(1,self._get_box_size(0))
        if display_progress == True:
            run_steps = ProgressBar(run_steps, title='Simulation running. Step',autohide=autohide_progress)

        for i in run_steps:
            self.step(callback=callback)

    def _run_slice(self, sliced , display_progress=True, autohide_progress=False, callback = None):
        """Simulate for _nt steps and return the resulting CoordinateNDArray."""
        
        import numpy as np
        from .progressbar import ProgressBar

        self.reset()

        if not isinstance(sliced,(list,tuple)):
            sliced = [sliced]
        sliced = list(sliced)
        if len(sliced) < self.ndim :
            raise ValueError('not all dimension coordinates specified')
        if len(sliced) != self.ndim + 1:
            sliced.append(slice(0,self._get_box_size(0)))

        sliced = [sliced[-1]] + list(sliced[:-1])
        box_size = self._get_nd_box_size()

        sliced_indices = [(s.indices(b) if isinstance(s,slice) else (s,s+1,1)) for b,s in zip(box_size,sliced)]

        box_size = [ (s[1] - s[0])/s[2] if s[2]!=1 else s[1] - s[0] for s in sliced_indices]
        box_size = [b for b in box_size if b != 1]

        sliced = tuple([slice(*s) if s[0]+1 != s[1] else s[0] for s in sliced_indices[1:]])
        def get_field():
            return self._get_field().__getitem__(sliced)

        test = get_field()
        for i in range(len(test.shape)):
            box_size[i+1] = get_field().shape[i]

        i = 0

        start,stop,step = sliced_indices[0]

        for i in range(start):
            self.step(callback=callback)

        field = np.zeros(box_size[::-1] , dtype = self.dtype).transpose()

        field[0] = get_field()

        run_steps = range(1,box_size[0])
        if display_progress == True:
            run_steps = ProgressBar(run_steps, title='Simulation running. Step',autohide=autohide_progress)

        for j in run_steps:
            for k in range(step):
                self.step(callback=callback)
                i += 1
            field[j] = get_field()

        for i in range(i+1,self._nt+1):
            self.step(callback=callback)

        return field.transpose(range(1,len(field.shape)) + [0])
        #res = CoordinateNDArray(field, bounds, axis, self._get_transform())
        #return res.transpose(res.axis[1:] + [res.axis[0]])
    
    
    
    
