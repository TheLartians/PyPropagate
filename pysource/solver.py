
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
    _downscale:                 int [optional, defaults to 1]
                                How much the result should be downscaled
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
        
    def _get_downscale(self):
        try:
            return self._downscale
        except AttributeError:
            return 1
        
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
        
    def _get_box_size(self,axis,downscaled=True):
        downscale = self._get_downscale()

        if axis == 0:
            return self._nt/downscale+1 if downscaled else self._nt
        if axis == 1:
            return self._nx/downscale if downscaled else self._nx
        if axis == 2:
            return self._ny/downscale if downscaled else self._ny
        if axis == 3:
            return self._nz/downscale if downscaled else self._nz
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
        #if self._get_downscale() != 1:
        #    res = rebin(res,self._get_nd_box_size()[1:])
        return CoordinateNDArray(res,self._get_nd_boundary()[1:],self._get_nd_axis_symbols()[1:],self._get_transform())
    
    def set_field(self,field):
        
        if isinstance(field,CoordinateNDArray):
            field = field.data
        
        for x,xi,xj in zip(self._get_nd_axis_symbols()[1:],self._get_nd_box_size(downscaled=False)[1:],field.shape):
            if xi != xj:
                raise ValueError('Field size in %s direction (%s) doesn\'t match size defined in settings: %s' % (x,xj,xi))
        self._set_field(field)
            
    def reset(self):
        self._reset()
        self._i = 0
    
    def step(self):
        for updater in self._updaters.values():
            updater(self._i,self._get_field())
        self._step()
        self._i += 1
        
    def run(self,*args,**kwargs):
        """Simulate for _nt steps and return the resulting CoordinateNDArray."""
        return self.run_slice(self._get_nd_axis_symbols(self.ndim)[1:],*args,**kwargs)

    def run_slice(self,axis,display_progress=True,autohide_progress=False,slice_positions=None):
        """Simulate for _nt steps and return the resulting CoordinateNDArray."""
        
        import numpy as np
        from progressbar import ProgressBar

        self.reset()
        
        if not isinstance(axis, (list,tuple)):
            axis = [axis]
        
        nd = self.ndim
        
        all_axis = self._get_nd_axis_symbols(nd)
        
        indices = [0]
                
        for ax in axis:
            if ax not in all_axis[1:]:
                raise ValueError('no space axis %s defined' % ax)
            for i in range(nd+1):
                if all_axis[i] == ax:
                    indices.append(i)
                    break
        
        indices = sorted(set(indices))
        
        box_size = [self._get_box_size(i) for i in indices]
        downscale = self._get_downscale()
        
        field = np.zeros(box_size , dtype = self.dtype)
        
        if slice_positions == None:
            slice_positions = []
            for i in range(1,nd+1):
                slice_positions.append(downscale * self._get_box_size(i)/2)
        
        if len(slice_positions) != nd:
            raise ValueError('dimension mismatch')
        
        args = [slice(None) if i in indices[1:] else slice_positions[i-1] for i in range(1,nd+1)]
                                        
        if downscale == 1:
            def get_field():
                return self._get_field()[args]
        else:
            def get_field():
                return rebin(self._get_field()[args],box_size[1:])
        
        field[0] = get_field()

        
        steps = range(1,box_size[0])
        
        if display_progress == True:
            steps = ProgressBar(steps, title='Simulation running. Step',autohide=autohide_progress)
        
        for i in steps:
            if downscale == 1:
                self.step()
                field[i] = get_field()
            else:
                for j in range(downscale):
                    self.step()
                    field[i] += get_field()
                field[i] /= downscale

        
        axis_symbols = [self._get_axis_symbol(i) for i in indices]
        boundary = [self._get_boundary(i) for i in indices]
        
        res = CoordinateNDArray(field,boundary,axis_symbols, self._get_transform())
        return res.transpose(res.axis[1:] + [res.axis[0]])    
    
    
    
    