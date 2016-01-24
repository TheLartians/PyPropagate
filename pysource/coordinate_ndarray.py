
import numpy as np

class CoordinateNDArray(object):
    """
    A wrapper for numpy arrays containing physical bound properties for the array. This enables slicing based on physical coordinates rather than array indices.
    
    Members / Initializer Parameters:
    ---------------------------------
    data:       numpy.ndarray (required)
                Array with the n-dimensional data
    bounds:     list of tuples [(min,max),(min,max),...] (optional, defaults to [(0,1)]*n)
                A list with n elements containing the physical boundaries of the array
    axis:       list (optional, defaults to range(n))
                A list with n elements containing the names of the axes (eg. ["x","y","z"])
    transform:  function-like (optional, defaults to lambda x:x)
                When accessing the array, arguments are transformed by this function
    
    Example:
    --------
    >>> data = numpy.linspace(-1,1,21)
    >>> arr = CoordinateNDArray(data,[[-1,1]],['x'])
    >>> arr[0:1].bounds
    [(0.0, 1.0)]
    >>> arr[0:1].data
    array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ])
    
    WARNING: 
    reverse iteration (slicing with a negative step size) is not supported yet
    """
    
    @staticmethod
    def __default_transform(x):
        return x
    
    def __init__(self,data,bounds = None,axis = None,transform=None):        
        self.data = np.array(data)
        n = len(self.data.shape)
        self.bounds = [(b[0],b[1]) for b in bounds] if bounds != None else [(0,1)]*n
        self.axis = [a for a in axis] if axis != None else range(n)
        self.transform = transform if transform != None else self.__default_transform
        if not len(self.data.shape) == len(self.bounds) == len(self.axis):
            raise ValueError("dimensions do not match")
        self.__dbounds = [(b[1] - b[0])/float(n-1) if n != 1 else 0 for b,n in zip(self.bounds,self.data.shape)]
    
    def __constrain_coordinate(self,val,axis):
        val = self.transform(val)
        return min(max(val,self.bounds[axis][0]),self.bounds[axis][1])
    
    def __test_slice(self,sliced,axis):
        try:
            for v in [sliced.start,sliced.stop,sliced.step]:
                if v:
                    int(self.transform(v)/self.__dbounds[axis])
        except:
            raise ValueError("Cannot slice array: cannot convert to int: %s/%s." % (self.transform(v),self.__dbounds[axis]))
    
    def __convert_slice(self,sliced,axis):
        s = [sliced.start,sliced.stop,sliced.step]
        if s[0] == None: s[0] = self.bounds[axis][0]
        else: s[0] = self.__constrain_coordinate(s[0],axis)
        if s[1] == None: s[1] = self.bounds[axis][1]
        else: s[1] = self.__constrain_coordinate(s[1],axis)
        if s[2] == None: s[2] = self.__dbounds[axis]
        else:
            if s[2]> 0: s[2] = max(min(self.bounds[axis][1] - self.bounds[axis][0],s[2]), self.__dbounds[axis])
            else: s[2] = min(max(-self.bounds[axis][1] + self.bounds[axis][0],s[2]), -self.__dbounds[axis])
        if s[2] < 0: raise RuntimeError("negative step size not implemented yet")
        return s
    
    def __get_numpy_slice(self,s,axis):
        start = int(((s[0] - self.bounds[axis][0])/self.__dbounds[axis] + 0.5))
        stop  = int(((s[1] - self.bounds[axis][0])/self.__dbounds[axis] + 1.5))
        step  = int((s[2]/self.__dbounds[axis] + 0.5))
        return slice(start,stop,step)
    
    def __get_coordinate_from_index(self,idx,axis):
        return self.bounds[axis][0] + idx * self.__dbounds[axis]
    
    def __repr__(self):
        return "<CoordinateNDArray, axis: %r, bounds %r, shape: %r, dtype: %s>" % (self.axis,self.bounds,self.data.shape,self.data.dtype)
    
    def __get_bounds_for_slice(self,sliced,axis):
        if not isinstance(sliced,slice): sliced = slice(sliced,sliced)
        self.__test_slice(sliced,axis)
        s = self.__convert_slice(sliced,axis)
        n = self.__get_numpy_slice(s,axis)
        n_data_end = n.start + (n.stop - n.start-1) # * n.step
        slice_bounds = [self.__get_coordinate_from_index(n.start,axis),self.__get_coordinate_from_index(n_data_end,axis)]
        return slice_bounds,n
    
    def __getitem__(self,sliced):
        if not isinstance(sliced,tuple): sliced = (sliced,)
        assert len(sliced) <= len(self.axis), "too many indices for array"
        data_slice = []
        new_bounds = []
        new_axis = []
        for axis,sliced_part in enumerate(sliced):
            slice_bounds,n = self.__get_bounds_for_slice(sliced_part,axis)
            if n.start+1 != n.stop:
                new_bounds.append(slice_bounds)
                new_axis.append(self.axis[axis])
                data_slice.append(n)
            else:
                data_slice.append(n.start)
        new_bounds += self.bounds[len(sliced):]
        new_axis += self.axis[len(sliced):]
        data_slice = tuple(data_slice)
        if len(new_axis) == 0: return self.data.__getitem__(data_slice)
        return CoordinateNDArray(self.data.__getitem__(data_slice),new_bounds,new_axis,self.transform)
    
    def copy(self):
        """Creates a full copy."""
        bounds_copy = [(b[0],b[1]) for b in self.bounds]
        axis_copy = [a for a in self.axis]
        return CoordinateNDArray(self.data.copy(),bounds_copy,axis_copy)
    
    def soft_copy(self):
        """Creates a copy referencing the same numpy data."""
        bounds_copy = [(b[0],b[1]) for b in self.bounds]
        axis_copy = [a for a in self.axis]
        return CoordinateNDArray(self.data,bounds_copy,axis_copy)

    def is_compatible(self,value):
        """Axis and bounds checks."""
        return self.axis == value.axis and self.data.shape == value.data.shape # and self.bounds == value.bounds Not possible due to sympy rounding errors
    
    def __convert_arg(self,arg):
        if isinstance(arg,CoordinateNDArray):
            if not self.is_compatible(arg): raise ValueError("incompatible coordinate arrays")
            else: return arg.data
        return arg
    
    def __mutated_copy(self,function,*args,**kwargs):
        args = [self.__convert_arg(arg) for arg in args]
        kwargs = {key:self.__convert_arg(kwargs[key]) for key in kwargs}
                
        if "axis" in kwargs:
            if kwargs["axis"] not in self.axis:
                raise ValueError("array does not contain axis: %s" % kwargs["axis"])
            idx = self.axis.index(kwargs["axis"])
            kwargs["axis"] = idx
            
        res = function(self.data, *args, **kwargs)
            
        if "axis" in kwargs and res.shape != self.data.shape:
            idx = kwargs["axis"]
            bounds_copy = [(b[0],b[1]) for i,b in enumerate(self.bounds) if i != idx]
            axis_copy = [a for i,a in enumerate(self.axis) if i != idx]
        else:
            bounds_copy = [(b[0],b[1]) for b in self.bounds]
            axis_copy = [a for a in self.axis]
        
        if not isinstance(res,np.ndarray):
            return res
        
        return CoordinateNDArray(function(self.data, *args, **kwargs),bounds_copy,axis_copy,self.transform)
    
    def transpose(self,axis = None):
        if axis == None:
            axis = list(reversed(self.axis))
        np_axes = [self.axis.index(ax) for ax in axis]
        axis_copy = [a for a in axis]
        bounds_copy = [(self.bounds[i][0],self.bounds[i][1]) for i in np_axes]
        return CoordinateNDArray(self.data.transpose(np_axes),bounds_copy,axis_copy,self.transform)
    
    def apply_numpy_function(self,function,*args,**kwargs):
        "Apply a numpy function to the data array."
        return self.__mutated_copy(function,*args,**kwargs)
    
    def __add__(self, value) :
        return self.__mutated_copy(np.add,value)
    def __sub__(self, value) :
        return self.__mutated_copy(np.subtract,value)
    def __mul__(self, value) :
        return self.__mutated_copy(np.multiply,value)
    def __div__(self, value) :
        return self.__mutated_copy(np.divide,value)
    def __pow__(self, value) :
        return self.__mutated_copy(np.power,value)
    def __lt__(self, value) :
        return self.__mutated_copy(np.less,value)
    def __le__(self, value) :
        return self.__mutated_copy(np.less_equal,value)
    def __gt__(self, value) :
        return self.__mutated_copy(np.greater,value)
    def __ge__(self, value) :
        return self.__mutated_copy(np.greater_equal,value)
    def __eq__(self, value) :
        return self.__mutated_copy(np.equal,value)
    def __ne__(self, value) :
        return self.__mutated_copy(np.not_equal,value)
    def __abs__(self) :
        return self.__mutated_copy(np.abs)

    def _numpy_function_wrapper(self,function):
        def wrapped_function(*args, **kwargs):
            return self._CoordinateNDArray__mutated_copy(function,*args, **kwargs)
        return wrapped_function
    
    def __getattr__(self, key):
        try:
            return self.__getattribute__(key)
        except AttributeError:
            # fallback to numpy attributes
            res = self.data.__getattribute__(key)
            if hasattr(res, '__call__'):
                def wrapped_member(self,*args,**kwargs):
                    return res(*args,**kwargs)
                def wrapped_function(*args,**kwargs):
                    return self.__mutated_copy(wrapped_member,*args,**kwargs)
                return wrapped_function
            return res
        
    def save_to_file(self,path,warn_transform = True):
        "Save the array to a given Path. Please note that custom transform functions cannot be saved."
        
        import cPickle as pickle
        
        with open(path, 'wb') as output:  
            as_array = [self.data,self.bounds,self.axis]
            pickle.dump(as_array, output, pickle.HIGHEST_PROTOCOL)
        
        if self.transform != self.__default_transform and warn_transform:
            import warnings
            warnings.warn('cannot save custom transform function')
        
    @classmethod
    def load_from_file(cls,path):
        "Load an array from a given Path."

        import cPickle as pickle

        with open(path, 'rb') as input:
            loaded = pickle.load(input)
        
        return CoordinateNDArray(*loaded)
            
def numpy_function_wrapper(function):
    def wrapped_function(arr=None,*args, **kwargs):
        if isinstance(arr,CoordinateNDArray):
            return arr._CoordinateNDArray__mutated_copy(function,*args, **kwargs)
        else: return function(arr,*args,**kwargs)
    return wrapped_function
        
class WrappedNumpy(object):
    
    """
    Numpy wrapper for CoordinateNDArrays. Wrappes functions of a module that accept numpy arrays to functions that also accept CoordinateNDArrays.
    
    Example:
    --------
    >>> import numpy
    >>> np = WrappedNumpy(numpy) # Wrapped numpy module that accepts CoordinateNDArrays

    WARNING: 
    submodules are not wrapped.
    """
    
    def __init__(self,module = np):
        self._module = module
        
    def __dir__(self):
        return self._module.__dict__.keys()
        
    def __getattr__(self, name):
        if not self._module.__dict__.has_key(name):
            raise AttributeError(name)
        numpy_function = self._module.__dict__.get(name)
        if type(numpy_function) == type:
            return numpy_function
        if hasattr(numpy_function, '__call__'):
            f = numpy_function_wrapper(numpy_function)
            f.__doc__ = "Wrapped numpy function " + numpy_function.__doc__ 
            return f
        else:
            return numpy_function
    
    
