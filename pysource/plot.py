
def get_metric_prefix(numbers):
    import numpy as np
                  
    if not isinstance(numbers,list):
        numbers = list(numbers)
        
    def get_exponent(number):
        return int(np.log10(np.abs(number))) if number != 0 else 0
    
    from units import metric_prefixes
    
    exponents = [get_exponent(number) for number in numbers]
    largest = max(exponents,key=lambda x:abs(x))
    
    closest = min(metric_prefixes, key=lambda x:abs(x[1]+1-largest))
    return (closest[2],10**closest[1])

def get_unitless_bounds(array):

    from .units import get_unit

    bounds = []
        
    for l,r in array.bounds:
        unit = get_unit(l)
        if unit == None:
            unit = get_unit(r)
        try:
            if unit == None:
                bounds.append((float(l),float(r),1)) 
            else:
                bounds.append((float(l/unit),float(r/unit),unit)) 
        except TypeError:
            raise ValueError('Cannot convert to unitless expression: %s with unit: %s' % ((l,r),unit))
    
    return bounds
    
def image_plot(carr,ax = None,figsize = None,title = None, **kwargs):
    import matplotlib.pyplot as plt

    # fix missing \text support
    from pycas import latex as rlatex
    latex = lambda x:rlatex(x).replace(r'\text',r'\mathrm')

    fig = None
    if ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    
    if title:
        ax.set_title(title)
    
    e = get_unitless_bounds(carr)
    xprefix,xfactor = get_metric_prefix(e[1][:2])
    yprefix,yfactor = get_metric_prefix(e[0][:2])
        
    extent = [float(e[1][0])/xfactor,float(e[1][1])/xfactor,float(e[0][1])/yfactor,float(e[0][0])/yfactor]
    image = ax.imshow(carr.data, extent= extent, aspect='auto', **kwargs )
    ax.set_ylabel("$%s$ [$%s %s$]" % (latex(carr.axis[0]),yprefix,latex(e[0][2])))
    ax.set_xlabel("$%s$ [$%s %s$]" % (latex(carr.axis[1]),xprefix,latex(e[1][2])))
    
    if fig:
        fig.colorbar(image)

    if ax == None:
        plt.show()

    return image
        
def line_plot(carr,ax = None,ylabel = None,figsize = None,title = None,**kwargs):
    import matplotlib.pyplot as plt
    import numpy as np

    # fix missing \text support
    from pycas import latex as rlatex
    latex = lambda x:rlatex(x).replace(r'\text',r'\mathrm')

    fig = None
    if ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    
    if title:
        ax.set_title(title)
    
    e = get_unitless_bounds(carr)[0]
    
    prefix,factor = get_metric_prefix(e[:2])
    
    lines = ax.plot(np.linspace(float(e[0])/factor,float(e[1])/factor,carr.data.shape[0]),carr.data, **kwargs)
    ax.set_xlabel("$%s$ [$%s %s$]" % (latex(carr.axis[0]),prefix,latex(e[2])))
    if ylabel: ax.set_ylabel(ylabel)

    if ax == None:
        plt.show()

    return lines[0]
        
def plot(carr,*args,**kwargs):
    """
    Simple plot function for 1D and 2D coordinate arrays. If the data is complex, the absolute square value of the data will be plottted.
    
    Parameters
    -----------
    carr: coordinate array
          the input data
    
    **kwargs: additional parameters to be passed to the plot functions
    
    Returns
    --------
    plot: output of ax.plot for 1D and ax.imshow for 2D arrays
    
    """
    
    import numpy as np
    if not np.can_cast(carr.data.dtype, np.float): carr = abs(carr)**2
    if len(carr.axis) == 1: return line_plot(carr,*args,**kwargs)
    elif len(carr.axis) == 2: return image_plot(carr,*args,**kwargs)
    else: raise ValueError("input array must be one or two dimensional")
    
    
  