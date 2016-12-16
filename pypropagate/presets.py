
def add_coordinate(settings,category,name):

    import expresso.pycas as pc

    x = category.create_key(name,pc.Symbol(name,type=pc.Types.Real),info='coordinate')
    xmin = category.create_key('%smin' % name,pc.Symbol('%s_min' % name),info='minimum value')
    xmax = category.create_key('%smax' % name,pc.Symbol('%s_max' % name),info='maximum value')
    xi = pc.Symbol('%s_i' % name,type=pc.Types.Natural)
    xif = category.create_key('%si' % name,pc.Function('%s_i' % name)(x),info='numerical index')
    Nx = category.create_key('N%s' % name,pc.Symbol('N_%s' % name,type=pc.Types.Natural),info='numerical steps')
    sx = category.create_key('s%s' % name,pc.Symbol('s_%s' % name),info='total size')
    dx = category.create_key('d%s' % name,pc.Symbol('Delta %s' % name),info='step size')

    setattr(category,name, xmin + xi * dx)
    setattr(category,'%si' % name,(x - xmin)/dx)

    setattr(category,'s%s' % name,xmax - xmin)
    setattr(category,'d%s' % name,sx/(Nx-1))

    category.lock(name,'defined by %si' % name)
    category.lock('s%s' % name,'defined by %smin and %smax' % (name,name))
    category.lock('d%s' % name,'defined by s%s and N%s' % (name,name))

    settings.unitless.add_key('%s_coordinate' % name,x)

def add_time_symbols(settings):
    settings.simulation_box.unlock()
    add_coordinate(settings,settings.simulation_box,'t')
    sb  = settings.simulation_box
    sb.tmin = 0
    sb.export(settings.symbols,warn=False)

    import units
    settings.unitless.create_key(None,units.s,10**10)

    if settings.has_category('wave_equation'):
        settings.wave_equation.create_function('u0',(sb.x,sb.y,sb.z,sb.t))
        settings.partial_differential_equation.u0 = settings.wave_equation.u0.subs(sb.t,0)

def add_simulation_box_category(settings,coords = ['x','y','z']):
    import expresso.pycas as pc

    sb = settings.create_category("simulation_box",info="parameters and dimensions of the simulation box")

    for c in coords:
        add_coordinate(settings,sb,c)

    class CoordinateAttrs:
        def __init__(self,sb,name):
            self.name = name
            self.symbol = getattr(sb,name)
            self.min = getattr(sb,name+'min')
            self.max = getattr(sb,name+'max')
            self.step = getattr(sb,'d'+name)
            self.steps = getattr(sb,'N'+name)
            self.size = getattr(sb,'s'+name)
            self.index =  pc.Symbol('%s_i' % name,type=pc.Types.Natural)
        def __repr__(self):
            return "<%s Attrs>" % self.name

    sb.add_attribute("coordinates", tuple(CoordinateAttrs(sb, c) for c in coords))
    sb.add_attribute("coordinate_dict", {getattr(sb, c):d for c, d in zip(coords, sb.coordinates)})

    sb.lock()

    def set_size(sb,axis_name,size):
        sb.unlock("%smin" % axis_name)
        sb.unlock("%smax" % axis_name)

        if axis_name in 't,z':
            setattr(sb,"%smin" % axis_name,0)
            setattr(sb,"%smax" % axis_name,size)
        else:
            setattr(sb,"%smin" % axis_name,-size/2)
            setattr(sb,"%smax" % axis_name,size/2)

        sb.lock("%smin" % axis_name,'defined by s%s' % axis_name)
        sb.lock("%smax" % axis_name,'defined by s%s' % axis_name)

    def set_vsize(sb,axis_name,size):
        setattr(sb,"N%s" % axis_name,size)

    def set_physical_size(sb,*sizes):
        if len(coords) != len(sizes):
            raise ValueError('number of arguments does not match coordinates %s' % coords)
        for c,s in zip(coords,sizes):
            sb.set_size(c,s)

    def set_voxel_size(sb,*sizes):
        if len(coords) != len(sizes):
            raise ValueError('number of arguments does not match coordinates %s' % coords)
        for c,s in zip(coords,sizes):
            sb.set_vsize(c,s)

    def set_method(sb,physical_size,voxel_size):
        """Sets the simulation box size using the physical_size and voxel_size arguments which are 3-tupels containing the simulation box dimensions in x, y, and z directions."""
        sb.set_physical_size(*physical_size)
        sb.set_voxel_size(*voxel_size)

    sb.add_method('set_size', set_size)
    sb.add_method('set_vsize', set_vsize)

    sb.add_method('set_physical_size', set_physical_size)
    sb.add_method('set_voxel_size', set_voxel_size)
    sb.add_method('set', set_method)

    def make_unitless(settings):
        sb = settings.simulation_box

        from units import get_unit
        defined = set()

        for s in settings.get_numeric(tuple(getattr(sb,'s'+c) for c in coords)):
            unit = get_unit(s)

            if unit is None or unit in defined:
                continue

            unit_name = str(unit)

            value = (unit/s).evaluate()

            if unit.function == pc.fraction:
                unit = unit.args[0]
                value = 1/value

            if  unit.is_function:
                continue

            if unit is None or unit in defined or unit.is_function:
                continue

            if not settings.unitless.has_name(unit_name):
                settings.unitless.create_key(unit_name,unit)

            setattr(settings.unitless,unit_name,value)

            defined.add(unit)

    settings.initializers['make_unitless'] = make_unitless

    return sb

def add_partial_differential_equation_category(settings,coordinates = None):
    import expresso.pycas as pc
    sb = settings.simulation_box
    pde = settings.create_category('partial_differential_equation',short_name='PDE',info="parameters of the partial differential equation")

    arg_attrs = [sb.coordinate_dict[x] for x in coordinates] if coordinates is not None else sb.coordinates
    pde.add_attribute('coordinates', arg_attrs)

    x,y,z = [a.symbol for a in arg_attrs]

    dx,dy,dz = [s.step for s in arg_attrs]
    args = (x,y,z)

    pde.create_function('A',args )
    pde.create_function('C',args ,pde.A)
    pde.create_function('F',args )

    pde.create_function('ra',args ,pde.A*dz/dx**2,info="finite difference paramter")
    pde.create_function('rc',args ,pde.C*dz/dy**2,info="finite difference paramter")
    pde.create_function('rf',args ,pde.F*dz/2,info="finite difference paramter")

    pde.create_function('u',args ,info='solution to the PDE')
    pde.lock('u')

    pde.add_attribute('equation', pc.equal(pc.derivative(pde.u, z), pde.A * pc.derivative(pc.derivative(pde.u, x), x) + pde.C * pc.derivative(pc.derivative(pde.u, y), y) + pde.F * pde.u))

    pde.create_key('u0',pc.Function('u_0_PDE')(*args ),info="field initial condition")
    pde.create_function('u_boundary',args ,info="field boundary condition");

    pde.create_key(arg_attrs[1].name+'0',pc.Symbol(y.name+'_0_PDE'),0,info='Value to which the %s is set for 1D solvers' % y.name)

    pde.lock()

    return pde

def create_paraxial_settings():
    from .settings import Settings
    settings = Settings('settings for solving the paraxial differential equation')
    add_simulation_box_category(settings)
    settings.simulation_box.export(settings.symbols)
    add_partial_differential_equation_category(settings)
    settings.partial_differential_equation.export(settings.symbols)
    return settings

def add_wave_equation_category(settings):
    import units
    from expresso.pycas import Function,Symbol,Types,pi

    s = settings.simulation_box

    c = settings.numerics.create_key('c',Symbol('c'),units.c,info='speed of light')
    h = settings.numerics.create_key('h',Symbol('h'),units.h,info='Planck\'s constant')
    hbar = settings.numerics.create_key('hbar',Symbol('hbar'),units.hbar,info='reduced Planck\'s constant')

    we = settings.create_category("wave_equation",info="parameters for solving the wave equation",short_name="WE")

    n = we.create_key("n",Function("n")(*[c.symbol for c in s.coordinates]))

    omega = we.create_key('omega',Symbol('omega'),info='angular wave frequency') if not hasattr(s,'omega') else we.add_key('omega',s.omega)
    wavelength = we.create_key('wavelength',Function("lambda")(omega),settings.numerics.c*2*pi/omega,info='vacuum wavelength')
    k = we.create_key("k",Function("k")(omega),omega/settings.numerics.c,info='wave number')
    E = we.create_key("E",Function("E")(omega),omega * settings.numerics.hbar,info='photon energy')

    we.lock('k','defined by omega')
    we.lock('E','defined by omega')

    def set_energy(we,value):
        we.unlock('E')
        we.omega = we.E / hbar
        we.lock('omega','defined by energy')
        we.E = value

    import types
    we.add_method('set_energy',set_energy)

    return we

def create_paraxial_wave_equation_settings(fresnel_compatible = False):

    from .settings import Settings
    import expresso.pycas as pc

    settings = Settings('settings for solving the paraxial differential equation')

    add_simulation_box_category(settings)
    settings.simulation_box.export(settings.symbols)

    add_wave_equation_category(settings)
    settings.wave_equation.export(settings.symbols)

    add_partial_differential_equation_category(settings)
    settings.partial_differential_equation.export(settings.symbols)

    from expresso.pycas import I

    pe = settings.partial_differential_equation
    s = settings.symbols
    pe.F = -1j*s.k*(s.n-1)

    if fresnel_compatible:
        pe.A = 1/(2j*s.k)
    else:
        pe.A = 1/(2j*s.k*s.n)

    pe.lock('F','defined by wave equation')
    pe.lock('A','defined by wave equation')
    s.remove_name("F")
    s.remove_name("A")

    return settings

def set_1D_boundary_condition(settings):
    from expresso.pycas import exp

    s = settings.simulation_box
    pe = settings.partial_differential_equation

    pe.u_boundary = pe.u0.subs(s.z,s.zmin) * exp(pe.F.subs(s.z,s.zmin)*s.z)


def set_plane_wave_initial_conditions(settings):
    """Sets the boundary conditions to a plane wave with intensity 1.
    The boundary are set to the index of refraction at z=0."""

    s = settings.simulation_box
    pe = settings.partial_differential_equation

    pe.u0 = 1
    set_1D_boundary_condition(settings)


def create_next_settings(old_settings):
    settings = old_settings.copy(copy_initializers = False,copy_updaters = False)

    sb = settings.simulation_box
    sb.unlock('zmin')
    sb.unlock('zmax')
    sb.unlock('sz')

    sb.zmin = old_settings.get_numeric(sb.zmax)
    sb.zmax = sb.zmin + sb.sz

    sb.lock('zmin','defined by previous simulation')
    sb.lock('zmax','defined by zmin and sz')

    return settings

def add_padding(array,factor,mode = 'edge',**kwargs):
    import numpy as np
    from coordinate_ndarray import CoordinateNDArray

    padding_points = [[int(x*factor)]*2 for x in array.data.shape]
    new_data = np.pad(array.data,padding_points,mode,**kwargs)

    extension = [d*p[0] for d,p in zip(array._dbounds,padding_points)]
    new_bounds = [(b-i,e+i) for i,(b,e) in zip(extension,array.bounds)]

    return CoordinateNDArray(new_data,new_bounds,array.axis,array.evaluate)

def set_initial(settings,initial_array):
    import expresso.pycas as pc
    from coordinate_ndarray import CoordinateNDArray

    if isinstance(initial_array,CoordinateNDArray):
        initial = pc.array("initial",initial_array.data)
    else:
        initial = pc.array("initial",initial_array)

    sb = settings.simulation_box

    if tuple(initial_array.axis) == (sb.x,):
        settings.partial_differential_equation.u0 = initial(sb.xi)
    elif tuple(initial_array.axis) == (sb.x,sb.y):
        settings.partial_differential_equation.u0 = initial(sb.yi,sb.xi)
        sb.Ny = initial_array.shape[1]
        if isinstance(initial_array,CoordinateNDArray):
            sb.unlock('ymin')
            sb.unlock('ymax')
            sb.unlock('sy')
            sb.ymin = initial_array.bounds[1][0]
            sb.ymax = initial_array.bounds[1][1]
            sb.sy = sb.ymax - sb.ymin
            sb.lock('ymin','defined by initial array')
            sb.lock('ymax','defined by initial array')
            sb.lock('sy','defined by ymin and ymax')

    else:
        raise ValueError('initial array axis must be (x,) or (x,y)')

    sb.Nx = initial_array.shape[0]

    if isinstance(initial_array,CoordinateNDArray):
        sb.unlock('xmin')
        sb.unlock('xmax')
        sb.unlock('sx')
        sb.xmin = initial_array.bounds[0][0]
        sb.xmax = initial_array.bounds[0][1]
        sb.sx = sb.xmax - sb.xmin
        sb.lock('xmin','defined by initial array')
        sb.lock('xmax','defined by initial array')
        sb.lock('sx','defined by xmin and xmax')

def get_refraction_indices(material,min_energy,max_energy,steps,density=-1,uniform_distance = False):

    if min_energy < 0 and max_energy < 0:
        return get_refraction_indices(material,abs(min_energy),abs(max_energy),steps,density,uniform_distance)

    if min_energy > max_energy:
        return get_refraction_indices(material,max_energy,min_energy,steps,density,uniform_distance)[::-1]

    max_steps = 499

    if steps > max_steps:
        import numpy as np
        dn = (max_energy - min_energy)/(steps - 1)
        current_max = min_energy + max_steps * dn
        missing = max(steps-max_steps,3)
        return np.append(get_refraction_indices(material,min_energy,current_max ,max_steps,density,uniform_distance), \
               get_refraction_indices(material,current_max + dn,current_max + dn * missing,missing,density,uniform_distance),axis=0)

    from mechanize import Browser
    br = Browser()

    br.open("http://henke.lbl.gov/optical_constants/getdb.html")

    br.select_form(nr=0)
    
    br.form[ 'Formula' ] = material
    br.form[ 'Density' ] = str(density)
    br.form[ 'Min' ] = str(min_energy) if not uniform_distance else str(min_energy-1)
    br.form[ 'Max' ] = str(max_energy) if not uniform_distance else str(max_energy+1)
    br.form[ 'Npts' ] = str(steps-1)
    br.form[ 'Output' ] = ['Text File']

    res = br.submit().read()
    
    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def get_numbers(line):
        return [float(v) for v in line.split(' ') if is_number(v)]
    try:
        betadelta = [get_numbers(line) for line in res.split('\n') if len(get_numbers(line)) == 3]
        E_values = [float(v[0]) for v in betadelta]
        n_values = [complex(1-float(v[1]),-float(v[2])) for v in betadelta]
    except:
        betadelta = []

    if len(betadelta) != steps:
        raise RuntimeError('error retrieving refractive index for %s (E from %s to %s in %s steps)\nserver response: %s' % (
            material,min_energy,max_energy,steps,res))
    
    if uniform_distance:
        from scipy.interpolate import interp1d
        import numpy as np
        int_f = interp1d(E_values,n_values)
        E_values = np.linspace(min_energy,max_energy,steps)
        interpolated = int_f(E_values)
        return interpolated
    
    return zip(E_values,n_values)

def create_material(name,settings,density=-1):
    '''
    density in gm/cm^3
    '''

    import expresso.pycas as pc

    nname = 'n_%s' % name

    if not settings.has_category('refractive_indices'):
        settings.create_category('refractive_indices')
    r = settings.refractive_indices

    def init_material(settings):
        import units
        import numpy as np

        sb = settings.simulation_box
        r = settings.refractive_indices
        omega = settings.wave_equation.omega

        try:
            N = settings.get_as(sb.Nomega,int)
            omegamin,omegamax = (sb.omegamin,sb.omegamax)

            EminExpr = omegamin * units.hbar / units.eV
            EmaxExpr = omegamax * units.hbar / units.eV

            omega_dependent = True
        except:
            N = 3
            E = (units.hbar * omega / units.eV)
            omega_i = 1
            omega_dependent = False
        try:
            Enum = settings.get_as(E,float)
            Emin = Enum-1
            Emax = Enum+1
        except:
            setattr(r,nname,None)
            return

        key = (nname,N,Emin,Emax,density)
        if not hasattr(r,'_cache'):
            r.add_attribute('_cache', {})
        else:
            if key in r._cache:
                setattr(r,nname,r._cache[key])
                return
        if omega_dependent:
            narr = pc.array(nname,np.array(get_refraction_indices(name,Emin,Emax,N,density,True)))
            setattr(r,nname,narr(sb.omegai))
            r._cache[key] = narr(sb.omegai)
        else:
            val = get_refraction_indices(name,Emax,Emin,3,density,True)[1]
            setattr(r,nname,val)
            r._cache[key] = val

    settings.initializers["init_" + nname] = init_material

    if r.has_name(nname):
        return getattr(r,nname)
    n = r.create_key(nname,pc.Function(nname)(settings.wave_equation.omega))

    settings.numerics.add_key(nname,n)
    return n

def create_2D_paraxial_frequency_settings():
    from .settings import Settings

    settings = Settings()
    sb = add_simulation_box_category(settings,['x','omega','z'])
    pde = add_partial_differential_equation_category(settings,(sb.x,sb.omega,sb.z))
    we = add_wave_equation_category(settings)

    sb.export(settings.symbols,warn=False)
    we.export(settings.symbols,warn=False)

    settings.symbols.add_key('u0',pde.u0)
    settings.symbols.add_key('u_boundary',pde.u_boundary)

    pde.A = 1/(2j*we.k*we.n)
    pde.C = 0
    pde.F = -1j*we.k*(we.n-1)

    return settings

def create_2D_paraxial_settings_with_parameter(parameter_name = 'a'):
    from .settings import Settings

    settings = Settings()
    sb = add_simulation_box_category(settings,['x',parameter_name,'z'])
    pde = add_partial_differential_equation_category(settings,(sb.x,getattr(sb,parameter_name),sb.z))
    we = add_wave_equation_category(settings)

    sb.export(settings.symbols,warn=False)
    we.export(settings.symbols,warn=False)

    settings.symbols.add_key('u0',pde.u0)
    settings.symbols.add_key('u_boundary',pde.u_boundary)

    pde.A = 1/(2j*we.k*we.n)
    pde.C = 0
    pde.F = -1j*we.k*(we.n-1)

    return settings

def create_2D_frequency_settings_from(settings):
    import expresso.pycas as pc
    import units
    from .plot import expression_to_array
    from .coordinate_ndarray import CoordinateNDArray
    import numpy as np

    settings.initialize()
    settings = settings.copy()

    sb = settings.simulation_box
    we = settings.wave_equation
    omega0 = settings.get_numeric(we.omega)
    we.unlock('omega')
    we.omega = None

    freq_settings = create_paraxial_settings()
    pde = freq_settings.partial_differential_equation

    freq_settings.simulation_box.unlock()
    freq_settings.simulation_box.remove_name('y')
    omega = freq_settings.simulation_box.create_key('y',pc.Symbol('omega'))

    fsb = freq_settings.simulation_box
    freq_settings.simulation_box.y = fsb.ymin + fsb.yi * fsb.dy

    n = settings.get_numeric(we.n.subs(sb.y,0)).subs(omega,abs(omega - omega0))

    p = freq_settings.create_category('paramters',short_name="p")
    p.create_function('n',(sb.x,sb.y,sb.z,omega),n)
    p.create_symbol('k',(omega - omega0)/units.c)
    p.create_symbol('k_0',omega0/units.c)

    pde.A = 1j/(2*p.k_0)
    pde.C = 0
    pde.F = 1j/(2*p.k_0) * (p.n**2 * p.k**2 - p.k_0**2)

    xmin = settings.get_numeric(sb.xmin)
    xmax = settings.get_numeric(sb.xmax)

    freq_settings.simulation_box.xmin = xmin
    freq_settings.simulation_box.xmax = xmax
    freq_settings.simulation_box.Nx = settings.get_as(sb.Nx,int)

    omegamin = settings.get_numeric(-2*pc.pi*sb.Nt/sb.st)
    omegamax = settings.get_numeric( 2*pc.pi*sb.Nt/sb.st)

    freq_settings.simulation_box.ymin = omegamin
    freq_settings.simulation_box.ymax = omegamax
    freq_settings.simulation_box.Ny = settings.get_as(sb.Nt,int)

    freq_settings.simulation_box.zmin = settings.get_numeric(sb.zmin)
    freq_settings.simulation_box.zmax = settings.get_numeric(sb.zmax)
    freq_settings.simulation_box.Nz = settings.get_as(sb.Nz,int)

    freq_settings.unitless.create_key('s',units.s,settings.get_as( 2/(omegamax - omegamin)/units.s , float ) )

    u0 = expression_to_array(settings.wave_equation.u0.subs([(sb.y, 0), (sb.z, sb.zmin)]), settings, axes=(sb.x, sb.t))

    u0hat = CoordinateNDArray( np.fft.fftshift( np.fft.fft(u0.data,axis=1) , axes=(1,) ) ,
                               [(xmin,xmax),( omegamin, omegamax )] ,
                               [sb.x,omega])

    set_initial(freq_settings,u0hat)
    freq_settings.partial_differential_equation.u_boundary = 0

    return freq_settings

def fourier_transform(array,axis,new_axis,inverse=False):
    import numpy as np
    from numpy.fft import fftshift,ifftshift
    from expresso.pycas import pi
    from .coordinate_ndarray import CoordinateNDArray

    try:
        from pyfftw.interfaces.numpy_fft import fft,ifft
    except ImportError:
        from numpy.fft import fft,ifft

    axi = array.axis.index(axis)

    if not inverse:
        new_data = fftshift(fft(array.data,axis=axi),axes=[axi])
        new_data *= 1/np.sqrt(2*np.pi)  
    else:
        new_data = ifft(ifftshift(array.data,axes=[axi]),axis=axi)
        new_data *= np.sqrt(2*np.pi) 

    sw = array.bounds[axi][1] - array.bounds[axi][0]
    tmin,tmax = array.evaluate((-(pi*array.shape[axi])/sw,
                                 (pi*array.shape[axi])/sw))

    new_bounds = [(b[0],b[1]) if i!=axi else (tmin,tmax) for i,b in enumerate(array.bounds)]
    new_axis = [a  if i!=axi else new_axis for i,a in enumerate(array.axis)]

    return CoordinateNDArray(new_data,new_bounds,new_axis,array.evaluate)

def inverse_fourier_transform(*args):
    return fourier_transform(*args,inverse=True)

def u_from_utilde(field,omega0,s=None,a=None):

    if a is None and s is None:
	raise ValueError("a or s needs to be specified")
    if a is not None and s is not None:
        raise ValueError("only one value for a/s can be specified")
    if a is not None:
        s = 1/(1-a)

    import numpy as np
    import expresso.pycas as pc
    import units

    omegamin,omegamax = field.bounds[1]
    zmin, zmax = field.bounds[2]
    sz = zmax - zmin

    ukmin = float(field.evaluate( (omegamin-omega0)/units.c*sz ))
    ukmax = float(field.evaluate( (omegamax-omega0)/units.c*sz ))
    uzmin = 0
    uzmax = 1

    #print ukmin
    #print ukmax

    nz,ik = np.meshgrid(np.linspace(uzmin,uzmax,field.shape[2]),np.linspace(1j*ukmin,1j*ukmax,field.shape[1]))
    factor = np.exp(-ik*nz/s)
    
    transform = field * factor
    del factor,nz,ik

    return fourier_transform(transform, field.axis[1], pc.Symbol('t'), inverse=True)



