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

    from .. import units
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
        if axis_name in 't,z':
            setattr(sb,"%smin" % axis_name,0)
            setattr(sb,"%smax" % axis_name,size)
        else:
            setattr(sb,"%smin" % axis_name,-size/2)
            setattr(sb,"%smax" % axis_name,size/2)

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

        from ..units import get_unit
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
    from ..settings import Settings
    settings = Settings('settings for solving the paraxial differential equation')
    add_simulation_box_category(settings)
    settings.simulation_box.export(settings.symbols)
    add_partial_differential_equation_category(settings)
    settings.partial_differential_equation.export(settings.symbols)
    return settings

def add_wave_equation_category(settings):
    from .. import units
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
        we.omega = value / hbar

    import types
    we.add_method('set_energy',set_energy)

    return we

def create_paraxial_wave_equation_settings(fresnel_compatible = False):

    from ..settings import Settings
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


def create_2D_paraxial_settings_with_parameter(parameter_name = 'a'):
    from ..settings import Settings

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
    from .. import units
    from ..plot import expression_to_array
    from ..coordinate_ndarray import CoordinateNDArray
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

def create_2D_paraxial_frequency_settings():
    from ..settings import Settings

    settings = Settings()
    sb = add_simulation_box_category(settings,['x','omega','z'])
    pde = add_partial_differential_equation_category(settings,(sb.x,sb.omega,sb.z))
    we = add_wave_equation_category(settings)

    sb.export(settings.symbols,warn=False)
    we.export(settings.symbols,warn=False)

    settings.symbols.add_key('u0',pde.u0)
    settings.symbols.add_key('u_boundary',pde.u_boundary)

    pde.omega0 = (sb.omegamax + sb.omegamin)/2

    pde.A = 1/(2j*we.k*we.n)
    pde.C = 0
    pde.F = -1j*we.k*(we.n-1)

    return settings
