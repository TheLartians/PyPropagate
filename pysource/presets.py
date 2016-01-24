
import units
import sympy
from .coordinate_ndarray import CoordinateNDArray

def set_plane_wave(settings):
    """Sets the intial conditions to a plane wave with intensity 1. 
    The boundary are set to the index of refraction at z=0."""
    
    from sympy import exp
    z = settings.simulation_box.z
    settings.wave_equation.u0 = 1
    settings.wave_equation.u0_boundary = exp(settings.finitedifferences.F.subs(z,0)*z)

def add_wave_equation_symbols(settings):
    
    from sympy import Symbol,Function,I,pi,exp
    s = settings.simulation_box
    
    settings.create_category("wave_equation",info="parameters for solving the wave equation")
    
    m = settings.base_units.create_key("m",units.m,units.m/s.sx)

    n = settings.wave_equation.create_key("n",Function("n")(s.x,s.y,s.z))
    k = settings.wave_equation.create_key("k",Symbol("k",complex = True))
    wavelength = settings.wave_equation.create_key("wavelength",Symbol(r"\lambda",real=True,positive=True))
    omega = settings.wave_equation.create_key("omega",Symbol(r"\omega",real=True,positive=True),2*pi/wavelength)
    E = settings.wave_equation.create_key("E",Symbol("E",real=True,positive=True))
    
    psi0 = settings.wave_equation.create_key("psi0",Function(r"\psi_0")(s.x,s.y,s.z))

    if settings.has_category('symbols'):
        settings.wave_equation.export(settings.symbols)
    
    hbar = settings.wave_equation.create_key("hbar",Symbol(r"\hbar",real=True),units.hbar)
    
    c = settings.wave_equation.create_key("c",Symbol("c",real=True),units.c)

    if settings.has_category('finitedifferences'):
        settings.finitedifferences.F = -I*k/2*(n**2-1)
        settings.finitedifferences.A = -I/(2*k)
        settings.finitedifferences.lock('F','defined by wave equation')
        settings.finitedifferences.lock('A','defined by wave equation')

        if settings.has_category('symbols'):
            settings.symbols.remove_name("F")
            settings.symbols.remove_name("A")
        
        settings.wave_equation.add_key('u0',settings.finitedifferences.u0,info = 'paraxial field')
        settings.wave_equation.add_key('u0_boundary',settings.finitedifferences.u0_boundary,info = 'paraxial field boundary condition')
        settings.wave_equation.psi0 = exp(-I*k*s.z)*settings.wave_equation.u0
        settings.wave_equation.lock('psi0','defined by u0')
        
    settings.wave_equation.k = E / (hbar*c)
    settings.wave_equation.wavelength = 2*pi/k
    
def add_finitedifference_symbols(settings):
            
    settings.create_category("finitedifferences",info="parameters of the differential equation")
    s = settings.simulation_box
    
    from sympy import Symbol,Function

    settings.finitedifferences.create_key("F",Function("F")(s.x,s.y,s.z),info="F(x,y,z) Parameter of the differential equation")
    settings.finitedifferences.create_key("A",Symbol("A",complex=True),info="A Parameter of the differential equation")
    settings.finitedifferences.create_key("u0",Function("u_0")(s.x,s.y,s.z),info="field initial condition")
    settings.finitedifferences.create_key("u0_boundary",Function("u_{0,\mathrm{boundary}}")(s.x,s.y,s.z),info="field boundary condition")
    
    settings.finitedifferences.lock()

    if settings.has_category('symbols'):
        settings.finitedifferences.export(settings.symbols)

def add_propagator_symbols(settings):
                
    from categorized_dictionary import Category
        
    settings.create_category("simulation_box",info="parameters and dimensions of the simulation box")

    from sympy import Symbol,symbols,Function,pi,I

    x = settings.simulation_box.create_key("x",Symbol("x",real = True),info="field coordinate for 1D simulations")
    y = settings.simulation_box.create_key("y",Symbol("y",real = True),info="second coordinate for 2D simulations")
    z = settings.simulation_box.create_key("z",Symbol("z",real = True),info="propagation direction")
    
    settings.simulation_box.x = Symbol("x",real = True)
    settings.simulation_box.y = Symbol("y",real = True)
    settings.simulation_box.z = Symbol("z",real = True)

    nx = settings.simulation_box.create_key("nx", Symbol("n_x",integer=True,positive=True),info="voxels in x direction")
    ny = settings.simulation_box.create_key("ny", Symbol("n_y",integer=True,positive=True),info="voxels in y direction")
    nz = settings.simulation_box.create_key("nz", Symbol("n_z",integer=True,positive=True),info="voxels in z direction")

    sx = settings.simulation_box.create_key("sx",Symbol("s_x",real=True,positive=True),info="simulation box size in x direction")
    sy = settings.simulation_box.create_key("sy",Symbol("s_y",real=True,positive=True),info="simulation box size in y direction")
    sz = settings.simulation_box.create_key("sz",Symbol("s_z",real=True,positive=True),info="simulation box size in z direction")

    settings.simulation_box.create_key("dx",Symbol("d_x",real=True,positive=True),info="voxel size in x direction")
    settings.simulation_box.create_key("dy",Symbol("d_y",real=True,positive=True),info="voxel size in y direction")
    settings.simulation_box.create_key("dz",Symbol("d_z",real=True,positive=True),info="voxel size in z direction")

    xmin = settings.simulation_box.create_key("xmin", Symbol("x_\mathrm{min}",real=True),info="smallest x value inside simulation box")
    xmax = settings.simulation_box.create_key("xmax", Symbol("x_\mathrm{max}",real=True),info="largest x value inside simulation box")
    ymin = settings.simulation_box.create_key("ymin", Symbol("y_\mathrm{min}",real=True),info="smallest y value inside simulation box")
    ymax = settings.simulation_box.create_key("ymax", Symbol("y_\mathrm{max}",real=True),info="largest y value inside simulation box")
    zmin = settings.simulation_box.create_key("zmin", Symbol("z_\mathrm{min}",real=True),info="smallest z value inside simulation box")
    zmax = settings.simulation_box.create_key("zmax", Symbol("z_\mathrm{max}",real=True),info="largest z value inside simulation box")
    
    settings.simulation_box.dx = sx/(nx-1)
    settings.simulation_box.dy = sy/(ny-1)
    settings.simulation_box.dz = sz/(nz-1)
    
    settings.simulation_box.sx = xmax - xmin
    settings.simulation_box.sy = ymax - ymin
    settings.simulation_box.sz = zmax - zmin
    settings.simulation_box.lock('sx','defined by xmin and xmax')
    settings.simulation_box.lock('sy','defined by ymin and ymax')
    settings.simulation_box.lock('sz','defined by zmin and zmax')
    
    settings.simulation_box.create_key("downscale", Symbol("\mathrm{downscale}"), 1,info="factor by which the resulting field of the simulation will be scaled down")
    
    settings.simulation_box.lock()
    
    if settings.has_category('symbols'):
        settings.simulation_box.export(settings.symbols)
    
    return settings

def create_finitedifference_settings():
    from settings import Settings
    settings = Settings('settings for a finite difference simulation')
    add_propagator_symbols(settings)
    add_finitedifference_symbols(settings)
    return settings
        
def as_wavelengths(field,settings):
    """Will convert field boundries to wavelengths."""
    s = settings.wave_equation
    bounds = [[settings.get(b/s.wavelength,substitute_numeric=True)*s.wavelength for b in B] for B in field.bounds]
    return CoordinateNDArray(field.data,bounds,field.axis,field.transform)
    
def as_meters(field,settings):
    """Will convert field boundries to meters."""
    field = field.soft_copy()
    field.bounds = [[settings.get(b/units.m,substitute_numeric=True)*units.m for b in B] for B in field.bounds]
    return field
