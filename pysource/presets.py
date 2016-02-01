


def add_simulation_box_symbols(settings):

    sb = settings.create_category("simulation_box",info="parameters and dimensions of the simulation box")

    from pycas import Symbol,symbols,Function,pi,I,Types

    x = sb.create_key("x",Symbol("x",type = Types.Real),info="field coordinate for 1D simulations")
    y = sb.create_key("y",Symbol("y",type = Types.Real),info="second coordinate for 2D simulations")
    z = sb.create_key("z",Symbol("z",type = Types.Real),info="propagation direction")

    sb.x = x
    sb.y = y
    sb.z = z

    nx = sb.create_key("nx", Symbol("n_x",type = Types.Integer,positive=True),info="voxels in x direction minus the boundary conditions")
    ny = sb.create_key("ny", Symbol("n_y",type = Types.Integer,positive=True),info="voxels in y direction minus the boundary conditions")
    nz = sb.create_key("nz", Symbol("n_z",type = Types.Integer,positive=True),info="voxels in z direction minus the boundary condition")

    Nx = sb.create_key("Nx", Symbol("N_x",type = Types.Integer,positive=True),info="voxels in x direction")
    Ny = sb.create_key("Ny", Symbol("N_y",type = Types.Integer,positive=True),info="voxels in y direction")
    Nz = sb.create_key("Nz", Symbol("N_z",type = Types.Integer,positive=True),info="voxels in z direction")

    sb.Nx = nx+2
    sb.Ny = ny+2
    sb.Nz = nz+1
    sb.lock('Nx','defined by nx')
    sb.lock('Nz','defined by ny')
    sb.lock('Nz','defined by nz')

    sx = sb.create_key("sx",Symbol("s_x",type = Types.Real,positive=True),info="simulation box size in x direction")
    sy = sb.create_key("sy",Symbol("s_y",type = Types.Real,positive=True),info="simulation box size in y direction")
    sz = sb.create_key("sz",Symbol("s_z",type = Types.Real,positive=True),info="simulation box size in z direction")

    sb.create_key("dx",Symbol("d_x",type = Types.Real,positive=True),info="voxel size in x direction")
    sb.create_key("dy",Symbol("d_y",type = Types.Real,positive=True),info="voxel size in y direction")
    settings.simulation_box.create_key("dz",Symbol("d_z",type = Types.Real,positive=True),info="voxel size in z direction")

    xmin = sb.create_key("xmin", Symbol("x_min",type = Types.Real),info="x value at the lower simulation box boundary")
    xmax = sb.create_key("xmax", Symbol("x_max",type = Types.Real),info="x value at the upper simulation box boundary")
    ymin = sb.create_key("ymin", Symbol("y_min",type = Types.Real),info="y value at the lower simulation box boundary")
    ymax = sb.create_key("ymax", Symbol("y_max",type = Types.Real),info="y value at the upper simulation box boundary")
    zmin = sb.create_key("zmin", Symbol("z_min",type = Types.Real),info="z value at the lower simulation box boundary")
    zmax = sb.create_key("zmax", Symbol("z_max",type = Types.Real),info="z value at the upper simulation box boundary")

    sb.dx = sx/(nx+1)
    sb.dy = sy/(ny+1)
    sb.dz = sz/nz
    sb.lock('dx','defined by sx and nx')
    sb.lock('dy','defined by sy and ny')
    sb.lock('dz','defined by sz and nz')

    sb.sx = xmax - xmin
    sb.sy = ymax - ymin
    sb.sz = zmax - zmin
    sb.lock('sx','defined by xmin and xmax')
    sb.lock('sy','defined by ymin and ymax')
    sb.lock('sz','defined by zmin and zmax')

    sb.create_key("fy",Symbol("fixed y",type = Types.Real),info="fixed y value for 2D simulations")
    sb.fy = (ymin + ymax)/2

    sb.create_key("downscale", Symbol("downscale"), 1,info="factor by which the resulting field of the simulation will be scaled down")

    sb.lock()

    def set_physical_size(sx,sy,sz):
        'Sets the physical box size of the simulation in x, y and z direction'
        sb.unlock('xmin')
        sb.unlock('xmax')
        sb.unlock('ymin')
        sb.unlock('ymax')
        sb.unlock('zmin')
        sb.unlock('zmax')
        sb.unlock('sx')
        sb.unlock('sy')
        sb.unlock('sz')

        sb.sx,sb.sy,sb.sz = (sx,sy,sz)
        sb.xmin = -sb.sx/2
        sb.ymin = -sb.sy/2
        sb.zmin = 0
        sb.xmax = sb.sx/2
        sb.ymax = sb.sy/2
        sb.zmax = sb.sz

        sb.lock('xmin','defined by sx')
        sb.lock('xmax','defined by sx')
        sb.lock('ymin','defined by sy')
        sb.lock('ymax','defined by sy')
        sb.lock('zmin','defined by sz')
        sb.lock('zmax','defined by sz')

        from units import get_unit
        defined = set()

        for s in (sx,sy,sz):
            unit = get_unit(s)
            if unit is None or unit in defined:
                continue
            defined.add(unit)
            unit_name = str(unit)
            if not settings.unitless.has_name(unit_name):
                settings.unitless.create_key(unit_name,unit)
            setattr(settings.unitless,unit_name,(2*unit/s).evaluate())


    def set_voxel_size(nx,ny,nz):
        'Sets the voxe size of the simulation in x, y and z direction'
        voxel_size = (nx,ny,nz)
        sb.nx,sb.ny,sb.nz = voxel_size

    def set_simulation_box(physical_size,voxel_size):
        """Sets the simulation box size using the physical_size and voxel_size arguments which are 3-tupels containing the simulation box dimensions in x, y, and z directions."""
        set_physical_size(*physical_size)
        set_voxel_size(*voxel_size)

    sb._set_attribute('set_physical_size', set_physical_size)
    sb._set_attribute('set_voxel_size', set_voxel_size)
    sb._set_attribute('set', set_simulation_box)


def add_paraxial_equation_symbols(settings):
    from pycas import Symbol,Function,Types

    pe = settings.create_category("paraxial_equation",info="parameters of the paraxial differential equation")
    s = settings.simulation_box

    pe.create_key("F",Function("F")(s.x,s.y,s.z),info="F(x,y,z) Parameter of the differential equation")
    pe.create_key("A",Symbol("A",type=Types.Complex),info="A Parameter of the differential equation")

    pe.create_key("u0",Function("u_0")(s.x,s.y,s.z),info="field initial condition")
    pe.create_key("u_boundary",Function("u_boundary")(s.x,s.y,s.z),info="field boundary condition")

    pe.lock()

def create_paraxial_settings():
    from .settings import Settings
    settings = Settings('settings for solving the paraxial differential equation')
    add_simulation_box_symbols(settings)
    settings.simulation_box.export(settings.symbols)
    add_paraxial_equation_symbols(settings)
    settings.paraxial_equation.export(settings.symbols)
    return settings

def add_wave_equation_symbols(settings):
    import units
    from pycas import Function,Symbol,Types,pi

    s = settings.simulation_box

    we = settings.create_category("wave_equation",info="parameters for solving the wave equation")

    n = we.create_key("n",Function("n")(s.x,s.y,s.z))
    k = we.create_key("k",Symbol("k",type = Types.Complex))
    wavelength = we.create_key("wavelength",Symbol(r"lambda",type = Types.Real,positive=True))

    we.wavelength = 2*pi/k
    omega = we.create_key("omega",Symbol(r"omega",type = Types.Real,positive=True),2*pi/wavelength)

    we.lock('wavelength','defined by wave number')
    we.lock('omega','defined by wavelength')

    def set_energy(value):
        if not we.has_name('E'):
            we.create_key("E",Symbol("E",type = Types.Real,positive=True),info='Wave energy')
            settings.symbols.add_key("E",we.E)
            we.k = we.E / (units.hbar*units.c)
            we.lock('k','defined by energy')
        we.E = value

    we._set_attribute('set_energy',set_energy)

def create_paraxial_wave_equation_settings():

    from .settings import Settings
    settings = Settings('settings for solving the paraxial differential equation')

    add_simulation_box_symbols(settings)
    settings.simulation_box.export(settings.symbols)

    add_wave_equation_symbols(settings)
    settings.wave_equation.export(settings.symbols)

    add_paraxial_equation_symbols(settings)
    settings.paraxial_equation.export(settings.symbols)

    from pycas import I

    pe = settings.paraxial_equation
    s = settings.symbols
    pe.F = -I*s.k/2*(s.n**2-1)
    pe.A = -I/(2*s.k)

    pe.lock('F','defined by wave equation')
    pe.lock('A','defined by wave equation')
    s.remove_name("F")
    s.remove_name("A")

    return settings

def set_plane_wave(settings):
    """Sets the intial conditions to a plane wave with intensity 1.
    The boundary are set to the index of refraction at z=0."""
    from pycas import exp

    s = settings.simulation_box
    pe = settings.paraxial_equation

    pe.u0 = 1
    pe.u_boundary = exp(pe.F.subs(s.z,0)*s.z)
