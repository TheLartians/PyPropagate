from unittest import TestCase


from pypropagate import *


def analytical_circular_waveguide(settings):

    from pypropagate.coordinate_ndarray import CoordinateNDArray

    import scipy.optimize
    import scipy.integrate
    from scipy.special import j0,j1,k0,k1
    import numpy as np
    from numpy import sqrt,ceil,pi,sign,cos,sin,exp,sum,abs
    import warnings

    s = settings.symbols

    kn = settings.get_as(s.k,float)
    n1 = settings.get_as(settings.waveguide.n_1,complex)
    n2 = settings.get_as(settings.waveguide.n_2,complex)
    a  = settings.get_as(settings.waveguide.r,float)

    # We will determine um for all guided modes as the roots of the characteristic equation:
    def char(um):
        wm = sqrt(V**2-um**2)
        return um*j1(um)/j0(um) - wm*k1(wm)/k0(wm)

    # Dimensionless waveguide paramter
    V = kn*a*sqrt(n1.real**2 - n2.real**2)

    # V determines the number of guided modes
    N = ceil(V/pi)

    # Now find all N roots of the characteristic equation
    u_values = []
    segments = N

    # Silence division warnings which might occur since we will search for sign changes in arbitrary intervals
    warnings.filterwarnings('ignore')

    while len(u_values) < N:
        segments = segments*2
        u_values = []
        for interval in [(xi,xi+V/segments) for xi in np.linspace(0,V*(1-1./segments),segments)]:
            va,vb = char(np.array(interval))
            # Check if values for interval boundaries are finite and have opposite sign
            if np.isfinite(va) and np.isfinite(vb) and sign(va) != sign(vb):
                # There might be a root in the interval. Find it using Brent's method!
                um,r = scipy.optimize.brentq(char,interval[0],interval[1],full_output=True)
                # Check that point converged and value is small (we might have converged to a pole instead)
                if r.converged and abs(char(um))<1:
                    u_values.append(um)
    u_values = np.array(u_values)

    # Reactivate warnings
    warnings.resetwarnings()

    # Define the guided modes
    def psi(um,r):
        wm = sqrt(V**2-um**2)
        return np.piecewise(r,[r<a,r>=a],[lambda r:j0(um*r/a)/j0(um),lambda r:k0(wm*r/a)/k0(wm)])

    # Normalize the modes by integrating the intensity of the guided modes over R^2
    B_values = [(2*pi*scipy.integrate.quad(lambda r:abs(psi(um,r))**2*r,0,np.inf)[0])**(-0.5) for um in u_values]

    # Project the modes onto the plain wave initial condition
    c_values = [2*pi*scipy.integrate.quad(lambda r:psi(um,r)*B*1*r,0,np.inf)[0] for um,B in zip(u_values,B_values)]

    # Determine propagation constants
    def mu(r):
        if r<a: return kn*2*n1.imag
        else: return kn*2*n2.imag

    # Integrate mu in two parts to aviod discontinuity
    mu_values = np.array([
            2*pi*scipy.integrate.quad(lambda r:abs(psi(um,r)*B)**2*mu(r)*r,0,a)[0] +
            2*pi*scipy.integrate.quad(lambda r:abs(psi(um,r)*B)**2*mu(r)*r,a,np.inf)[0]
                         for um,B in zip(u_values,B_values)])

    beta_values = sqrt(n1.real**2*kn**2-u_values**2/a**2)

    # The full solution is the superposition of the guided modes
    def field(x,r):
        solution = None
        for i in range(len(u_values)):
            mode = c_values[i] * psi(u_values[i],np.abs(r)) * B_values[i] * exp( (-1j*beta_values[i]-mu_values[i]/2) * x)
            if solution is None: solution = mode
            else: solution += mode
        return solution

    x_values = np.linspace(*settings.get_as((s.zmin,s.zmax,s.Nz),float))
    r_values = np.linspace(*settings.get_as((s.xmin,s.xmax,s.Nx),float))

    data = np.conjugate(field(*np.meshgrid(x_values,r_values)) * exp(1j*kn*x_values))

    sx = settings.get_numeric(s.sx)
    sz = settings.get_numeric(s.sz)

    res = CoordinateNDArray(data,[(-sx/2,sx/2),(0,sz)],(s.x,s.z),settings.get_numeric_transform())
    return res


class TestWaveguide(TestCase):

    def test_analytical(self):
        settings = presets.create_paraxial_wave_equation_settings()
        presets.set_plane_wave_initial_conditions(settings)
        s = settings.symbols

        wg = settings.create_category('waveguide')
        wg.create_symbol('n_1')
        wg.n_1 = 1

        wg.create_symbol('n_2')
        wg.n_2 = 1-6.006E-06+6.32E-07j

        wg.create_symbol('r')
        wg.r = 24*units.nm

        wg.create_symbol('l')
        wg.l = 0.6*units.mm

        s.n = pc.piecewise((wg.n_1,s.x**2+s.y**2<wg.r**2),(wg.n_2,True))

        settings.wave_equation.set_energy(12*units.keV)
        settings.simulation_box.set((0.1*units.um,0.1*units.um,0.8*units.mm),(500,500,1000))

        propagator = propagators.FiniteDifferencesPropagator2D(settings)
        an_field = analytical_circular_waveguide(settings)
        fd_field = propagator.run_slice()[:,0]

        self.assertLess(abs(an_field - fd_field)[:,0.1*units.mm:].max(),0.1)
