

def get_refractive_indices(formula,density,min_energy,max_energy,steps):
    """
get refractive indices for a compound material determined by the chemical formula and density [g/cm^3] and energy range [keV].
For room temperature density of elements set density to None.
    """
    import xraylib
    import numpy as np
   
    if density == None: 
        element_numbers = {'Ru': 44, 'Re': 75, 'Rf': 104, 'Rg': 111, 'Ra': 88, 'Rb': 37, 'Rn': 86, 'Rh': 45, 'Be': 4, 'Ba': 56, 'Bh': 107, 'Bi': 83, 'Bk': 97, 'Br': 35, 'Uuh': 116, 'H': 1, 'P': 15, 'Os': 76, 'Es': 99, 'Hg': 80, 'Ge': 32, 'Gd': 64, 'Ga': 31, 'Uub': 112, 'Pr': 59, 'Pt': 78, 'Pu': 94, 'C': 6, 'Pb': 82, 'Pa': 91, 'Pd': 46, 'Cd': 48, 'Po': 84, 'Pm': 61, 'Hs': 108, 'Uuq': 114, 'Uup': 115, 'Uus': 117, 'Uuo': 118, 'Ho': 67, 'Hf': 72, 'K': 19, 'He': 2, 'Md': 101, 'Mg': 12, 'Mo': 42, 'Mn': 25, 'O': 8, 'Mt': 109, 'S': 16, 'W': 74, 'Zn': 30, 'Eu': 63, 'Zr': 40, 'Er': 68, 'Ni': 28, 'No': 102, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10, 'Np': 93, 'Fr': 87, 'Fe': 26, 'Fm': 100, 'B': 5, 'F': 9, 'Sr': 38, 'N': 7, 'Kr': 36, 'Si': 14, 'Sn': 50, 'Sm': 62, 'V': 23, 'Sc': 21, 'Sb': 51, 'Sg': 106, 'Se': 34, 'Co': 27, 'Cm': 96, 'Cl': 17, 'Ca': 20, 'Cf': 98, 'Ce': 58, 'Xe': 54, 'Lu': 71, 'Cs': 55, 'Cr': 24, 'Cu': 29, 'La': 57, 'Li': 3, 'Tl': 81, 'Tm': 69, 'Lr': 103, 'Th': 90, 'Ti': 22, 'Te': 52, 'Tb': 65, 'Tc': 43, 'Ta': 73, 'Yb': 70, 'Db': 105, 'Dy': 66, 'Ds': 110, 'I': 53, 'U': 92, 'Y': 39, 'Ac': 89, 'Ag': 47, 'Uut': 113, 'Ir': 77, 'Am': 95, 'Al': 13, 'As': 33, 'Ar': 18, 'Au': 79, 'At': 85, 'In': 49}
        density = xraylib.ElementDensity(element_numbers[formula])

    return np.array([xraylib.Refractive_Index(formula, E, density).conjugate() for E in np.linspace(min_energy,max_energy,steps)])

def get_henke_refractive_indices(material ,min_energy ,max_energy ,steps ,density=-1 ,uniform_distance = False):

    if min_energy < 0 and max_energy < 0:
        return get_refraction_indices(material ,abs(min_energy) ,abs(max_energy) ,steps ,density ,uniform_distance)

    if min_energy > max_energy:
        return get_refraction_indices(material ,max_energy ,min_energy ,steps ,density ,uniform_distance)[::-1]

    max_steps = 499

    if steps > max_steps:
        import numpy as np
        dn = (max_energy - min_energy ) /(steps - 1)
        current_max = min_energy + max_steps * dn
        missing = max( steps -max_steps ,3)
        return np.append \
            (get_refraction_indices(material ,min_energy ,current_max ,max_steps, density, uniform_distance), \
            get_refraction_indices(material, current_max + dn, current_max + dn * missing, missing, density,
                                   uniform_distance), axis=0)

    from mechanize import Browser
    br = Browser()

    br.open("http://henke.lbl.gov/optical_constants/getdb.html")

    br.select_form(nr=0)

    br.form['Formula'] = material
    br.form['Density'] = str(density)
    br.form['Min'] = str(min_energy) if not uniform_distance else str(min_energy - 1)
    br.form['Max'] = str(max_energy) if not uniform_distance else str(max_energy + 1)
    br.form['Npts'] = str(steps - 1)
    br.form['Output'] = ['Text File']

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
        n_values = [complex(1 - float(v[1]), -float(v[2])) for v in betadelta]
    except:
        betadelta = []

    if len(betadelta) != steps:
        raise RuntimeError(
            'error retrieving refractive index for %s (E from %s to %s in %s steps)\nserver response: %s' % (
                material, min_energy, max_energy, steps, res))

    if uniform_distance:
        from scipy.interpolate import interp1d
        import numpy as np
        int_f = interp1d(E_values, n_values)
        E_values = np.linspace(min_energy, max_energy, steps)
        interpolated = int_f(E_values)
        return interpolated

    return zip(E_values, n_values)


def create_material(name, settings, density = None):
    '''
    density in g/cm^3 or None for room temperature of element
    '''

    import xraylib
    import expresso.pycas as pc

    nname = 'n_%s' % name

    if not settings.has_category('refractive_indices'):
        settings.create_category('refractive_indices')
    r = settings.refractive_indices

    def init_material(settings):
        from .. import units
        import numpy as np

        sb = settings.simulation_box
        r = settings.refractive_indices
        omega = settings.wave_equation.omega

        try:
            N = settings.get_as(sb.Nomega, int)
            omegamin, omegamax = (sb.omegamin, sb.omegamax)

            EminExpr = omegamin * units.hbar / units.keV
            EmaxExpr = omegamax * units.hbar / units.keV
            
            Emin = settings.get_as(EminExpr,float)
            Emax = settings.get_as(EmaxExpr,float)

            omega_dependent = True
        except:
            N = 1
            E = (units.hbar * omega / units.keV)
            omega_i = 1
            omega_dependent = False
            try:
                Enum = settings.get_as(E, float)
                Emin = Enum
                Emax = Enum
            except:
                setattr(r, nname, None)
                return

        key = (nname, N, Emin, Emax, density)
        if not hasattr(r, '_cache'):
            r.add_attribute('_cache', {})
        else:
            if key in r._cache:
                setattr(r, nname, r._cache[key])
                return
        if omega_dependent:
            narr = pc.array(nname, np.array(get_refractive_indices(name, density, Emin, Emax, N)))
            setattr(r, nname, narr(sb.omegai))
            r._cache[key] = narr(sb.omegai)
        else:
            val = get_refractive_indices(name, density, Emin, Emax, N)[0]
            setattr(r, nname, val)
            r._cache[key] = val

    settings.initializers["init_" + nname] = init_material

    if r.has_name(nname):
        return getattr(r, nname)
    n = r.create_key(nname, pc.Function(nname)(settings.wave_equation.omega))

    settings.numerics.add_key(nname, n)
    return n

