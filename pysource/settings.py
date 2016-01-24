from .categorized_dictionary import Category,CategorizedDictionary
import collections
from .timeout import timeout,TimeoutError

class Settings(CategorizedDictionary):
    """
    Stores the parameters for a FiniteDifference simulation. This class inherits from CategorizedDictionary. Initially, it contains four categories named finite_differences, simulation_box, numerics and base_units. See the info(name) member function for more information about the categories and keys. 
    There is an additional symbols member containing aliases to all keys. Adding a key with the same name will overwrite the name in the symbols member.
    The `get` function provides a method of transforming an expression using all known substitutions defined by the key/value pairs of this dictionary.
    """
    
    def __init__(self,create_categories = True):
        
        self.initializers = collections.OrderedDict()
        self.updaters = collections.OrderedDict()
        
        super(Settings,self).__init__()
        
        self.__doc__ = "Settings of a FiniteDifference simulation"
        
        if create_categories :
            self.create_category("numerics",info="pure numerical values")
            self.create_category("base_units",info="conversions to unitless coordinates")
            self.create_category("symbols",info="symbol aliases to frequently used symbols")

    def initialize(self):
        for initializer in self.initializers.values():
            initializer(self)
        
    def get_transform(self):
        return lambda x:self.get_numeric(x)

    def copy(self):
        "Return a soft copy of the Settings object."
        copy = Settings(create_categories = False)
        super(Settings,self).copy(copy = copy)
        copy.initializers = (self.initializers.copy())
        copy.updaters = (self.updaters.copy())
        return copy
                
    def _is_numeric(self,value):
        "Tests if a value should be classified as numeric"
        try:
            complex(value)
            return True
        except:
            pass
        
        from sympy.core.basic import Basic,preorder_traversal
        import units
        
        u = set([units.m,units.kg,units.s,units.A,units.K] + self.base_units.keys())
        
        if isinstance(value,Basic):
            for arg in preorder_traversal(value):
                if arg in u:
                    return True
                
        return False
                
    def _set_value(self,key,value):
        super(Settings,self)._set_value(key,value)
        
        is_numeric = self._is_numeric(value)
        
        if key not in self.numerics.keys() and is_numeric:
            self.numerics.add_key(None,key)
        elif key in self.numerics.keys() and not is_numeric:
            self.numerics.remove_key(key)
        
    @staticmethod
    def substitute_all(expr,substitutions,ignore = None,debug=False,simplify_timeout=1):
        """Recursively substitute all expressions from `substitutions` into `expr`"""
        
        from sympy import S,Function,preorder_traversal,postorder_traversal
                
        if ignore == None:
            ignore = set()
        
        if expr in ignore:
            return expr
        
        functions = {type(F):(F,substitutions[F]) for F in substitutions if isinstance(F,Function)}
        expr = S(expr)
    
        stop = False
        
        final_subs = []
        
        ignore.add(expr)
        for arg in preorder_traversal(expr):
            #print 'test: %s' % arg
            rep = None
            if arg not in final_subs:
                #print 'test 2: %s' % arg
                if arg in substitutions:
                    rep = substitutions[arg]
                elif isinstance(arg,Function) and type(arg) in functions:
                    #print 'test 3: %s' % arg
                    F,rep = functions[type(arg)]
                    f_subs = {F.args[i]:arg.args[i] for i in range(len(arg.args)) if not F.args[i]==arg.args[i]}
                    if debug==True: print "substituting %s in\t%s" % (f_subs,S(rep))
                    rep = S(rep).subs(f_subs)
                if rep != None:
                    if debug==True: print "substituting %s with %s in\t%s" % (arg,rep,expr)
                    final_subs.append( [arg,Settings.substitute_all(rep,substitutions,ignore,debug)] )
        
        res = Settings.substitute_all(expr.subs(final_subs),substitutions,ignore,debug)
        
        ignore.remove(expr)
                
        if len(ignore) == 0 and simplify_timeout > 0:
            try:
                @timeout(simplify_timeout)
                def simplify_expr(expr):
                    from sympy import simplify
                    return simplify(expr)
                expr = simplify_expr(res)
            except TimeoutError:
                expr = res
                pass
        else:
            expr = res
                
        return expr
    
    def _get_substitutions_without_units(self,substitute_numeric = True):
        res = self.dictionary()
        for key in res.keys():
            if res[key]==key:
                del res[key]
        for symbol in self.base_units.keys():
            if symbol in res:
                del res[symbol]
        if not substitute_numeric:
            for key in self.numerics.keys():
                if key in res:
                    del res[key]
        return res
    
    def _get_substitutions_with_safe_units(self):
        res =  self._get_substitutions_without_units()
        tmp_dict = res.copy()
        for symbol in self.base_units.keys():
            res[symbol] = self.substitute_all(self.data[symbol],tmp_dict)
        return res
    
    def get(self,expr,type = None,substitute_numeric = False,substitute_units = False,**kwargs):
        """Returns the expression substituted with all known substitutions. If the optional parameter `substitute_numeric` is True, also substitutions from the 'numerics' category will be used. If the parameter `substitute_units` is True, substitutions from the category `base_units` will also be used. If a `type` is provided as second argument, substitutions from all categories will be used and the result will be casted to `type`."""
        
        if isinstance(expr,(list,tuple)):
            res = []
            for e in expr:
                res.append(self.get(e,type,substitute_numeric,substitute_units,**kwargs))
            if isinstance(expr,tuple):
                res = tuple(res)
            return res
        
        if type != None or substitute_units:
            main_dict = self._get_substitutions_with_safe_units()
        else:
            main_dict = self._get_substitutions_without_units(substitute_numeric=substitute_numeric)
        
        expr = self.substitute_all(expr,main_dict,**kwargs)
        
        if type != None: 
            try:
                return type(expr)
            except TypeError:
                raise TypeError('cannot convert expression to float: %s' % expr)
        else: return expr
    
    def get_numeric(self,expr,**kwargs):
        """Same as get(expr,substitute_numeric=True,**kwargs)"""
        return self.get(expr,substitute_numeric=True,**kwargs)
    
    def get_unitless(self,expr,**kwargs):
        """Same as get(expr,substitute_numeric=True,substitute_units=True,**kwargs)"""
        return self.get(expr,substitute_numeric=True,substitute_units=True,**kwargs)
    
    def set_physical_size(self,sx,sy,sz):
        'Sets the physical box size of the simulation in x, y and z direction'
        s = self.simulation_box
        s.unlock('xmin')
        s.unlock('xmax')
        s.unlock('ymin')
        s.unlock('ymax')
        s.unlock('zmin')
        s.unlock('zmax')
        s.unlock('sx')
        s.unlock('sy')
        s.unlock('sz')
        
        s.sx,s.sy,s.sz = (sx,sy,sz)
        s.xmin = -s.sx/2
        s.ymin = -s.sy/2
        s.zmin = 0
        s.xmax = s.sx/2
        s.ymax = s.sy/2
        s.zmax = s.sz
        
        s.lock('xmin','defined by sx')
        s.lock('xmax','defined by sx')
        s.lock('ymin','defined by sy')
        s.lock('ymax','defined by sy')
        s.lock('zmin','defined by sz')
        s.lock('zmax','defined by sz')
        
    def set_voxel_size(self,nx,ny,nz):
        'Sets the voxe size of the simulation in x, y and z direction'
        voxel_size = (nx,ny,nz)
        self.simulation_box.nx,self.simulation_box.ny,self.simulation_box.nz = voxel_size
        
    def set_simulation_box(self,physical_size,voxel_size):
        """Sets the simulation box size using the physical_size and voxel_size arguments which are 3-tupels containing the simulation box dimensions in x, y, and z directions."""
        self.set_physical_size(*physical_size)
        self.set_voxel_size(*voxel_size)
