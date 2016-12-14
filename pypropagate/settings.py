from .categorized_dictionary import Category,CategorizedDictionary
import sys


class Settings(CategorizedDictionary):
    """
    Stores the parameters for a simulation. This class inherits from CategorizedDictionary. Initially, it contains four categories named finite_differences, simulation_box, numerics and base_units. See the info(name) member function for more information about the categories and keys.
    There is an additional symbols member containing aliases to all keys. Adding a key with the same name will overwrite the name in the symbols member.
    The `get` function provides a method of transforming an expression using all known substitutions defined by the key/value pairs of this dictionary.
    """

    def __init__(self,create_categories = True):
        import collections

        self._cache = dict()
        self._initialized = True
        self._initializing = False
        self._initializers = collections.OrderedDict()
        self.updaters = collections.OrderedDict()

        super(Settings,self).__init__()
        
        self.__doc__ = "Settings of a FiniteDifference simulation"
        
        if create_categories :
            self.create_category("numerics",info="pure numerical values")
            self.create_category("unitless",info="conversions to unitless coordinates")
            self.create_category("symbols",info="symbol aliases to frequently used symbols")

    def create_category(self,cat_name,*args,**kwargs):
        short_name = kwargs.pop('short_name',cat_name)

        import expresso.pycas as pc

        def add_symbol_creation_to_category(cat,short_name,cat_path=[]):
            cat_path = cat_path + [short_name]
            def create_symbol(name,value=None,info=None,**kwargs):
                prefix = '_'.join(cat_path)
                return cat.create_key(name,pc.Symbol("%s_%s" % (name,prefix),**kwargs),value,info)
            def create_function(name,args,value=None,info=None,**kwargs):
                prefix = '_'.join(cat_path)
                return cat.create_key(name,pc.Function("%s_%s" % (name,prefix),**kwargs)(*args),value,info)
            old_create_category = cat.create_category
            def create_category(cat_name,*args,**kwargs):
                short_name = kwargs.pop('short_name',cat_name)
                inner_cat = old_create_category(cat_name,*args,**kwargs)
                add_symbol_creation_to_category(inner_cat,short_name,cat_path)
                return inner_cat

            cat.add_attribute('create_symbol', create_symbol)
            cat.add_attribute('create_function', create_function)
            cat.add_attribute('create_category', create_category)

        cat = super(Settings, self).create_category(cat_name,*args,**kwargs)
        add_symbol_creation_to_category(cat,short_name)

        return cat

    @property
    def initializers(self):

        class Initializers(object):
            def __init__(self,parent):
                self.parent = parent
            def __getitem__(self, item):
                return self.parent._initializers[item]
            def __setitem__(self, key, value):
                self.parent._initialized = self.parent._initializing
                self.parent._initializers.__setitem__(key,value)

        return Initializers(self)

    def initialize(self):
        if self._initialized or self._initializing:
            return

        self._initializing = True

        try:
            for initializer in self._initializers.values():
                initializer(self)
        except:
            self._initialized = False
            self._initializing = False
            raise

        self._initialized = True
        self._initializing = False

    def copy(self,copy_initializers = True,copy_updaters = True):
        "Return a soft copy of the Settings object."
        copy = Settings(create_categories = False)
        super(Settings,self).copy(copy = copy)
        if copy_initializers:
            copy._initializers = (self._initializers.copy())
        if copy_updaters:
            copy.updaters = (self.updaters.copy())
        return copy

    def _is_numeric(self,value):
        "Tests if a value should be classified as numeric"
        if isinstance(value,(bool,int,long,float,complex)):
            return value not in [1,0,-1,1j,-1j]

        import expresso.pycas as pc
        if isinstance(value,pc.Expression):
            try:
                value.N()
                return True
            except:
                pass

            import units
            if units.contains_unit(value):
                return True

        return False

    def _set_value(self,key,value):
        import expresso.pycas

        if isinstance(value,expresso.pycas.Expression):
            value = value.evaluate(cache = self.get_cache())

        if key in self.data and self.data[key] == value:
            return

        super(Settings,self)._set_value(key,value)
        is_numeric = self._is_numeric(value)
        
        if key not in self.numerics.keys() and is_numeric:
            self.numerics.add_key(None,key)
        elif key in self.numerics.keys() and not is_numeric:
            self.numerics.remove_key(key)

        self.clear_cache()
        self._initialized = self._initializing

    def _get_evaluator(self,numeric = False,unitless = False):
        key = (numeric,unitless)

        try:
            return self._cache[key]
        except:
            pass

        from expresso.pycas import Expression,RewriteEvaluator,ReplaceEvaluator,MultiEvaluator,Wildcard,S

        replacement_evaluator = ReplaceEvaluator(recursive=True)
        rule_evaluator = RewriteEvaluator(recursive=True)
        evaluator = MultiEvaluator(recursive=True)

        evaluator.add_evaluator(replacement_evaluator)
        evaluator.add_evaluator(rule_evaluator)

        numeric_keys = self.numerics.keys()
        unitless_keys = self.unitless.keys()

        for s,r in self.data.iteritems():

            if (not numeric and (s in numeric_keys)) or (not unitless and (s in unitless_keys)):
                continue

            if isinstance(s,Expression) and isinstance(r,(Expression,int,float,complex)):
                sr = S(r)
                if s==sr:
                    continue
                if s.is_function and not s.function.is_operator:
                    wc_args = {arg:Wildcard(arg.name) for arg in s.args}
                    rule_evaluator.add_rule(s.subs(wc_args),sr.subs(wc_args))
                else:
                    replacement_evaluator.add_replacement(s,r)

        self._cache[key] = evaluator

        return evaluator

    def clear_cache(self,*args):
        self._cache = {}

    def get_cache(self, *args):
        return self._cache

    def get(self,expr,numeric = False,unitless = False,evaluate = True):
       
        self.initialize()
        evaluator = self._get_evaluator(numeric, unitless)
        res = evaluator(expr, cache = self._cache)
        if evaluate == True:
            res = res.evaluate(cache = self.get_cache()) 
        return res

    def get_numeric(self,expr,**kwargs):
        return self.get(expr,numeric=True,**kwargs)

    def get_unitless(self,expr,**kwargs):
        return self.get(expr,numeric=True,unitless=True,**kwargs)

    def get_coordinates(self,vars):
        if not isinstance(vars,(tuple,list)):
            pass

    def get_optimized(self,expr,**kwargs):
        from expresso.pycas.evaluators.optimizers import optimize_for_compilation
        return optimize_for_compilation(self.get_unitless(expr,**kwargs),cache = self.get_cache())

    def get_definition(self,expr):
        self.initialize()
        return self.data[expr]

    def get_as(self,expr,cast):
        import expresso.pycas as pc
        res = self.get_unitless(expr)
        if res.function == pc.Tuple:
            return [cast(e) for e in res]
        return cast(res)

    def get_numeric_transform(self):
        copy = self.copy()
        return lambda x: copy.get_numeric(x)

