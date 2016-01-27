from .categorized_dictionary import Category,CategorizedDictionary

class Settings(CategorizedDictionary):
    """
    Stores the parameters for a simulation. This class inherits from CategorizedDictionary. Initially, it contains four categories named finite_differences, simulation_box, numerics and base_units. See the info(name) member function for more information about the categories and keys.
    There is an additional symbols member containing aliases to all keys. Adding a key with the same name will overwrite the name in the symbols member.
    The `get` function provides a method of transforming an expression using all known substitutions defined by the key/value pairs of this dictionary.
    """
    
    def __init__(self,create_categories = True):
        import collections

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
        
    def copy(self):
        "Return a soft copy of the Settings object."
        copy = Settings(create_categories = False)
        super(Settings,self).copy(copy = copy)
        copy.initializers = (self.initializers.copy())
        copy.updaters = (self.updaters.copy())
        return copy

    def _is_numeric(self,value):
        "Tests if a value should be classified as numeric"
        if isinstance(value,(bool,int,long,float,complex)):
            return True
        try:
            from pycas import S
            S(value).N()
            return True
        except:
            pass
        return False

    def _set_value(self,key,value):
        super(Settings,self)._set_value(key,value)
        
        is_numeric = self._is_numeric(value)
        
        if key not in self.numerics.keys() and is_numeric:
            self.numerics.add_key(None,key)
        elif key in self.numerics.keys() and not is_numeric:
            self.numerics.remove_key(key)

    def get(self,expr,numeric = False):
        from pycas import Expression,RewriteEvaluator,ReplaceEvaluator,MultiEvaluator,Wildcard,S
        import units

        replacement_evaluator = ReplaceEvaluator(recursive=True)
        rule_evaluator = RewriteEvaluator(recursive=True)
        evaluator = MultiEvaluator(recursive=True)

        evaluator.add_evaluator(replacement_evaluator)
        evaluator.add_evaluator(rule_evaluator)

        numeric_keys = self.numerics.keys()

        for s,r in self.data.iteritems():

            if not numeric and (s in numeric_keys):
                continue

            if isinstance(s,Expression) and isinstance(r,(Expression,int,float,complex)):
                sr = S(r)

                if s==sr or (not numeric and units.contains_unit(sr)):
                    continue

                if s.is_function:
                    wc_args = {arg:Wildcard(arg.name) for arg in s.args}
                    rule_evaluator.add_rule(s.subs(wc_args),sr.subs(wc_args))
                else:
                    replacement_evaluator.add_replacement(s,r)

        return evaluator(expr).evaluate()

    def get_numeric(self,expr):
        return self.get(expr,numeric=True)



















