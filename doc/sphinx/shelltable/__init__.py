from periodictable.core import default_table, delayed_load

# Delayed loading of the element discoverer information
def _load():
    """
    The name of the person or group who discovered the element.
    """
    from . import shelltable
    shelltable.init(default_table())
delayed_load(
    ['shells'], _load, isotope=True, element=False)
