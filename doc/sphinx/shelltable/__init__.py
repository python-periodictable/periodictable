import periodictable.core

# Delayed loading of the element discoverer information
def _load():
    """
    The name of the person or group who discovered the element.
    """
    from . import shelltable
    shelltable.init(periodictable.core.default_table())
periodictable.core.delayed_load(
    ['shells'], _load, isotope=True, element=False)
