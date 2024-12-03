import periodictable.core

# Delayed loading of the element discoverer information
def _load_discoverer():
    """
    The name of the person or group who discovered the element.
    """
    from . import discoverer
    discoverer.init(periodictable.core.default_table())
periodictable.core.delayed_load(['discoverer'], _load_discoverer)
