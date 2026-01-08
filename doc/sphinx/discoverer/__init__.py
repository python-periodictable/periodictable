from periodictable.core import default_table, delayed_load

# Delayed loading of the element discoverer information
def _load_discoverer():
    """
    The name of the person or group who discovered the element.
    """
    from . import discoverer
    discoverer.init(default_table())
delayed_load(['discoverer'], _load_discoverer)
