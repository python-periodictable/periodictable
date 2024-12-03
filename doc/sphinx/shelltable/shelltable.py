"""
Example of isotope specific extensions to the periodic table.
"""
from periodictable.core import Isotope

def init(table, reload=False):
    if 'shells' in table.properties and not reload: return
    table.properties.append('shells')

    # Set the default.  This is required, even if it is only
    # setting it to None.  If the attribute is missing then
    # the isotope data reverts to the element to supply the
    # value, which is almost certainly not what you want.
    Isotope.shells = None

    # Load the data
    for symbol,eldata in data.items():
        el = table.symbol(symbol)
        for iso,isodata in eldata.items():
            el[iso].shells = isodata

# Define the data
data = dict(
    Fe = {56: "56-Fe shell info",
          58: "58-Fe shell info",
         },
    )
