"""
Partial table of element discoverers.

From http://en.wikipedia.org/wiki/Discoveries_of_the_chemical_elements.
"""

import periodictable.core

def init(table, reload=False):
    if 'discoverer' in table.properties and not reload: return
    table.properties.append('discoverer')

    # Set the default, if any
    periodictable.core.Element.discoverer = "Unknown"

    # Not numeric, so no discoverer_units

    # Load the data
    for name,person in data.items():
        el = table.name(name)
        el.discoverer = person

data = dict(
 arsenic="Jabir ibn Hayyan",
 antimony="Jabir ibn Hayyan",
 bismuth="Jabir ibn Hayyan",
 phosphorus="H. Brand",
 cobalt="G. Brandt",
 platinum="A. de Ulloa",
 nickel="A.F. Cronstedt",
 magnesium="J. Black",
)
