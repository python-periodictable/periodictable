.. _extending:

*****************
Adding properties
*****************

The periodic table is extensible.  Third party packages can
add attributes to the table, and they will appear in all of
the elements.

In order to add a new property to the table, you need to define
a python package which contains the required information, and can
attach the information to the periodic table so that it is
available on demand.  This is done with the function ``init(table)`` in
your table extension.

This example adds the attribute ``discoverer`` to each element.  First
create the file ``discoverer/discoverer.py``:

.. literalinclude:: /discoverer/discoverer.py
   :language: python

Now that we have defined the ``init(table)`` function, we need a way to call it.
The simplest solution is to load it directly when your package is imported.
In the current example, this could be done by adding the following
line to the end of the file::

    init(periodictable.core.default_table())

This would be fine for the current example because the table size is
small and load time is fast.  For large tables, you may wish to
delay loading the table data until it is needed.  To do this, we
use the :func:`delayed_load <periodictable.core.delayed_load>` function in our
package init file ``discoverer/__init__.py``:

.. literalinclude:: /discoverer/__init__.py
   :language: python

The first argument to delayed_load is the list of all attributes that will
be defined when the module is loaded.  The second argument is the loading
function, whose docstring will appear as the attribute description for
each attribute in the first list.

Check that it works:

.. doctest::

    >>> import discoverer
    >>> import periodictable
    >>> print(periodictable.magnesium.discoverer)
    J. Black


Isotope and ion specific data is also supported.  In this case we
need a data table that contains information for each isotope of
each element.  The following example uses a dictionary of elements,
with a dictionary of isotopes for each.  It adds the ``shells``
attribute to Fe[56] and Fe[58].

Define ``shelltable/shelltable.py``:

.. literalinclude:: /shelltable/shelltable.py
   :language: python

Again, we are going to initialize the table with delayed loading.
In this case it is very important that we set the ``isotope=True``
keyword in the ``delayed_load`` call.  If we don't, then the magic we
use to return the correct value after loading the new table information
fails.  Since unknown attributes are delegated to the underlying
element, the value for the natural abundance will be returned
instead.  On subsequent calls the isotope specific value will be
returned.

This is demonstrated in ``shelltable/__init__.py``:

.. literalinclude:: /shelltable/__init__.py
   :language: python

Check that it works:

.. doctest::

    >>> import shelltable
    >>> import periodictable
    >>> print(periodictable.Fe[56].shells)
    56-Fe shell info
    >>> print(periodictable.Ni[58].shells)
    None


Ion specific data is more complicated, particularly because of the
interaction with isotopes.  For example, ``Ni[58].ion[3]`` should have
the same neutron scattering factors as ``Ni[58]`` (the neutron is
only sensitive to the nucleus), but different scattering factors
from ``Ni.ion[3]``.  X-rays are sensitive to the electronic structure
of the atom and not the nucleus, so ``Ni[58].ion[3].xray.f0(Q)``
and ``Ni.ion[3].xray.f0(Q)`` are the same but different
from ``Ni.xray.f0(Q)``.

Current support for ion dependent properties is awkward.  The X-ray
table :mod:`periodictable.xsf` creates a specialized structure
for each ion as it is requested.  The magnetic form factor table
:mod:`periodictable.magnetic_ff` does not try to support ``ion.magnetic_ff``
directly, but instead the user must request ``ion.magnetic_ff[ion.charge]``.
Support for ion mass, which is isotope mass adjusted for the number of
electrons is built into the Ion class.  There are not yet any examples of
extension tables that depend on both isotope and ion.

Be sure to use the ``ion=True`` keyword for ``delayed_load`` when the
table extension contains ion specific information.

