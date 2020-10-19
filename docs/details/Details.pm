# Details

This documentation will show examples of how the documentation pipeline operates.  We shorten the references to
doxygen to DOX and sphinx generated intermediate resturctured text to RST.

## Fortran fortran and subroutine argument documentation

## Images

For images, the caption must be on a continuous line and double quoted.

*doxygen*

A DOX file might contain:
```
\anchor remap3
\image html remapping3.png "The final state after remapping."
```

In DOX, to create a custom caption link for a reference:
```
\ref remap3 "custom caption text"
```

NOTE: Default caption text is replaced by a internally tracked reference number for latex generated documents.

*sphinx*

The resulting RST:
```
.. _remap3:

.. figure:: /images/remapping3.png

   The final state after remapping.
```

The link can now be used as a reference (`:ref:`).  If the above reference is used, by default is to show the caption text.

To utilize the link and default caption text:

```
As shown in figure :ref:`remap3`.
```

If your writing RST, to supply your own caption text for the same link:
```
As shown in figure :ref:`custom caption text<remap3>`.
```

