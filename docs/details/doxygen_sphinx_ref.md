# Doxygen vs Sphinix References and Labels

## Reference Labels

Reference labels are not case sensitive in Sphinx.

The following characters are translated to dashes (`-`): colon(`:`) underscore(`_`)
Here are examples:
* 1.8.13: `\eqref{eq:Coriolis_abcd}` => `equation-eq-corilolis-abcd`

## Equations

*Syntax:* `\eqref{tag}`

A [custom command](https://www.doxygen.nl/manual/custcmd.html) is created to process the tag for use in latex and html.  For latex, `\eqref` is intercepted and changed to `\ref`.  For html, it is unchanged and intercepted by [MathJax](#mathjax).

**Why?** For latex, all references to equations use the `\ref` tag.  For MathJax, once the html page is loaded, the javascript hunts for `\eqref` tags and replaces them with links to equations if the autonumbering feature is turned on.  For the MOM6 manual, this feature is turned off and renumbering is performed by an external program.

### Complicated Equation References

These refer to multi-line named equations that need to be referred to from other html pages.  If you do not need to link to these
equations from other html pages, then you do not need to use this custom command.  You can use the `\eqref{tag}` instead.

*Syntax:* `\eqref{tagA,tagB,tagC}`

Large formulas with column formatting using `&` generally are wrapped in a `\f{eqnarray}`:
```
\f{eqnarray}
\label{html:ale-equations}\notag \\
h^\dagger &= h^{(n)} - \Delta t \left[ \nabla_r \cdot \left( h \, \mathbf{u} \right) \right]
&\mbox{thickness} \label{eq:ale-thickness-equation} \\
\theta^\dagger \, h^\dagger &= \theta^{(n)} \, h^{(n)} - \Delta t \left[ \nabla_r \cdot \left( \theta h \, \mathbf{u} \right) - h \boldsymbol{\mathcal{N}}_\theta^\gamma + \delta_r J_\theta^{(z)} \right]
&\;\;\;\;\mbox{potential temp} \label{eq:ale-temperature-equation} \\
h^{(n+1)} &= h^\dagger - \Delta t \, \delta_r \left( z_r \dot{r} \right)
&\mbox{move grid} \label{eq:ale-new-grid} \\
\theta^{(n+1)} h^{(n+1)} &= \theta^\dagger h^\dagger - \Delta t \, \delta_r \left( z_r \dot{r} \, \theta^\dagger  \right)
&\mbox{remap temperature.} \label{eq:ale-remap-temperature}
\f}
```

A *special* label is added right after the `\f{eqnarray}` opening.  This helps produce an implicit label and allows the formula to be linked to from external pages.  Notice that each line of the formula has its own `\mbox` and `\label`.   The combination of the html link and the mbox and label tags form the arguments to the `\eqref` command.

To create a reference link to the first line of the formula, you would use:
`\eqref{eq:ale-thickness-equation,ale-equations,thickness}`.  For latex, this is translated into `\ref{eq:ale-thickness-equation}` and all the references work out as usual.  For html, the second argument should match the special label placed after the opening `\f{eqnarray}`.  The allows generation of a html link back to the entire formula block.  The last argument informs the reader which line of the formula is of interest.

**Why?** MathJax is unable to maintain a list of formula references that span multiple html pages. Restructured text also only supports one label per large `:math:` block.  Restructured text cannot uniquely number equations across pages.

**Sphinx**: NOTE: In the example above, the reference `ale-equations` is translated by sphinx into `equation-ale-equations`.

NOTE: A post-processor, `postProcessEquations.py` has been written to renumber equations for sphinx and doxygen generated html.


## Image Captions and References

For images, the caption must be on a continuous line and double quoted.  The double quotes are stripped
prior to generation of html and pdf.  See doxygen [image](https://www.doxygen.nl/manual/commands.html#cmdimage) command for other options for controlling width and height of the image.

*doxygen*

A DOX file might contain:
```
\anchor remap3
\image html remapping3.png "The final state after remapping."
```

In DOX, to create a custom caption link for a reference:
```
\ref remap3 "Reference text"
```

Using `\ref remap3` will show up as `remap3` in rendered html and pdf.  In Sphinx, the rendered html and pdf will utilize the caption which will be "The final state after remapping."

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

If your writing RST, to supply your own caption text or override the default text for the reference:
```
As shown in figure :ref:`custom caption text<remap3>`.
```

## Tables

To reference tables, be sure the `id=` argument is in quotes.  Example:

```
<caption id="scale_factors">
```

A reference to this table is now possible using `\ref scale_factors`.
