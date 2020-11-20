# Doxygen troubleshooting

# Travis

The VM requires latex to provide infrastructure to compile the bibliographic references for use in
generated html and pdf documents.

We currently support the use of math and formulas in figure captions.
The will look similar to:

```
_Discrete_PG.dox:39: warning: Illegal command \f as part of a image
```

Backslash (escape) f is ok, other flagged illegal commands that are latex or math based need to be
double backslashed.  See below.

```
_Discrete_PG.dox:39: warning: Illegal command \Phi as part of a image
```

The above is produced if this was used: `\f$\Phi\f$`

The fix is to put a double backslash Phi to: `\f$\\Phi\f$`

For latex, just provide the regular math equation without `\f`.  For latex, just use `$\Phi$`.

Full examples for html and latex:
```
\image html shao3.png "Diagram of sublayer thickness for the sublayer bounded by surfaces \f$\\gamma_n\f$ and \f$\\gamma_{n+1}\f$."
\imagelatex{shao3.png,Diagram of sublayer thickness for the sublayer bounded by surfaces $\gamma_n$ and $\gamma_{n+1}$.,\includegraphics[width=\textwidth\,height=\textheight/2\,keepaspectratio=true]}
```

File, shao3.png, has to be added to the doxygen configuration file using `LATEX_EXTRA_FILES`.

# Warnings

## anchor

Using the `\anchor` prior to a paragraph has the side effect of turning the paragraph into a
block paragraph.

## arguments to functions or subroutines

This snippet will not properly document the MEKE argument.
```
!> Integrates forward-in-time the MEKE eddy energy equation.
!! See \ref section_MEKE_equations
subroutine step_forward_MEKE(MEKE, h, visc, dt, G, CS)
type(MEKE_type),                       pointer       :: MEKE !< MEKE
```

We also discovered that a warning is generated if unknown is used or the same word is
repeated for multiple arguments.

```
type(MEKE_type),                       pointer       :: MEKE !< unknown
```

**This is better**
```
type(MEKE_type),                       pointer       :: MEKE !< MEKE data
```

It is best to provide as much text about the argument as possible.

# Workarounds

## Footnotes

*Syntax:* `\footnote{text}`

The footnote command is split between latex and html.  For latex, it is passed through unchanged and latex generates footnotes appropriately.  For html(doxygen), a superscript `[*]` is created with a title attrbute equalling the text.  For html(sphinx), the footnote is translated into footnote RST syntax.  Footnotes are automatically numbered and numbering is restarted for individual pages.

NOTE: Currently the footnote text is not processed.  Any embedded citations currently do not function properly for html(doxygen).  They were recently fixed in html(sphinx) and pdf(sphinx).
