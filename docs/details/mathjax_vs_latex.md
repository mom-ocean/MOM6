# MathJax vs Latex

MathJax tends to me more forgiving about syntax errors than latex.   This guide will cover
several common issues.  The documentation infrastructure is in continuous development, so
support or workarounds may become available in the future.

Good locations to test equations for both latex and MathJax:
- [LaTex Base](https://latexbase.com/)
- [MathJax](https://www.mathjax.org/#demo)

## Unsupported syntax

The `\mathbold` command is not supported.

Red highlighting is currently unsupported.

## eqnarray

Hints on the use of `\eqnarray`:
- Use of `\mbox{}` may require surrounding braces as in {\mbox{}}
- If a formula needs formatting using `&` you must use eqnarray
- MathJax does not handle backslashes (`\`) within `\mbox{}`
  - Wrong (ok in latex): `\mbox{nonpen\_SW}`
  - Expand the string by moving the backslash out between two mbox commands
  - Correct: `\mbox{nonpen}\_\mbox{SW}`

MathJax automatically switches to text mode within an mbox{}.  If you need to
embed a formula within an mbox{} escape it with a `$`.
