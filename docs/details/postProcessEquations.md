# Post Processing Equations

Several tasks are performed with the postProcessingEquations.py python script.  The majority of
the work is to post process equations within the HTML.  Additional details are provided below.

NOTE: Some of this processing may be handled by internals of sphinx!

# Doxygen

## Equations with image captions

Here is an example of how to enable equations in an image caption.  You can use the normal `\image latex`
if you are not using equations in the caption.  Double escapes (`\\`) are needed for any latex command in the
`\image html` caption.  e.g. `\\theta`

```
\image html PG_loop.png "Schematic of the finite volume used for integrating the \f$u\f$-component of momentum. The thermodynamic variables \f$\\theta\f$ and \f$s\f$ reside on the sides of the depicted volume and are considered uniform for the vertical extent of the volume but with linear variation in the horizontal. The volume is depicted in \f$(x, p)\f$ space so \f$p\f$ is linear around the volume but \f$\\Phi\f$ can vary arbitrarily along the edges."

\imagelatex{PG_loop.png,Schematic of the finite volume used for integrating the $u$-component of momentum. The thermodynamic variables $\theta$ and $s$ reside on the sides of the depicted volume and are considered uniform for the vertical extent of the volume but with linear variation in the horizontal. The volume is depicted in $(x\, p)$ space so $p$ is linear around the volume but $\Phi$ can vary arbitrarily along the edges.,\includegraphics[width=\textwidth\,height=\textheight/2\,keepaspectratio=true]}
```

This is enabled through the use of the doxygen `ALIASES` command.

```
ALIASES += imagelatex{3}="\latexonly\begin{DoxyImage}\n\3{\1}\n\doxyfigcaption{\2}\n\end{DoxyImage}\endlatexonly\xmlonly<image type=\"latex\" name=\"\1\">\2</image>\endxmlonly"
```

For doxygen 1.8.13, we need the `\xmlonly` command to trigger inclusion of the images in sphinx renderings of html and pdf.

Arguments:
* Argument 1: The image file name
* Argument 2: The image caption
* Argument 3: Image resizing control via latex

The second argument contains the caption text.  If the caption text contains commas (`,`), these must be escaped.

Each image also needs to be included in the doxygen configuration file under `LATEX_EXTRA_FILES`.

## Complex equations blocks

For equations with multiple lines and labels, there is a special equation reference that we provide.  In doxygen,
there is a special command with three arguments (see ALIASES `eqref{3}`).

Ex: \eqref{eq:ale-thickness-equation,ale-equations,thickness}
* Argument 1: The normal latex reference for the specific equation
* Argument 2: The equation reference for the entire equation block (html)
* Argument 3: The label for the equation (html)

The first argument is simply passed through to latex for normal processing.  For html, a `\eqref2{}` is
passed through to the html and processed by `postProcessEquations.py`.

# Sphinx
