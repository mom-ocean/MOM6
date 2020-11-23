# Latex

What we know about the construction of latex documents based on provided document files.

## Doxygen

The section order of the generated PDF seems to honor the file order provided by the `INPUT=` tag in
the configuration file.

## Sphinx

The generated PDF will follow the general structure of the `index.rst` file.  

For sections of the `index.rst` that link to auto generated pages, the order may appear random.

It turns out that the order of auto generated pages is alphabetical by tag name, `{tag}`, as provided in DOX files.

Here is the doxygen command:
`\page {tag} {name}`
