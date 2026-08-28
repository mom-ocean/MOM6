import re

from . import get_doxygen_root, get_doxygen_id_index


def flatten(xmlnode):
    # <xmlnode>this.text<child0>child0.text</child0>child0.tail...</xmlnode>

    t = ''

    # text of this node
    if xmlnode.text is not None:
        t += xmlnode.text

    # process all children recursively
    for n in xmlnode:
        t += ' '
        t += flatten(n)
        if n.tail is not None:
            t += ' '
            t += n.tail

    return t

def format_xml_paragraph(xmlnode,build_mode,nsOrig=None,verbosity=0):
    """Format an Doxygen XML segment (principally a detaileddescription)
    as a paragraph for inclusion in the rst document

    Parameters
    ----------
    xmlnode

    Returns
    -------
    lines
        A list of lines.
    """
    # Here we are operating on the entire document for the template
    # This helps support \footnotes{}
    if nsOrig is not None:
        xmlParagraphFormatter = _DoxygenXmlParagraphFormatter()
        xmlParagraphFormatter.setNS(nsOrig)
        xmlParagraphFormatter.setVerbosity(verbosity)
        xmlParagraphFormatter.generic_visit(xmlnode,build_mode=build_mode)
        xmlParagraphFormatter.ns['text'] = [l.rstrip() for l in xmlParagraphFormatter.lines]
        return xmlParagraphFormatter.ns
    else:
        # Return processing for typically ns['text'] only
        # Expand to allow setting of options
        xmlParagraphFormatter = _DoxygenXmlParagraphFormatter()
        xmlParagraphFormatter.setVerbosity(verbosity)
        xmlParagraphFormatter.generic_visit(xmlnode,build_mode=build_mode)
        return [l.rstrip() for l in xmlParagraphFormatter.lines]

class _DoxygenXmlParagraphFormatter(object):
    # This class follows the model of the stdlib's ast.NodeVisitor for tree traversal
    # where you dispatch on the element type to a different method for each node
    # during the traverse.

    # It's supposed to handle paragraphs, references, preformatted text (code blocks), and lists.

    def __init__(self):
        self.ns = {}
        self.lines = ['']
        self.continue_line = False
        # We need to track specified math labels and place them prior to the ".. math::" blocks
        self.math_labels = []
        self.build_mode = None
        self.verbosity = 0
        self.indent = -1
        self.options = []

    # new
    def setNS(self, ns):
        self.ns = ns

    def setVerbosity(self, verbosity):
        self.verbosity = verbosity
        if self.verbosity > 0: print("[debug] verbosity = %s" % (self.verbosity))

    def visit_latexonly(self, node):
        if not(self.build_mode in ('latexpdf','latex')):
            return

        text = node.text
        if text == None:
            return

        # Convert \\ref{tag} to :ref:` ` and the sphinx latex processor
        # converts it to a proper label reference.
        ref_match = re.search('\\\\ref{(.*?)}', text)
        if ref_match is not None:
            tag_string = ref_match.groups()[0]
            #val = [' :ref:`%s`' % tag_string]
            val = [':latex:`\\ref{%s}`' % tag_string]
            #self.lines[-1] += ''.join(val)
            self.concat_text(val[0])
            return

        # If we have <image type="latex"> then skip DoxyImage provided material
        if 'skipDoxyImage' in self.options:
            if text.find('DoxyImage') >= 0:
                return

        # At this point, just pass everything through to latex
        self.concat_text(':latex:`%s`' % (text))

        return

    # new
    # Newer versions of doxygen utilize <htmlonly> tag in XML
    # Doxygen 1.8.13 leaves all this in <para> see: para_eqref
    def visit_htmlonly(self, node):
        if self.build_mode != 'html':
            return

        text = node.text
        if text == None:
            return

        # Check for \eqref2{tag,txt} and convert to :ref:`tag`_
        eqref_match = re.search('\\\\eqref2{(.*?)}', text)
        if eqref_match is not None:
            tag_string = eqref_match.groups()[0]
            if tag_string.find(',') >= 0:
                fc = tag_string.find(',')
                val = [':math:numref:`%s` - %s' % (tag_string[0:fc],tag_string[fc+1:])]
            else:
                val = [':math:numref:`%s`' % tag_string]
            #self.lines[-1] += ''.join(val)
            self.concat_text(val[0])
            return

        # This supports \footnotes{}
        if text.find('title=') >= 0:
            text = text.replace('\n',' ')
            title_match = re.search('title="(.*)"', text)
            if title_match:
                title_string = title_match.groups()[0]
                # Recover \cite that have been converted to @cite to :cite:`%s`
                if title_string.find('@cite') >= 0:
                    citeCommand = '@cite ([\w\-\_]+)'
                    m = re.search(citeCommand, title_string)
                    while m:
                        replStr = title_string[m.start():m.end()]
                        newStr = ':cite:`%s`' % (m.groups()[0])
                        title_string = title_string.replace(replStr, newStr)
                        m = re.search(citeCommand, title_string)
                if 'footnotes' in self.ns:
                    self.ns['footnotes'].append(title_string)
                else:
                    self.ns['footnotes'] = [title_string]

                val = ["[#]_"]
                #self.lines[-1] += ''.join(val)
                self.concat_text(val[0])
                return

        # Check for \eqref{ replace with :ref:`tag`_
        # Post processing of equations will place a link into the HTML
        eqref_match = re.search('\\\\eqref{(.*?)}', text)
        if eqref_match is not None:
            tag_string = eqref_match.groups()[0]
            val = [':math:numref:`%s`' % tag_string]
            #self.lines[-1] += ''.join(val)
            self.concat_text(val[0])
            return

        # undefined
        if self.verbosity > 0:
            print("[debug] WARNING: Uncaptured htmlonly string (%s)" % text)

    # new
    # reStructured text only permits one label per math:: block
    def emit_math_labels(self):
        if len(self.math_labels) == 0:
            return

        if self.verbosity > 0: print("[debug] inserting math labels")

        math_block_idx = -1
        for idx in range(len(self.lines)-1,0,-1):
            if self.lines[idx].startswith('.. math::'):
                math_block_idx = idx
                break

        # Add new label right after the math:: block
        if math_block_idx >=0:
            new_lines = self.lines[0:math_block_idx+1]
            new_label = "   :label: %s" % (self.math_labels[0])
            new_lines.append(new_label)
            #new_lines.append('')
            self.blank_line()
            new_lines = new_lines + self.lines[math_block_idx+1:]
            self.lines = new_lines

        self.math_labels = []

    # Add appropriate implicit labels from anchors
    def visit_anchor(self, node):
        if self.verbosity > 0:
            print("[debug] anchor id(%s)" % (node.get('id')))
        citeID = node.get('id')
        if citeID.find('_1CITE') == 0:
            citeID = "citeref_%s" % (citeID)
        implicitLink = '.. _%s:' % (citeID)
        self.lines.append(implicitLink)
        #self.lines.append('')
        self.blank_line()

    # Original
    def visit(self, node):
        method = 'visit_' + node.tag
        if self.verbosity > 0: print("[debug] method=%s" % (method))
        if len(self.math_labels) > 0 and node.tag != 'formula':
          self.emit_math_labels()
        visitor = getattr(self, method, self.generic_visit)
        return visitor(node)

    def generic_visit(self, node, build_mode=None):
        if build_mode:
            if self.verbosity > 2: print("[debug] Setting build mode: %s" % (build_mode))
            self.build_mode = build_mode
            # Perform a scan for htmlonly or latexonly to prevent double processing of
            # references
            if not('scanned' in self.options):
                self.options.append('scanned')
                self.scanNode(node)
        for child in node.getchildren():
            self.visit(child)
        return self

    # Scan the node and set appropriate options
    def scanNode(self, node):
        # NOTE: these XPath expressions must use './/' rather than '//'.
        # In XPath, '//foo' is an abbreviation for /descendant-or-self::node()/foo
        # starting from the *document root*, not from `node`. Because our
        # autodoc_doxygen extension concatenates every doxygen XML file into a
        # single merged tree (see set_doxygen_xml), every `node` passed here is
        # a small subtree (e.g. a single <detaileddescription>) whose owner
        # document is the *entire* MOM6 doxygen output. Using '//' here
        # therefore scans the whole merged tree on every call, which made
        # scanNode the dominant cost of `make html` -- 75% of single-threaded
        # build time at full MOM6 input scale, quadratic in the tree size.
        # Using './/' scans only descendants of the actual node, which is what
        # was intended and makes each call O(local subtree size).
        xp = node.xpath('.//latexonly')
        if len(xp) > 0:
            self.options.append('latexonly')
        xp = node.xpath('.//htmlonly')
        if len(xp) > 0:
            self.options.append('htmlonly')

        if 'latexonly' in self.options:
            xp = node.xpath('.//image[@type="latex"]')
            if len(xp) > 0:
                self.options.append('skipDoxyImage')

    def visit_ref(self, node):
        refid = node.get('refid')
        name_node = None
        ream_name = None
        kind = None

        # O(1) lookup via the lazily-built id -> element index.
        # The previous findall('.//*[@id=X]') on the merged tree was the
        # single largest cost in `make html` under XML_PROGRAMLISTING=YES
        # -- see docs/REMAINING_TASKS.md / profile notes.
        hit = get_doxygen_id_index().get(refid)
        if self.verbosity > 0: print("[debug] refid(%s) kindref(%s) ref(%s)" %
            (refid, node.get('kindref'), hit))
        if hit is not None:
            ref = hit
            kind = ref.get('kind')
            if self.verbosity > 0: print("[debug] ref(%s)" % ref.items())
            if ref.tag == 'memberdef':
                parent = ref.xpath('./ancestor::compounddef/compoundname')[0].text
                name = ref.find('./name').text
                real_name = parent + '::' + name
            elif ref.tag in ('compounddef', 'enumvalue'):
                if kind == 'page':
                    # :ref: works, but requires an explicit tag placed at the top of pages
                    # that generates an INFO message.  FIX LATER.
                    val = [':ref:`%s`' % ref.get('id')]
                    #self.lines[-1] += ''.join(val)
                    self.concat_text(val[0])
                    return
                name_node = ref.find('./name')
                real_name = name_node.text if name_node is not None else ''
            elif ref.tag in ('anchor','sect1','sect2','sect3','sect4'):
                # If _1CITEREF_ this is a doxygen processed citation
                if refid.find('_1CITEREF_') >= 0:
                    citation = refid[18:]
                    val = [':cite:`%s`' % (citation)]
                    #self.lines[-1] += ''.join(val)
                    self.concat_text(val[0])
                    return
                # Capture sectional links

                # Treat the rest of these as general links
                if refid.find('_1') >= 0:
                    reftext = node.text
                    reftext = reftext.strip()
                    refid2 = refid[refid.find('_1')+2:]
                    if reftext != '' and reftext != refid2:
                        if self.verbosity > 0: print("[debug] refid2(%s) reftext(%s)" % (refid2,reftext))
                        val = [':ref:`%s<%s>`' % (reftext,refid)]
                    else:
                        if self.verbosity > 0: print("[debug] refid(%s)" % (refid))
                        val = [':ref:`%s`' % refid]
                    #self.lines[-1] += ''.join(val)
                    self.concat_text(val[0])
                    return
                else:
                    print('[error] Unimplemented anchor tag: %s' % (ref.tag))
                    raise NotImplementedError(ref.tag)
            else:
                print('[error] Unimplemented tag: %s' % (ref.tag))
                raise NotImplementedError(ref.tag)
        else:
            real_name = None


        # Older doxygen support 1.8.13 for citation references
        if node.get('kindref') == 'member' and refid.find('_1CITEREF_') >= 0:
            citation = refid[18:]
            val = [':cite:`%s`' % (citation)]
            #self.lines[-1] += ''.join(val)
            self.concat_text(val[0])
            return

        # if kind='file' treat as file references
        if kind == 'file':
            # for now treat these as text
            # TODO: references to code
            val = ['``%s``' % node.text]
            #self.lines[-1] += ''.join(val)
            self.concat_text(val[0])
            return

        #debug
        code_type = 'f'
        if code_type == 'f':
            val = [':%s:func:`%s' % (code_type, node.text)]
        else:
            val = [':%s:any:`' % code_type, node.text]
        if real_name:
            val.extend((' <%s>`' % (real_name)))
        else:
            val.append('`')
        if node.tail is not None:
            val.append(node.tail)

        if self.verbosity > 0: print("[debug] kind(%s) real_name(%s) node_name(%s)" %
            (kind, real_name, name_node))
        #self.lines[-1] += ''.join(val)
        self.concat_text(''.join(val))

    # add visit_ulink
    def visit_ulink(self, node):
        self.para_text('`%s <%s>`_' % (node.text, node.get('url')))

    # add visit_emphasis
    def visit_emphasis(self, node):
        self.para_text('*%s*' % node.text)

    # add role_text
    def role_text(self, node, role):
        # Is this even used?
        if self.verbosity > 0:
            print("[debug] role_text")
        # XXX we should probably escape preceeding whitespace...
        # but there's no backward equivalent of `tail`
        text = ' :%s:`%s`' % (role, node.text)

        if node.tail is not None and not node.tail.startswith(' '):
            # escape following whitespace
            text += '\\'

        text += ' ' # interpretered text needs surrounding whitespace
        self.para_text(text)

    # add visit_image
    def visit_image(self, node):

        # Filter activity based on build type and type of image
        image_type = node.get('type')
        if image_type == 'html' and self.build_mode != 'html':
            return
        if image_type == 'latex' and not(self.build_mode in ('latexpdf','latex')):
            return

        if self.verbosity > 0:
            print("[debug] image type(%s) mode(%s)" % (image_type, self.build_mode))

        # node.text is None for an empty <image/> element (no caption text);
        # treat that the same as an empty caption and emit `.. image::` rather
        # than `.. figure::`. The fork's original code crashed with
        # AttributeError on these. Doxygen produces empty <image> elements for
        # cases like an image referenced from a `\image` command with no caption.
        if node.text and node.text.strip():
            type = 'figure'
        else:
            type = 'image'

        self.lines.append('.. %s:: /images/%s' % (type, node.get('name')))

        if type in 'figure':
            # NOTE: Escaped math equations do not play nicely with "literal strings" in python!

            caption = node.text

            # Detect simple math commands and replace them with sphinx :math: directives
            mathCommand = '\\\\f\$(.*?)\\\\f\$'
            m = re.search(mathCommand, caption)
            while m:
                replStr = caption[m.start():m.end()]
                newStr = ':math:`%s`' % (m.groups()[0])
                caption = caption.replace(replStr, newStr)
                m = re.search(mathCommand, caption)

            # Only html needs to be double escaped
            if image_type == 'html':
                caption = node.text.replace('\\','\\\\')

            #if caption.find('Phi') >= 0:

            # Search for $[math]$ and convert to \([math]\) for html
            # Use :math: for latex
            mathCommand = '\$(.*?)\$'
            m = re.search(mathCommand, caption)
            while m:
                replStr = caption[m.start():m.end()]
                if image_type == 'html':
                    newStr = '\\\\(%s\\\\)' % (m.groups()[0])
                else:
                    newStr = ':math:`%s`' % (m.groups()[0])
                caption = caption.replace(replStr, newStr)
                m = re.search(mathCommand, caption)

            # For html, scan for \\f and remove that too
            if image_type == 'html':
                caption = caption.replace('\\f','')

            if self.verbosity > 0:
                # For math we have to double the number of escapes so we pass an
                # escape from RST to HTML.
                print("[debug] caption text(%s)" % (caption))
            self.lines.extend(['', "   %s" % (caption), ''])

    # add visit_superscript
    def visit_superscript(self, node):
        self.role_text(node, 'superscript')

    # add visit_subscript
    def visit_subscript(self, node):
        self.role_text(node, 'subscript')

    # add visit_sup
    # Support for doxygen 1.8.13 as it passes everything to <para>
    # Support for doxygen \footnote{}
    def visit_sup(self, node):

        # Skip if we detect htmlonly or latexonly
        if self.para_ignore():
            return

        title_string = node.get('title')
        if title_string:
            citeCommand = '@cite ([\w\-\_]+)'
            m = re.search(citeCommand, title_string)
            while m:
                replStr = title_string[m.start():m.end()]
                newStr = ':cite:`%s`' % (m.groups()[0])
                title_string = title_string.replace(replStr, newStr)
                m = re.search(citeCommand, title_string)

            if 'footnotes' in self.ns:
                self.ns['footnotes'].append(title_string)
            else:
                self.ns['footnotes'] = [title_string]

            val = ["[#]_"]
            #self.lines[-1] += ''.join(val)
            self.concat_text(val[0])

    # Ignore duplicates provided by xmlonly if we detect latexonly or htmlonly
    def para_ignore(self):

        if 'latexonly' in self.options or 'htmlonly' in self.options:
            return True
        return False

    # add replace any references of \eqref, \eqref2, \eqref4
    # Doxygen 1.8.13
    # html: use eqref2; remove eqref4
    # latex: use eqref4; remove eqref2
    # with appropriate replacements
    # Remove duplicates here if latexonly or htmlonly is detected
    def para_eqref(self, text):

        chg = True
        while text.find('\\\\eqref2') >= 0 and chg:
            chg = False
            m = re.search('\\\\eqref2{(.*?)}', text)
            if m:
                ref = m.groups()[0]
                fullRef = '\\\\eqref2{%s}' % (ref)
                if ref.find(',') >= 0:
                    i = ref.find(',')
                    sphinxRef = ':math:numref:`%s` - %s' % (ref[0:i],ref[i+1:])
                    if self.build_mode in ('latexpdf','latex'):
                        sphinxRef = ''
                    if self.para_ignore():
                        sphinxRef = ''
                    text = text.replace(fullRef, sphinxRef)
                    chg = True

        chg = True
        while text.find('\\\\eqref') >= 0 and chg:
            chg = False
            m = re.search('\\\\eqref{(.*?)}', text)
            if m:
                ref = m.groups()[0]
                fullRef = '\\\\eqref{%s}' % (ref)
                if self.build_mode in ('latexpdf','latex'):
                    sphinxRef = ':latex:`\\ref{%s}`' % ref
                else:
                    sphinxRef = ':math:numref:`%s`' % (ref)
                if self.para_ignore():
                    sphinxRef = ''
                text = text.replace(fullRef, sphinxRef)
                chg = True

        chg = True
        while text.find('\\\\eqref4') >= 0 and chg:
            chg = False
            m = re.search('\\\\eqref4{(.*?)}', text)
            if m:
                ref = m.groups()[0]
                fullRef = '\\\\eqref4{%s}' % (ref)
                sphinxRef = ':latex:`\\ref{%s}`' % (ref)
                if self.build_mode in ('html'):
                    sphinxRef = ''
                if self.para_ignore():
                    sphinxRef = ''
                text = text.replace(fullRef, sphinxRef)
                chg = True

        return text

    # Assistant for ensuring there is blank lines between directives
    # It makes sure we do not overly add blank lines
    def blank_line(self):
        if len(self.lines) == 0:
            return

        if self.lines[-1] == '':
            return

        self.lines.append('')

    # Assistant for putting sentences together
    def concat_text(self, text):
        if len(self.lines) == 0:
            self.lines.append(text)
            return

        lastLine = self.lines[-1]

        if len(lastLine) == 0:
            self.lines[-1] = text
            return

        lastChar = lastLine[-1]
        newText = text
        if len(newText) == 0:
            return

        firstChar = newText[0]

        # Emphasis
        if lastChar == "*" or firstChar == "*":
            newText = " %s" % (newText)
            firstChar = " "

        # whitespace after :cite:`tag`
        if lastChar == '`':
            if (firstChar >= 'a' and firstChar <= 'z') or (firstChar >= 'A' and firstChar <= 'Z') or firstChar in ['(','[','{']:
                newText = " %s" % (newText)
                firstChar = " "

        # whitespace before :commands:
        if firstChar == ':':
            if (lastChar >= 'a' and lastChar <= 'z') or (lastChar >= 'A' and lastChar <= 'Z') or lastChar in [',','.','=']:
                newText = " %s" % (newText)
                firstChar = " "

        # Footnotes and any items that end with _
        if newText == '[#]_':
            newText = " %s" % (newText)
        if lastChar == '_':
            if len(lastLine) > 3:
                if lastLine[-4:] == "[#]_" and firstChar != '.':
                    newText = " %s" % (newText)
                    firstChar = " "
                else:
                    newText = " %s" % (newText)

        # Inline text check for space before (``)
        if len(newText) >= 2:
            if newText[0:2] == "``" and lastChar != ' ':
                newText = " %s" % (newText)
                firstChar = " "

        self.lines[-1] += newText
        return

    # add para_text parser
    # Doxygen 1.8.13 support for \eqref \eqref2
    def para_text(self, text):

        if text is not None:
            if text.find('Some time later') >= 0:
                a = 0

            if text.find('eqref') >= 0:
                text = self.para_eqref(text)
            if self.continue_line:
                if len(self.lines) >= 1:
                    # If we are in a continue_line situation but already
                    # have a linefeed, do an append instead
                    if self.lines[-1] == '':
                        self.lines.append(text)
                        return
                self.concat_text(text)
            else:
                self.lines.append(text)

    def visit_para(self, node):

        self.para_text(node.text)

        # visit children and append tail
        for child in node.getchildren():
            self.visit(child)
            self.continue_line = True

            if child.tail is not None:
                self.para_text(child.tail.lstrip())

        # replaced
        #if node.text is not None:
        #    if self.continue_line:
        #        self.lines[-1] += node.text
        #    else:
        #        self.lines.append(node.text)
        #self.generic_visit(node)

        self.continue_line = False
        #self.lines.append('')
        self.blank_line()

    # add visit_formula
    def visit_formula(self, node):
        text = node.text

        # Remove the faked link for pdf version
        if self.build_mode in ('latexpdf','latex'):
            label_match = re.search(' \\\\label{(html:.*?)}.*?\\\\\\\\', text)
            if label_match:
                replace_string = label_match.group()
                text = text.replace(replace_string,'')

        # detect inline or block math
        if text.startswith('\\[') or not text.startswith('$'):
            if text.startswith('\\['):
                text = text[2:-2]

            # if we are emitting a math block and we have
            # pending math labels, go back and emit those
            # first.
            if len(self.math_labels) > 0:
                self.emit_math_labels()

            self.blank_line()
            if '\n' in text:
                self.lines.append('.. math::')
                self.lines.append('')
                for mathline in text.split('\n'):
                    self.lines.append('   ' + mathline)
            else:
                self.lines.append('.. math:: ' + text)
            self.blank_line()
            # Math blocks require an explicit blank line as well?
            #self.lines.append('')
            self.continue_line = False
        else:
            inline = ':math:`' + node.text.strip()[1:-1].strip() + '`'
            if self.continue_line:
                #self.lines[-1] += inline
                self.concat_text(inline)
            else:
                self.lines.append(inline)

            self.continue_line = True

        # detect \label{html:tag} blocks
        if text.find('\\label') >= 0:
            # If we have a big block of equations, supply one label
            label_matches = re.findall('\\\label{html:(.*?)?}',text)
            if len(label_matches) > 0:
                [self.math_labels.append(i) for i in label_matches]
            else:
                label_matches = re.findall('\\\label{(.*?)?}',text)
                if len(label_matches) > 0:
                    [self.math_labels.append(i) for i in label_matches]
                    if self.verbosity > 0:
                        # For math we have to double the number of escapes so we pass an
                        # escape from RST to HTML.
                        print("[debug] math_labels(%s)" % (label_matches))

    def visit_parametername(self, node):
        if 'direction' in node.attrib:
            direction = '[%s] ' % node.get('direction')
        else:
            direction = ''

        param_name = node.text or ''

        # Look up the parameter's Fortran type from the parent
        # memberdef's <param> elements. The <parameterlist> inside
        # <detaileddescription> only carries names and descriptions;
        # the types live on the <param> siblings of <detaileddescription>.
        param_type = ''
        memberdefs = node.xpath('./ancestor::memberdef')
        if memberdefs:
            for p in memberdefs[0].findall('param'):
                defname = p.find('defname')
                if defname is not None and defname.text == param_name:
                    type_el = p.find('type')
                    if type_el is not None:
                        param_type = ''.join(type_el.itertext()).strip()
                    break

        # Prepend the type as literal text in the description rather
        # than using :param type name: (sphinx-fortran's regex can't
        # handle commas in the type) or :type name: (sphinx-fortran's
        # xref resolver crashes on % in dimension expressions).
        if param_type:
            self.lines.append(':param %s: ``%s`` %s' % (param_name, param_type, direction))
        else:
            self.lines.append(':param %s: %s' % (param_name, direction))
        self.continue_line = True

    def visit_parameterlist(self, node):
        lines = [l for l in type(self)().generic_visit(node).lines if l != '']
        # replaced
        #self.lines.extend([':parameters:', ''] + ['* %s' % l for l in lines] + [''])
        self.lines.extend([''] + lines + [''])

    # TODO: Doxygen generates a simplesect for functions with
    # a specified return argument.  For now, we leave as
    # :returns undefined:
    # marker so we can fix up the document using flint.
    # Supports doxygen /sa or /see command
    def visit_simplesect(self, node):
        if self.verbosity > 0:
            print("[debug] simplesect kind(%s)" % (node.get('kind')))

        # Do nothing for \note for now

        # fortran function handling
        if node.get('kind') == 'return':
            self.lines.append(':returns undefined: ')
            self.continue_line = True
            self.generic_visit(node)

        # Add bold text psudo section for \see, \sa roughly acts like doxygen
        if node.get('kind') in ('see', 'sa'):
            see_also_label = "See also"
            #self.lines.append('')
            self.blank_line()
            self.lines.append('**%s**' % (see_also_label))
            #self.lines.append('')
            self.blank_line()
            #self.lines.append('')
            self.generic_visit(node)
    # add

    def visit_sect(self, node, char):
        """Generic visit section"""
        title_node = node.find('title')
        if title_node is not None:
            title = title_node.text
            # Filter html data (possibly if we see a <, / and >)
            if self.verbosity > 0:
                print("[debug] visit_sect id(%s) title(%s)" % (node.get('id'),title))
            if title.find('<') >=0 and title.find('>') >=0 and title.find('/') >=0:
                html_match = False
                # Filter <tt> => ``
                if title.find("<tt>") >= 0:
                    title = title.replace('<tt>','``')
                    title = title.replace('</tt>','`` ')
                    html_match = True
                if not(html_match) and self.verbosity > 0:
                    print("[debug] unmatched html (%s)" % (title))
            # Add a implicit lable for the sections
            implicitLink = '.. _%s:' % (node.get('id'))
            self.lines.append(implicitLink)
            #self.lines.append('')
            self.blank_line()
            self.lines.append(title)
            self.lines.append(len(title) * char)
            #self.lines.append('')
            self.blank_line()

        self.generic_visit(node)

    def visit_sect1(self, node):
        self.visit_sect(node, '=')

    def visit_sect2(self, node):
        self.visit_sect(node, '-')

    def visit_sect3(self, node):
        self.visit_sect(node, '^')

    def visit_sect4(self, node):
        self.visit_sect(node, '"')

    # add end

    # allows us to handle nested ordered lists
    def visit_orderedlist(self, node):
        self.indent = self.indent + 1
        self.generic_visit(node)
        #self.lines.append('')
        self.blank_line()
        self.indent = self.indent - 1

    # allows us to handle nested itemized lists
    def visit_itemizedlist(self, node):
        self.indent = self.indent + 1
        self.generic_visit(node)
        #self.lines.append('')
        self.blank_line()
        self.indent = self.indent - 1

    # Source of citation and numbering
    def visit_listitem(self, node):
        #char = '*' if node.getparent().tag == 'itemizedlist' else '#.'
        if node.getparent().tag == 'itemizedlist':
            #self.lines.append('')
            char = '*'
        else:
            char = '#.'
        if self.verbosity > 1: print("[debug] listitem indent = %s" % (self.indent))
        self.lines.append(' '*(self.indent*2) + char + ' ')
        # replaced
        #self.lines.append('   - ')
        self.continue_line = True
        # recursion
        self.generic_visit(node)

    # add
    def preformat_text(self, lines):
        self.lines.extend(('::', ''))
        self.lines.extend(['  ' + l for l in lines])
        #self.lines.append('')
        self.blank_line()

    def visit_preformatted(self, node):
        segment = [node.text if node.text is not None else '']
        for n in node.getchildren():
            segment.append(n.text)
            if n.tail is not None:
                segment.append(n.tail)

        lines = ''.join(segment).split('\n')
        # add line
        self.preformat_text(lines)
        # extra? no effect
        #self.lines.extend(('.. code-block:: C++', ''))
        #self.lines.extend(['  ' + l for l in lines])

    # add 
    def visit_programlisting(self, node):
        lines = []
        for n in node.getchildren():
            lines.append(flatten(n))
        self.preformat_text(lines)

    #add
    def visit_verbatim(self, node):
        self.visit_preformatted(node)

    def visit_computeroutput(self, node):
        c = node.find('preformatted')
        if c is not None:
            return self.visit_preformatted(c)
        # add
        # I don't think we can put links inside
        # computeroutput text...
        #self.lines[-1] += '``' + flatten(node) + '`` '
        self.concat_text('``' + flatten(node) + '``')
        # omitted
        #return self.visit_preformatted(node)

    def visit_xrefsect(self, node):
        if node.find('xreftitle').text == 'Deprecated':
            sublines = type(self)().generic_visit(node).lines
            self.lines.extend(['.. admonition:: Deprecated'] + ['   ' + s for s in sublines])
            return
        # add - if not depricated
        title = node.find('xreftitle').text
        sublines = type(self)().generic_visit(node).lines
        self.lines.extend(['.. admonition:: %s' % title] + ['   ' + s for s in sublines])
        #else:
        #    raise ValueError(node)

    def visit_subscript(self, node):
        #self.lines[-1] += '\ :sub:`%s` %s' % (node.text, node.tail)
        self.concat_text(':sub:`%s` %s' % (node.text, node.tail))

    def visit_table(self, node):
        # save the number of columns
        cols = int(node.get('cols'))
        table = []
        # save the current output
        lines = self.lines

        # get width of each column
        widths = [0] * cols

        # build up the table contents
        for row_node in node.findall('row'):
            row = []
            for i, entry in enumerate(row_node.getchildren()):
                self.lines = ['']
                self.generic_visit(entry)
                row.append(self.lines)

                # find width of this entry (including leading and trailing space)
                widths[i] = max(widths[i], max([len(line) for line in self.lines]) + 2)

            table.append(row)

        def append_row(row):
            # find number of lines in row
            num_lines = max([len(e) for e in row])
            lines = []

            for k in range(num_lines):
                line = '|'
                for i, e in enumerate(row):
                    if k < len(e):
                        # this is a valid line
                        line += ' ' + e[k]
                        # pad rest of line
                        line += ' ' * (widths[i] - len(e[k]) - 1)
                    else:
                        # invalid line, just fill with spaces
                        line += ' ' * widths[i]

                    line += '|'

                lines.append(line)

            return lines

        self.lines = lines
        # start with a blank
        #self.lines.append('')
        self.blank_line()

        # usual separator line
        sep = '+'
        for width in widths:
            sep += '-' * width
            sep += '+'

        self.lines.append(sep)

        # header row
        self.lines.extend(append_row(table[0]))
        # header separator uses '=' instead of '-'
        self.lines.append(sep.replace('-', '='))

        # loop over body rows
        for row in table[1:]:
            self.lines.extend(append_row(row))
            self.lines.append(sep)

        # end with a blank
        #self.lines.append('')
        self.blank_line()
