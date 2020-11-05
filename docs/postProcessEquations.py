import os, sys, pathlib, re
import itertools
from lxml import html
import argparse

# NOTES:
# &nbsp; = \xa0 = &#160;

class equationRenumber:

    def __init__(self, rootDir, buildType):
        '''
        ref, target: targets that require checking or a transformation
        checkFiles: list of htmlFiles that need to be checked

        eref, etarget: existing references and targets

        Keep a list of files with each reference tag
        eref[#tag] = [htmlFile,htmlFile,....]

        etargets[#tag] = htmlFile location of the target
          There should only be one unique target key amongst
          any number of files.

        fixanchor:
        fixtarget:
        fixcaption:
        '''
        if os.path.isdir(rootDir):
            os.chdir(rootDir)
        self.htmlFiles = []
        self.parsedFiles = []
        self.meta = {
                'ref': {},
                'target': {},
                'targets': [],
                'eref': {},
                'etarget': {},
                'fixdiv': {},
                'checkFiles': [],
                'fixanchor': {},
                'fixtarget': {},
                'fixcaption': {},
                }
        self.rootDir = os.getcwd()
        self.curDir = os.getcwd()
        self.buildType = buildType
        self.breakpoint = False
        self.verbose = False
        self.updates = False

    def getHtmlFiles(self):
        for root, dirs, files in os.walk(self.rootDir):
            for name in files:
                if pathlib.Path(name).suffix == ".html":
                    self.htmlFiles.append(os.path.join(root, name))
        if self.verbose: print("Scanned %s (%d)" % (self.rootDir,len(self.htmlFiles)))

    def walkDoc(self, tree):
        '''This walks the document finding all the <a href=""></a> elements.'''

        hrefs = tree.xpath('//a[@href!=""]')

        # no links found, end of path
        if len(hrefs) == 0:
            return

        for node in hrefs:
            # Make sure we are in the same directory as the tree element to get
            # relative paths correct.
            if self.curDir != os.path.dirname(tree.docinfo.URL):
                self.curDir = os.path.dirname(tree.docinfo.URL)
                os.chdir(self.curDir)
            page = node.get('href')
            if pathlib.Path(page).suffix == ".html":
                #htmlFile = os.path.join(self.rootDir, pathlib.Path(page).resolve())

                htmlFile = pathlib.Path(page).resolve().as_posix()
                #if htmlFile == '/var/www/html/index.html':
                #    import pdb; pdb.set_trace()


                #if self.breakpoint == True:
                #    import pdb; pdb.set_trace()

                # URL for current file/tree being parsed: tree.docinfo.URL
                #import pdb; pdb.set_trace()
                #if self.verbose: print(htmlFile)
                if pathlib.Path(htmlFile).is_file():
                    if not(htmlFile in self.parsedFiles):
                        #print(htmlFile)
                        #if htmlFile == '/var/www/html/MOM6/esmg/_build/html/api/modules.html':
                        #    self.breakpoint = True
                        self.parsedFiles.append(htmlFile)
                        newTree = html.parse(htmlFile)
                        self.walkDoc(newTree)
                #import pdb; pdb.set_trace()

    def htmlWalk(self, start):
        self.htmlFiles.sort()
        #print(self.rootDir,":",len(self.htmlFiles))
        #print(self.htmlFiles[0])

        startFile = os.path.join(self.rootDir, start)
        #if self.verbose:
        #    print("Root:",self.rootDir)
        #    print("Start:",startFile)
        if startFile in self.htmlFiles:
            self.parsedFiles.append(startFile)
            newTree = html.parse(startFile)
            self.walkDoc(newTree)

    def showUnresolved(self):
        '''List unresolved files'''
        for htmlFile in self.htmlFiles:
            if not(htmlFile) in self.parsedFiles:
                print(" unresolved>",htmlFile)

    def collectEquationLabels(self):
        '''
        Label equations sequentially within the tree that is mapped, then
        sequentially label equations in apparently orphaned content.
        '''

        # Fix links/references (doxygen/sphinx)
        for htmlFile in self.parsedFiles:
            tree = html.parse(htmlFile)

            # Search for \eqref and \eqref2 in html (doxygen)
            nodes = tree.xpath("//*[text()]")
            if len(nodes) > 0:
                #if self.verbose:
                #    print(htmlFile,len(nodes))
                ct = 0
                for node in nodes:
                    fullText = "%s%s" % (node.text, node.tail)
                    for m in re.finditer('\\\eqref2{(.*?)}', fullText):
                        if self.verbose:
                            ct = ct + 1
                            if ct == 1:
                                print(os.path.basename(htmlFile),len(nodes))
                            print('  found eqref2>%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                        tag_string = m.groups()[0]
                        fc = tag_string.find(',')
                        if fc >= 0:
                            tag = tag_string[0:fc]
                            if not(htmlFile in self.meta['ref'].keys()):
                                self.meta['ref'][htmlFile] = []
                            if not(tag in self.meta['ref'][htmlFile]):
                                self.meta['ref'][htmlFile].append(tag)
                    for m in re.finditer('\\\eqref{(.*?)}', fullText):
                        if self.verbose:
                            ct = ct + 1
                            if ct == 1:
                                print(os.path.basename(htmlFile),len(nodes))
                            print('  found eqref>%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                        tag = m.groups()[0]
                        # Doxygen 1.8.13
                        tag = tag.lower()
                        tag = tag.translate(tag.maketrans(':_','--'))
                        if not(htmlFile in self.meta['ref'].keys()):
                            self.meta['ref'][htmlFile] = []
                        if not(tag in self.meta['ref'][htmlFile]):
                            self.meta['ref'][htmlFile].append(tag)

            # fixtarget: Search for lone labels in formulas (doxygen)
            # Add <center> node prior to formula as a target
            # <center id="equation-ale-equations" class="math notranslate nohighlight">
            # <span class="eqno">(1)<a class="headerlink" href="#equation-ale-equations"
            # title="Permalink to this equation">&#182;</a></span></center>
            nodes = tree.xpath("//p[@class='formulaDsp']")
            if len(nodes) > 0:
                #if self.verbose:
                #    print(htmlFile,len(nodes))
                ct = 0
                for node in nodes:
                    # skip nodes that already have a <center id='' class='math'> before the <p>
                    prevNode = node.getprevious()
                    if prevNode == None:
                        continue
                    if prevNode.tag == 'center':
                        hasId = prevNode.get('id')
                        hasClass = prevNode.get('class')
                        if hasId and hasClass.find('math') != -1:
                            continue
                        #import pdb; pdb.set_trace()

                    fullText = "%s%s" % (node.text, node.tail)
                    for m in re.finditer('\\\\label{(.*?)}', fullText):
                        # If first label begins with html: then skip
                        tag_string = m.groups()[0]
                        if tag_string.find('html:') == 0:
                            break
                        if self.verbose:
                            ct = ct + 1
                            if ct == 1:
                                print(os.path.basename(htmlFile), len(nodes), 'dox')
                            print('  fixtarget>%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                        tag = tag_string
                        if not(tag in self.meta['fixtarget'].keys()):
                            self.meta['fixtarget'][tag] = []
                        if not(htmlFile in self.meta['fixtarget'][tag]):
                            self.meta['fixtarget'][tag].append(htmlFile)

            #import pdb; pdb.set_trace()

            # fixanchor: Convert lone anchors in front of equations to formal targets
            # TODO: We do not have a full operational example yet.  This is unfinished.
            # TODO: We do not know if we need to change the id.  Assume we need to add equation- prefix.
            # Rule: If the anchor is in a <p>, the next <p> must be the formula
            #       for this to work.
            # Replace lone anchors with a <center> node
            # <center id="equation-ale-equations" class="math notranslate nohighlight">
            # <span class="eqno">(1)<a class="headerlink" href="#equation-ale-equations"
            # title="Permalink to this equation">&#182;</a></span></center>
            # NOTE: Do not fix the tag name here as we will use it later to
            # quickly grab the node we are going to fix.
            nodes = tree.xpath("//p/a[@class='anchor']")
            if len(nodes) > 0:
                ct = 0
                for node in nodes:
                    tag = node.get('id')
                    nextParentNode = node.getparent().getnext()
                    if nextParentNode.tag == 'p' and nextParentNode.get('class') == 'formulaDsp':
                        #if tag.find('equation-') != 0:
                        #    tag = "equation-%s" % (tag)
                        if self.verbose:
                            ct = ct + 1
                            if ct == 1:
                                print(os.path.basename(htmlFile),len(nodes),'dox')
                            print('  fixanchor>%s' % (tag))
                        if not(tag in self.meta['fixanchor'].keys()):
                            self.meta['fixanchor'][tag] = []
                        if not(htmlFile in self.meta['fixanchor'][tag]):
                            self.meta['fixanchor'][tag].append(htmlFile)

            # Search for \label{html: (doxygen/sphinx)
            nodes = tree.xpath("//*[text()]")
            if len(nodes) > 0:
                #if self.verbose:
                #    print(htmlFile,len(nodes))
                ct = 0
                for node in nodes:
                    fullText = "%s%s" % (node.text, node.tail)
                    for m in re.finditer('\\\\label{(html:.*?)}\\\\notag', fullText):
                        if self.verbose:
                            ct = ct + 1
                            if ct == 1:
                                print(htmlFile,len(nodes))
                            print('%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                        tag_string = m.groups()[0]
                        fc = tag_string.find(":")
                        if fc >= 0:
                            tag = tag_string[fc+1:]
                            if not(tag in self.meta['target'].keys()):
                                self.meta['target'][tag] = htmlFile

            # Collect existing equation labels and references (sphinx)
            # These do not need modification other than renumbering

            # Link/Reference
            # <a class="reference internal" href="General_Coordinate.html#equation-h-equations">(2)</a>
            # Some of these references will not match with actual math references due to the use of equation-
            # as part of the tag.
            nodes = tree.xpath("//a[contains(@class,'reference') and contains(@href,'#equation')]")
            if len(nodes) > 0:
                for node in nodes:
                    href = node.get('href')
                    href = href.split('#')
                    if len(href) > 1:
                        tag = href[1]
                        if not(tag in self.meta['eref']):
                            self.meta['eref'][tag] = []
                            if self.verbose:
                                print("check eref>%s (%s)" % (tag,os.path.basename(htmlFile)))
                        if not(htmlFile in self.meta['eref'][tag]):
                            self.meta['eref'][tag].append(htmlFile)
                            if not(htmlFile in self.meta['checkFiles']):
                                self.meta['checkFiles'].append(htmlFile)

            # Anchors (sphinx)
            # <div class="math notranslate nohighlight" id="equation-eq-coriolis-abcd">
            # <span class="eqno">(1)<a class="headerlink" href="#equation-eq-coriolis-abcd"
            # title="Permalink to this equation">¶</a></span>\[**MATH**\]</div>

            # NOTE: Some anchors do not have <span> as above.
            # <div class="math" id="tag">\[**MATH**\]</div>
            # For sphinx, transform to <div><span><a></a></span>\[**MATH**\]</div>
            nodes = tree.xpath("//div[contains(@class,'math')]")
            if len(nodes) > 0:
                ct = 0
                for node in nodes:
                    tag = node.get('id')
                    if tag is None:
                        continue
                    # Check for equation- prefix and that there is a span node with class="eqno"
                    #import pdb; pdb.set_trace()
                    children = node.getchildren()
                    fixDiv = False
                    if len(children) == 0:
                        fixDiv = True
                    else:
                        if children[0].tag != 'span':
                            fixDiv = True
                        else:
                            if children[0].get('class') != 'eqno':
                                fixDiv = True

                    if fixDiv:
                        # If we don't have any children, then we need to fix this div possibly
                        if not(tag in self.meta['fixdiv'].keys()):
                            self.meta['fixdiv'][tag] = []
                        if not(htmlFile in self.meta['fixdiv'][tag]):
                            self.meta['fixdiv'][tag].append(htmlFile)
                        # if we do not have an equation- prefix, add one
                        if tag.find('equation-') != 0:
                            tag = "equation-%s" % (tag)

                    #import pdb; pdb.set_trace()
                    if not(tag in self.meta['targets']):
                        self.meta['targets'].append(tag)
                    if not(tag in self.meta['etarget'].keys()):
                        #self.meta['etarget'][tag] = []
                        ct = ct + 1
                        if self.verbose:
                            if ct == 1:
                                print(os.path.basename(htmlFile))
                            print("target>%s %s" % (tag,fixDiv))
                    if not(tag in self.meta['etarget'].keys()):
                        self.meta['etarget'][tag] = htmlFile
                        if not(htmlFile in self.meta['checkFiles']):
                            self.meta['checkFiles'].append(htmlFile)

            # Search for fixed anchors in case the program is re-run for renumbering
            # Doxygen only
            nodes = tree.xpath("//center[contains(@class,'math')]")
            if len(nodes) > 0:
                for node in nodes:
                    tag = node.get('id')
                    if tag is None:
                        continue
                    if not(tag in self.meta['targets']):
                        self.meta['targets'].append(tag)
                    #if not(tag in self.meta['etarget']):
                    if tag in self.meta['etarget'].keys():
                        print("ERROR: Duplicate target found in %s (%s)" % (htmlFile, tag))
                        sys.exit(1)
                        #self.meta['etarget'][tag] = []
                    else:
                        self.meta['etarget'][tag] = htmlFile
                        if not(htmlFile in self.meta['checkFiles']):
                            self.meta['checkFiles'].append(htmlFile)

            # fix captions: doxygen only
            # Scan figure captions for math and allow them to be shown
            # in MathJax
            # These are fixed on the fly since we are not dealing with targets and
            # references
            # From: \f$h(x)\f$
            # To  : <span class="math notranslate nohighlight">\(h(x)\)</span>
            self.updates = False
            refPattern = {}
            nodes = []
            if self.buildType == 'doxygen':
                nodes = tree.xpath("//div[@class='caption']")
                if len(nodes) == 0:
                    # Older doxygen: <span class="caption-text">
                    nodes = tree.xpath("//span[@class='caption-text']")
            if self.buildType == 'sphinx':
                # Older doxygen: <span class="caption-text">
                nodes = tree.xpath("//span[@class='caption-text']")
                #if len(nodes) == 0: This was XML...
                #    # Even older doxygen: <image type="html" name="Newton_PPM.png">
                #    nodes = tree.xpath("//image[@type='html']")
            if len(nodes) > 0:
                for node in nodes:
                    txt = ""
                    if node.text != None:
                        txt = node.text
                    if node.tail != None:
                        txt = "%s%s" % (txt, node.tail)
                    refPattern[0] = '(\\\\f\$.*?\\\\f\$)'
                    refPattern[1] = '\\\\f\$(.*?)\\\\f\$'
                    m = re.search(refPattern[0],txt)
                    if m:
                        self.fixCaptionMath(refPattern, node)
                    else:
                        refPattern[0] = '(\$.*?\$)'
                        refPattern[1] = '\$(.*?)\$'
                        m = re.search(refPattern[0],txt)
                        if m:
                            self.fixCaptionMath(refPattern, node)

            # Write the tree out if it was modified
            if self.updates:
                # Write tree back out to file
                output = html.tostring(tree)
                #import pdb; pdb.set_trace()
                fo = open(htmlFile, 'wb')
                fo.write(output)
                fo.close()

    def fixCaptionMath(self, refPattern, node):

        doFix = True
        # For the div with the caption, we want to keep scanning
        # until we rewrite all the equations we find.
        #import pdb; pdb.set_trace()
        while doFix:
            doFix = False

            #print("S>",html.tostring(node))
            nodes = [node]
            if len(node.getchildren()) > 0:
                nodes = nodes + node.getchildren()
            #import pdb; pdb.set_trace()
            for x in nodes:
                span = html.Element("span")
                span.set("class","math notranslate nohighlight")
                #print("T>",x.tag)

                #import pdb; pdb.set_trace()
                # Check text
                m = None
                if x.text:
                    m = re.search(refPattern[0],x.text)
                if m:
                    doFix = True
                    self.updates = True
                    txhead = x.text[0:m.start()]
                    txtail = x.text[m.end():]
                    m2 = re.search(refPattern[1],m.groups()[0])
                    #import pdb; pdb.set_trace()
                    span.text = "\\(%s\\)" % (m2.groups()[0])
                    # Do insert
                    x.insert(len(x.getchildren()),span)
                    if len(txhead) > 0:
                        x.text = txhead
                    else:
                        x.text = ''
                    if txtail and len(txtail)>0:
                        span.tail = txtail
                    #import pdb; pdb.set_trace()
                    #aatext = 0

                # Check tail
                m = None
                if x.tail:
                    m = re.search(refPattern[0],x.tail)
                if m:
                    doFix = True
                    self.updates = True
                    txhead = x.tail[0:m.start()]
                    txtail = x.tail[m.end():]
                    m2 = re.search(refPattern[1],m.groups()[0])
                    span.text = "\\(%s\\)" % (m2.groups()[0])
                    # We have to add to the node and shift text around
                    node.insert(len(node.getchildren()),span)
                    x.tail = txhead
                    span.tail = txtail
                    #import pdb; pdb.set_trace()
                    #aatail = 0

        #print("E>",html.tostring(node))
        #import pdb; pdb.set_trace()

    def insertRefNode(self, refPattern, node, aNode):

        nodect = 0
        nodetotal = len([x for x in node.iter()])

        # We want to iterate over all text within the node no
        # matter how deep
        for x in node.iter():

            # For DOM, we have to convert
            # from: <p>Before \eqref{aref} and after</p>
            #   to: <p>Before <a href="">(X)</a> and after</p>
            #
            # For the TEXT part of a node
            # DOM from:
            #  <p>.text = "Before \eqref{aref} and after"
            #  <p>.tail = None
            #  len(<p>.getchildren()) = 0
            #
            # DOM to:
            #  <p>.text = "Before"
            #  <p>.tail = None
            #  len(<p>.getchildren()) = 1
            #  <p>.getchildren()[0] = <a>
            #  <a>.text = "(X)"
            #  <a>.tail = " and after"

            # For the TAIL part of a node
            # the manipulation is similar
            # 
            # from: <p>Before <a href="">(X)</a> and after \eqref{bref} iterations</p>
            #   to: <p>Before <a href="">(X)</a> and after <a href="">(Y)</a> iterations</p>

            # DOM from:

            # DOM to:

            # Check the text portion of the node
            # In this scenario, we always insert to the front of the
            # node: x.insert(0, aNode), and update the tail of the new
            # child which resides at slot [0].
            m = re.search(refPattern,x.text)
            if m:
                if self.verbose:
                    print('%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                txhead = x.text[0:m.start()]
                txtail = x.text[m.end():]
                if aNode.tail and len(aNode.tail) > 0:
                    txtail = "%s%s" % (aNode.tail, txtail)

                #print('text>',html.tostring(node))
                #import pdb; pdb.set_trace()

                x.text = txhead
                x.insert(0, aNode)
                cNode = x.getchildren()[0]

                #print(html.tostring(node))
                #import pdb; pdb.set_trace()

                # We have to adjust the tail of the child
                # we just inserted
                if len(txtail) > 0:
                    cNode.tail = txtail

                #print(html.tostring(node))
                #import pdb; pdb.set_trace()

                return

            # For the tail portion of a node, this assumes we are already
            # a child of the node we are updating.
            # We have to get the position of the child in the node.
            # Insert a child after the child we found.
            # Adjust the tail of the inserted child.
            m = re.search(refPattern,x.tail)
            if m:
                if self.verbose:
                    print('%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                # For this match, we split the text between the two tails
                txhead = x.tail[0:m.start()]
                txtail = x.tail[m.end():]
                if aNode.tail and len(aNode.tail) > 0:
                    txtail = "%s%s" % (aNode.tail, txtail)
                #print('tail>',html.tostring(x.getparent()))
                #import pdb; pdb.set_trace()
                x.tail = txhead

                # We have to determine the position of this child in the parent
                childpos = -1
                if x in x.getparent().getchildren():
                    childpos = x.getparent().getchildren().index(x)
                x.getparent().insert(childpos+1, aNode)

                #print(html.tostring(x.getparent()))
                #import pdb; pdb.set_trace()

                # Update the new childNode
                childNode = x.getparent().getchildren()[childpos+1]
                if len(txtail) > 0 :
                    childNode.tail = txtail

                #print('after>',html.tostring(x.getparent()))
                #import pdb; pdb.set_trace()

    def updateTarget(self, node, m, fn):
        # We do not need to do crazy things here as this is typically
        # not nested within the same node.
        if self.verbose:
            print('  updateTarget>%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
        repl = m.group(0)
        # The real item to replace is "\\label{}\\notag \\\\ "
        replStr = "%s \\\\ " % (repl)
        node.text = node.text.replace(replStr, '')
        self.updates = True

        # Sphinx updates
        if self.buildType == 'sphinx':
            node.tail = node.tail.replace(replStr, '')
            tag_string = m.groups()[0]
            fc = tag_string.find(":")
            if fc >= 0:
                tag = "equation-%s" % (tag_string[fc+1:])
                if not(tag in self.meta['targets']):
                    self.meta['targets'].append(tag)
                    self.meta['target'][tag] = fn
                try:
                    eqno = self.meta['targets'].index(tag)+1
                except:
                    eqno = 0
            # Only update references that have a formula target
            if eqno > 0:
                eqstr = "(%d)" % (eqno)
                node.text = eqstr
            else:
                if self.verbose: print("WARNING: no target for %s" % (tag))
            #import pdb; pdb.set_trace()

        # Doxygen update targets
        if self.buildType == 'doxygen':
            tag_string = m.groups()[0]
            fc = tag_string.find(":")
            if fc >= 0:
                tag = "equation-%s" % (tag_string[fc+1:])
                # Doxygen 1.8.13 & to match sphinx syntax for tags
                tag = tag.lower()
                tag = tag.translate(tag.maketrans(':_','--'))
                if not(tag in self.meta['targets']):
                    self.meta['targets'].append(tag)
                    self.meta['target'][tag] = fn
                # This is now fixed and should be added to etarget now
                if not(tag in self.meta['etarget'].keys()):
                    self.meta['etarget'][tag] = fn
                try:
                    eqno = self.meta['targets'].index(tag)+1
                except:
                    eqno = 0

                # Don't set the equation node, set the parent id as in below
                #node.set('id',tag)
                # this also needs some node manipulation to add the
                # number and a permalink tag
                prevNode = node.getprevious()
                ele = html.Element("center")
                ele.set("id",tag)
                ele.set("class",'math notranslate nohighlight')
                prevNode.addnext(ele)
                eleNode = prevNode.getnext()

                span = html.Element("span")
                span.set("class","eqno")
                span.text = "(%d)" % (eqno)
                eleNode.insert(0, span)

                spanNode = eleNode.getchildren()[0]
                a = html.Element("a")
                a.set('class','headerlink')
                a.set('href','#%s' % (tag))
                a.set('title','Permalink to this equation')
                a.text = '¶'
                spanNode.insert(0, a)

    def checkSphinxLinks(self, target):
        # sphinx
        # Target
        # <div class="math notranslate nohighlight" id="equation-ale-equations">
        # <span class="eqno">(1)<a class="headerlink" href="#equation-ale-equations"
        # title="Permalink to this equation">¶</a></span>
        # Ref/Link
        # <a class="reference internal" href="api/generated/pages/Specifics.html#hydrostatic-balance"
        #>Hydrostatic balance</a>
        # NOTE: Only check/change text() target if it starts with '('.
        htmlFiles = self.meta['eref'][target]
        for htmlFile in htmlFiles:
            tree = html.parse(htmlFile)
            self.updates = False
            node = tree.xpath("//a[@class='reference internal' and contains(@href,'%s')]" % (target))
            if node:
                node = node[0]
                txt = node.text
                if txt == None:
                    continue
                if txt.find('(') == 0:
                    try:
                        eqno = self.meta['targets'].index(target)+1
                    except:
                        eqno = 0
                        #import pdb; pdb.set_trace()
                    if eqno > 0:
                        eqnoStr = "(%s)" % (eqno)
                        if eqnoStr != node.text:
                            #import pdb; pdb.set_trace()
                            node.text = eqnoStr
                            self.updates = True
                if self.updates:
                    # Write tree back out to file
                    output = html.tostring(tree)
                    #import pdb; pdb.set_trace()
                    fo = open(htmlFile, 'wb')
                    fo.write(output)
                    fo.close()


    def checkSphinxTargets(self, target):
        # sphinx
        # <div class="math notranslate nohighlight" id="equation-ale-equations">
        # <span class="eqno">(1)<a class="headerlink" href="#equation-ale-equations"
        # title="Permalink to this equation">¶</a></span>

        # TODO: This is always a list of one, rework this routine
        htmlFiles = [self.meta['etarget'][target]]
        for htmlFile in htmlFiles:
            tree = html.parse(htmlFile)
            self.updates = False
            node = tree.xpath("//div[@id='%s']/span" % (target))
            if node:
                node = node[0]
                try:
                    eqno = self.meta['targets'].index(target)+1
                except:
                    eqno = 0
                    #import pdb; pdb.set_trace()
                if eqno > 0:
                    eqnoStr = "(%s)" % (eqno)
                    if eqnoStr != node.text:
                        node.text = eqnoStr
                        self.updates = True
            if self.updates:
                # Write tree back out to file
                output = html.tostring(tree)
                #import pdb; pdb.set_trace()
                fo = open(htmlFile, 'wb')
                fo.write(output)
                fo.close()

    def fixEquationTargets(self):
        # fixdiv
        if self.verbose:
            print("fixdiv")
        targets = list(self.meta['fixdiv'].keys())
        for target in targets:
            htmlFiles = self.meta['fixdiv'][target]
            ntag = target
            if ntag.find('equation-') != 0:
                ntag = "equation-%s" % (ntag)
            for htmlFile in htmlFiles:
                tree = html.parse(htmlFile)
                self.updates = False
                if self.verbose:
                    print(" - %s fixdiv tag %s" % (os.path.basename(htmlFile), ntag))
                nodes = tree.xpath("//div[@id='%s']" % (target))
                #import pdb; pdb.set_trace()
                for node in nodes:
                    spanNode = html.Element('span')
                    spanNode.set('class','eqno')
                    try:
                        eqno = self.meta['targets'].index(ntag)+1
                    except:
                        eqno = 0
                    spanNode.text = "(new)"
                    if eqno > 0:
                        eqnoStr = "(%s)" % (eqno)
                        if eqnoStr != spanNode.text:
                            spanNode.text = eqnoStr
                    else:
                        print("Warning: %s tag not found." % (ntag))

                    #import pdb; pdb.set_trace()
                    node.insert(0, spanNode)

                    permNode = node.getchildren()[0]
                    a = html.Element("a")
                    a.set('class','headerlink')
                    a.set('href','#%s' % (ntag))
                    a.set('title','Permalink to this equation')
                    a.text = '¶'
                    permNode.insert(0, a)

                    self.updates = True
                    #import pdb; pdb.set_trace()

                if self.updates:
                    # Write tree back out to file
                    output = html.tostring(tree)
                    #import pdb; pdb.set_trace()
                    fo = open(htmlFile, 'wb')
                    fo.write(output)
                    fo.close()

        # fixtarget: Search for lone labels in formulas (doxygen)
        # Add <center> node prior to formula as a target
        # <center id="equation-ale-equations" class="math notranslate nohighlight">
        # <span class="eqno">(1)<a class="headerlink" href="#equation-ale-equations"
        # title="Permalink to this equation">&#182;</a></span></center>
        if self.verbose:
            print("fixtarget")
        targets = list(self.meta['fixtarget'].keys())
        for target in targets:
            htmlFiles = self.meta['fixtarget'][target]
            ntag = target
            tag = ''
            for htmlFile in htmlFiles:
                if self.verbose:
                    print("  fixtarget> %s %s" % (ntag, os.path.basename(htmlFile)))
                tree = html.parse(htmlFile)
                self.updates = False
                nodes = tree.xpath("//p[@class='formulaDsp']")
                # We have to hunt for the node with the tag we need to update
                for node in nodes:
                    fullText = "%s%s" % (node.text, node.tail)
                    for m in re.finditer('\\\label{(.*?)}', fullText):
                        tag = m.groups()[0]
                        if tag == ntag:
                            break
                    if tag == ntag:
                        break

                # We found the tag in the appropriate formula
                # Add <center> node and update the file
                # Just before the <p> node
                if tag == ntag:
                    # Maybe newer doxygen?
                    #if tag.find('eq:') == 0:
                    #    tag = tag.replace('eq:','equation-')
                    # Doxygen 1.8.13 & to match sphinx syntax for tags
                    tag = tag.lower()
                    tag = tag.translate(tag.maketrans(':_','--'))
                    if tag.find('equation-') != 0:
                        tag = "equation-%s" % (tag)
                    if not(tag in self.meta['targets']):
                        self.meta['targets'].append(tag)
                        if not(tag in self.meta['target'].keys()):
                            self.meta['target'][tag] = htmlFile
                    # This is now fixed and should be added to etarget now
                    if not(tag in self.meta['etarget'].keys()):
                        self.meta['etarget'][tag] = htmlFile

                    #prevNode = node.getprevious()
                    ele = html.Element("center")
                    ele.set("id",tag)
                    ele.set("class",'math notranslate nohighlight')
                    #prevNode.addnext(ele)
                    node.addprevious(ele)
                    eleNode = node.getprevious()

                    #import pdb; pdb.set_trace()
                    try:
                        eqno = self.meta['targets'].index(tag)+1
                    except:
                        eqno = 0
                    span = html.Element("span")
                    span.set("class","eqno")
                    span.text = "(%d)" % (eqno)
                    eleNode.insert(0, span)

                    #import pdb; pdb.set_trace()
                    spanNode = eleNode.getchildren()[0]
                    a = html.Element("a")
                    a.set('class','headerlink')
                    a.set('href','#%s' % (tag))
                    a.set('title','Permalink to this equation')
                    a.text = '¶'
                    spanNode.insert(0, a)
                    self.updates = True

                # Update the file if there were updates
                if self.updates:
                    # Write tree back out to file
                    output = html.tostring(tree)
                    #import pdb; pdb.set_trace()
                    fo = open(htmlFile, 'wb')
                    fo.write(output)
                    fo.close()

        # fixanchors
        # Replace lone anchors with a <center> node
        # Prev:
        # <a class="anchor" id="ht-equation"></a>
        # Replaced:
        # <center id="equation-ht-equation" class="math notranslate nohighlight">
        # <span class="eqno">(1)<a class="headerlink" href="#equation-ht-equation"
        # title="Permalink to this equation">&#182;</a></span></center>
        # NOTE: We fix the tag name here (if needed)
        if self.verbose:
            print("fixanchor")
        targets = list(self.meta['fixanchor'].keys())
        for target in targets:
            htmlFiles = self.meta['fixanchor'][target]
            tag = target
            for htmlFile in htmlFiles:
                if self.verbose:
                    print("  fixanchor> %s %s" % (tag, os.path.basename(htmlFile)))
                tree = html.parse(htmlFile)
                self.updates = False
                nodes = tree.xpath("//a[@class='anchor' and @id='%s']" % (tag))

                # Doxygen 1.8.13 & to match sphinx syntax for tags
                tag = tag.lower()
                tag = tag.translate(tag.maketrans(':_','--'))
                if tag.find('equation-') != 0:
                    tag = "equation-%s" % (tag)
                if not(tag in self.meta['targets']):
                    self.meta['targets'].append(tag)
                    if not(tag in self.meta['target'].keys()):
                        self.meta['target'][tag] = htmlFile
                # This is now fixed and should be added to etarget now
                if not(tag in self.meta['etarget'].keys()):
                    self.meta['etarget'][tag] = htmlFile

                for node in nodes:
                    # Convert found node to center and then setup
                    # the rest
                    node.tag = 'center'
                    node.set("id",tag)
                    node.set("class",'math notranslate nohighlight')
                    eleNode = node

                    #import pdb; pdb.set_trace()
                    try:
                        eqno = self.meta['targets'].index(tag)+1
                    except:
                        eqno = 0
                    span = html.Element("span")
                    span.set("class","eqno")
                    span.text = "(%d)" % (eqno)
                    eleNode.insert(0, span)

                    #import pdb; pdb.set_trace()
                    spanNode = eleNode.getchildren()[0]
                    a = html.Element("a")
                    a.set('class','headerlink')
                    a.set('href','#%s' % (tag))
                    a.set('title','Permalink to this equation')
                    a.text = '¶'
                    spanNode.insert(0, a)
                    self.updates = True

                # Update the file if there were updates
                if self.updates:
                    # Write tree back out to file
                    output = html.tostring(tree)
                    #import pdb; pdb.set_trace()
                    fo = open(htmlFile, 'wb')
                    fo.write(output)
                    fo.close()

                #import pdb; pdb.set_trace()

    def updateEquationLinks(self):
        # Update targets (doxygen)
        # Update links/references (sphinx)
        # Create a list so we can update the dictionary on the fly
        targets = list(self.meta['target'].keys())
        for target in targets:
            fn = self.meta['target'][target]
            self.updates = False
            tree = html.parse(fn)
            nodes = tree.xpath("//*[text()]")
            if len(nodes) > 0:
                if self.verbose:
                    print(fn,len(nodes))
                for node in nodes:
                    # Keep looping until there are no more updates to find
                    scanText = True
                    while scanText:
                        fullText = "%s%s" % (node.text, node.tail)
                        # Set to False unless an update is detected
                        scanText = False
                        m = re.search('\\\\label{(html:.*?)}\\\\notag', fullText)
                        if m:
                            self.updateTarget(node, m, fn)
                            scanText = True
                        # Check sphinx links/references
                        if node.tag == 'a' and node.get('class') == 'reference internal':
                            tag = node.get('href')
                            if tag.find('#equation') == 0:
                                tag = tag[1:]
                                try:
                                    eqno = self.meta['targets'].index(tag)+1
                                except:
                                    eqno = 0

                                # Only update references that have a formula target
                                if eqno > 0:
                                    eqnoStr = "(%s)" % (eqno)
                                    if eqnoStr != node.text:
                                        node.text = eqnoStr
                                        scanText = True
                                else:
                                    if self.verbose: print("WARNING: no target for %s" % (tag))

            # After all the updates are done, if the flag to update is
            # set, rewrite the html file.
            if self.updates:
                # Write tree back out to file
                output = html.tostring(tree)
                #import pdb; pdb.set_trace()
                fo = open(fn, 'wb')
                fo.write(output)
                fo.close()

        # Check Sphinx targets
        for target in self.meta['etarget'].keys():
            self.checkSphinxTargets(target)

        # Check Sphinx references
        for ref in self.meta['eref'].keys():
            self.checkSphinxLinks(ref)

        # Update links (doxygen & sphinx?)
        refPattern = '\\\eqref{(.*?)}'
        refPattern2 = '\\\eqref2{(.*?)}'
        for ref in self.meta['ref'].keys():
            fn = ref
            self.updates = False
            tree = html.parse(fn)
            # Manipulation of the document is very hard.  We need to rescan
            # the entire document each time there is a change in structure
            # to find all the updates.   To seed this loop, we pass a -1.
            updatesFound = -1
            npass = 0
            while updatesFound != 0:
                updatesFound = 0
                npass = npass + 1
                if self.buildType in ('doxygen','sphinx'):
                    nodes = tree.xpath("//*[text()]")
                    if len(nodes) > 0:
                        if self.verbose:
                            print(fn,npass,len(nodes))
                        for node in nodes:
                            if node.text == None:
                                continue
                            if node.tail:
                                fullText = "%s%s" % (node.text, node.tail)
                            else:
                                fullText = node.text

                            m = re.search(refPattern, fullText)
                            # \eqref
                            if m:
                                if self.verbose:
                                    print('  eqref>%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                                repl = m.group(0)
                                tag = m.groups()[0]
                                # Newer doxygen?
                                #if tag.find('eq:') != -1:
                                #    tag = tag.replace('eq:','')
                                # Doxygen 1.8.13
                                tag = tag.lower()
                                tag = tag.translate(tag.maketrans(':_','--'))
                                fullTag = tag
                                if tag.find('equation-') != 0:
                                    fullTag = "equation-%s" % (tag)

                                try:
                                    computePath = os.path.relpath(os.path.dirname(self.meta['etarget'][fullTag]),os.path.dirname(fn))
                                except:
                                    print("WARNING: Target tag not found: %s" % (fullTag))
                                    #import pdb; pdb.set_trace()
                                    continue

                                # If we are in the same directory, do not specify a path
                                if computePath == ".":
                                    computePath = ""
                                else:
                                    computePath = "%s/" % (computePath)
                                computePath = "%s#%s" % (os.path.basename(self.meta['etarget'][fullTag]),fullTag)
                                # Replace \\eqref{} with <a href=""></a> that looks like sphinx
                                try:
                                    eqno = self.meta['targets'].index(fullTag)+1
                                except:
                                    self.meta['targets'].append(tag)
                                    eqno = self.meta['targets'].index(fullTag)+1

                                aNode = html.Element('a')
                                aNode.set('href',computePath)
                                aNode.set('class','reference internal')
                                aNode.text = "(%d)" % (eqno)
                                aNode.tail = ''
                                self.insertRefNode(refPattern, node, aNode)
                                updatesFound = updatesFound + 1
                                self.updates = True

                                # This is now a cleaned up link to a target
                                if (not fullTag in self.meta['eref'].keys()):
                                    self.meta['eref'][fullTag] = []
                                if (not fn in self.meta['eref'][fullTag]):
                                    self.meta['eref'][fullTag].append(fn)


                            # \eqref2
                            m = re.search(refPattern2, fullText)
                            if m:
                                if self.verbose:
                                    print('  eqref2>%03d-%03d: %s' % (m.start(), m.end(), m.group(0)))
                                repl = m.group(0)
                                tag_string = m.groups()[0]
                                fc = tag_string.find(',')
                                if fc >= 0:
                                    tag = tag_string[0:fc]
                                    fullTag = "equation-%s" % (tag)
                                    try:
                                        computePath = os.path.relpath(os.path.dirname(self.meta['target'][tag]),os.path.dirname(fn))
                                    except:
                                        print("WARNING: Target tag not found: %s" % (tag))
                                        #import pdb; pdb.set_trace()
                                        continue

                                    # If we are in the same directory, do not specify a path
                                    if computePath == ".":
                                        computePath = ""
                                    else:
                                        computePath = "%s/" % (computePath)
                                    computePath = "%s#%s" % (os.path.basename(self.meta['target'][tag]),fullTag)

                                    # Replace \\eqref2{} references that translate to something similar to:
                                    # <a href="General_Coordinate.html#equation-h-equations" 
                                    #   class="reference internal">(7)</a> - momentum</p>
                                    try:
                                        eqno = self.meta['targets'].index(fullTag)+1
                                    except:
                                        eqno = 0

                                    aNode = html.Element('a')
                                    aNode.set('href',computePath)
                                    aNode.set('class','reference internal')
                                    aNode.text = "(%d)" % (eqno)
                                    aNode.tail = " - %s" % (tag_string[fc+1:])
                                    self.insertRefNode(refPattern2, node, aNode)
                                    updatesFound = updatesFound + 1
                                    self.updates = True

                                    # This is now a cleaned up link to a target
                                    if (not fullTag in self.meta['eref'].keys()):
                                        self.meta['eref'][fullTag] = []
                                    if (not fn in self.meta['eref'][fullTag]):
                                        self.meta['eref'][fullTag].append(fn)

            if self.updates:
                # Write tree back out to file
                #import pdb; pdb.set_trace()
                output = html.tostring(tree)
                fo = open(fn, 'wb')
                fo.write(output)
                fo.close()

        # Collect a list of files to scan and update existing targets and references
        #targets = self.meta['checkFiles']
        #for target in targets:
        #    import pdb; pdb.set_trace()
        #    continue

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", action="store_true", help="Turn on verbosity")
parser.add_argument("-d", "--dir", type=str, help="Root directory")
parser.add_argument("-p", "--project", type=str, help="Project directory")
parser.add_argument("-b", "--build", type=str, help="Build type")
parser.add_argument("-s", "--start", type=str, help="Start file")
parser.add_argument("-l", "--links", action="store_true", help="Show targets, links and references")

args = parser.parse_args()

rootDirectory = '/var/www/html'
verbose = False
showLinks = False
startFile = 'index.html'
projectDirectory = None
buildType = 'doxygen'

if args.verbose:
    verbose = True

if args.links:
    showLinks = True

if args.dir:
    rootDirectory = args.dir

if args.project:
    projectDirectory = args.project

if args.start:
    startFile = args.start

if args.build:
    buildType = args.build

if os.path.isdir(rootDirectory):
    #if verbose: print("Found root:",rootDirectory)
    os.chdir(rootDirectory)
    if os.path.isdir(projectDirectory):
        #if verbose: print("Found project:",projectDirectory)
        mathProc = equationRenumber(projectDirectory, buildType)
        mathProc.verbose = verbose
        if mathProc.verbose:
            print("** getHtmlFiles")
        mathProc.getHtmlFiles()
        if mathProc.verbose:
            print("** htmlWalk")
        mathProc.htmlWalk(startFile)
        if mathProc.verbose:
            print("** collectEquationLabels")
        mathProc.collectEquationLabels()
        if mathProc.verbose:
            print("** fixEquationLinks")
        mathProc.fixEquationTargets()
        if mathProc.verbose:
            print("** updateEquationLinks")
        mathProc.updateEquationLinks()
    else:
        print("ERROR: Project directory not found (%s). Exiting." % (projectDirectory))
        sys.exit(1)
else:
    print("ERROR: Root directory not found (%s). Exiting." % (rootDirectory))
    sys.exit(1)

if showLinks:
    #import pdb; pdb.set_trace()
    mathProc.meta['targets'].sort()
    print()
    msg = "Pages with formula targets and references"
    print(msg)
    print('-' * len(msg))
    for target in mathProc.meta['targets']:
        # Not all targets are formulas
        # Remove equation- prefix
        tag = target
        # There should be only one target per tag
        pages = []
        if tag in mathProc.meta['etarget'].keys():
            pages.append(mathProc.meta['etarget'][tag])

        for page in pages:
            print("target>%s: %s" % (tag, page))

            linkedPages = []
            if tag in mathProc.meta['ref'].keys():
                for lpage in mathProc.meta['ref'][tag]:
                    if not(lpage in linkedPages):
                        linkedPages.append(page)
            if tag in mathProc.meta['eref'].keys():
                for lpage in mathProc.meta['eref'][tag]:
                    if not(lpage in linkedPages):
                        linkedPages.append(page)
            for lPage in linkedPages:
                print("  %s" % (lPage))

        #import pdb; pdb.set_trace()
