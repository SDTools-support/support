from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst import roles

from sphinx.application import Sphinx

class helloworlddirective(Directive):
    def run(self):
        paragraph_node = nodes.paragraph(text='Hello World Directive!')
        return [paragraph_node]

def helloworldrole(role, rawtext, text, lineno, inliner,
                       options=None, content=None):
    text='Hello World Role : <font color=”green”>' + text + '</font>'
    node=nodes.Text(text)
    return [node], []
       
def setup(app):
    app.add_directive('helloworlddirective', helloworlddirective)
    app.add_role('helloworldrole', helloworldrole)