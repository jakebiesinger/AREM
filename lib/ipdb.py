"""Description: Ipython Debugger
Author: Godefroid Chapelle <gotcha at bubblenet be>
Home Page: https://github.com/gotcha/ipdb
Keywords: pdb ipython
License: GPL
Package Index Owner: gotcha
DOAP record: ipdb-0.3.xml
"""

import sys
from IPython.Debugger import Pdb
from IPython.Shell import IPShell
from IPython import ipapi

shell = IPShell(argv=[''])


def set_trace():
    ip = ipapi.get()
    def_colors = ip.options.colors
    Pdb(def_colors).set_trace(sys._getframe().f_back)


def post_mortem(tb):
    ip = ipapi.get()
    def_colors = ip.options.colors
    p = Pdb(def_colors)
    p.reset()
    while tb.tb_next is not None:
        tb = tb.tb_next
    p.interaction(tb.tb_frame, tb)


def pm():
    post_mortem(sys.last_traceback)
