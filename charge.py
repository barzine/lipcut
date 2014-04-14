# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2004 by Charles Moad <cmoad@indiana.edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
 
 
import tkSimpleDialog
import tkMessageBox
from pymol import cmd
import sys
 
def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                            'Charge',
                            label = 'Charge Lipase/Cutinase',
                            command = lambda s=self : chargePDBDialog(s))
 
def charge(pdbCode):
    pdbCode = pdbCode.upper()
    try:
        cmd.do("run "+pdbCode+"_stride_sol_1.py")
        cmd.color("deepteal",pdbCode+"_stride")
        cmd.do("recolore "+pdbCode+"_stride")
        #if (pymol.d1LPBS_stride_sol_1) :
        #    cmd.do("superpose "+pdbCode+"_stride, 1LPBS_stride")
        cmd.do("center")
        
    except:
        print "Unexpected error:", sys.exc_info()[0]
        tkMessageBox.showerror(' Code invalide',
                           'Cette entree n existe pas dans la base resultat :' + pdbCode)
 
def chargePDBDialog(app):
    pdbCode = tkSimpleDialog.askstring('Charge Lipase/Cutinase',
                                      'ID_PDB code:',
                                      parent=app.root)
 
    charge(pdbCode)
 
cmd.extend('charge', charge)
