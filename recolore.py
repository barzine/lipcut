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
 
from pymol import cmd
import tkSimpleDialog
import tkMessageBox
import sys
import pymol 


#d2LIP_stride_sol_1=('2LIP_stride',
#{
#'Ser' : '2LIP_stride&(c;A&i;87&n;ca)',
#'B5' : (('3:3:0', 'VAL', 'A.0', '81', 'HIS', 'A.0', '86', '2', 1),['2LIP_stride&(c;A&i;83&n;ca)','2LIP_stride&(c;A&i;84&n;ca)','2LIP_stride&(c;A&i;85&n;ca)','2LIP_stride&(c;A&i;86&n;ca)']),
#'B6' : (('4:4:0', 'VAL', 'A.0', '104', 'ILE', 'A.0', '110', '2', 1),['2LIP_stride&(c;A&i;107&n;ca)','2LIP_stride&(c;A&i;108&n;ca)','2LIP_stride&(c;A&i;109&n;ca)','2LIP_stride&(c;A&i;110&n;ca)']),
#'B7' : (('6:6:0', 'ASN', 'A.0', '202', 'GLY', 'A.0', '211', '2', 1),['2LIP_stride&(c;A&i;206&n;ca)','2LIP_stride&(c;A&i;207&n;ca)','2LIP_stride&(c;A&i;208&n;ca)','2LIP_stride&(c;A&i;209&n;ca)']),
#'B8' : (('9:9:0', 'GLN', 'A.0', '276', 'TYR', 'A.0', '282', '2', 1),['2LIP_stride&(c;A&i;276&n;ca)','2LIP_stride&(c;A&i;277&n;ca)','2LIP_stride&(c;A&i;278&n;ca)','2LIP_stride&(c;A&i;279&n;ca)']),
#'B4' : (('1:1:0', 'ILE', 'A.0', '11', 'VAL', 'A.0', '14', '2', 1),['2LIP_stride&(c;A&i;11&n;ca)','2LIP_stride&(c;A&i;12&n;ca)','2LIP_stride&(c;A&i;13&n;ca)','2LIP_stride&(c;A&i;14&n;ca)']),
#'B3' : (('2:2:0', 'VAL', 'A.0', '44', 'ALA', 'A.0', '47', '2', 1),['2LIP_stride&(c;A&i;44&n;ca)','2LIP_stride&(c;A&i;45&n;ca)','2LIP_stride&(c;A&i;46&n;ca)','2LIP_stride&(c;A&i;47&n;ca)']),
#'His' : ['2LIP_stride&(c;A&i;286&n;ca)'],
#'Acide' : ['2LIP_stride&(c;A&i;264&n;ca)','2LIP_stride&(c;A&i;288&n;ca)','2LIP_stride&(c;A&i;289&n;ca)']
#})




def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                            'Recolore',
                            label = 'Recolore',
                            command = lambda s=self : RecoloreDialog(s))




def colore(liste, couleur):
    for elmt in liste:
        cmd.color(couleur,elmt)


def Recolore(nom):
    
    try :
        exec 'struct =pymol.d%s_sol_1[1]' %nom #permet de mÈmoriser le contenu de la variable "nom" de pymol dans struct
    except:
        print "Erreur inatendue:", sys.exc_info()[0]
        tkMessageBox.showerror('ID invalide',
                             'l objet specifie n existe pas '+nom)
        
    #coloration de la structure si recuperation rÈussie
    print "Recoloration ...."
    cmd.color("chocolate", struct['Ser'])
    cmd.show("spheres",struct['Ser'])
    for his in struct['His']:
        cmd.show("spheres",his)
        cmd.color("lightorange",his)
    for aa in struct['Acide']:
        cmd.show("spheres",aa)
        cmd.color("wheat",aa)
    
    dic={'B2' : "lightpink",'B3' : "magenta",'B4' : "hotpink",'B5' : "red",'B6' : "orange",'B7' : "tv_orange",'B8' : "brightorange"}
    cle=struct.keys()
    
    for c in cle:
        if c in ['B2','B3','B4','B5','B6','B7','B8']:
            #le tuple considÈrÈ est celui d'un brin, on rÈcupere la liste pour colorer
            colore(struct[c][1],dic[c])
    print "finie"
        
        

def RecoloreDialog(app):
    nom = tkSimpleDialog.askstring('ID structure',
                                      'Entrez le nom de l objet SVP:',
                                      parent=app.root)
 
    Recolore(nom)


cmd.extend('colore',colore)
cmd.extend('Recolore', Recolore)
