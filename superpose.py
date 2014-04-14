#!python
# -*-coding:Latin-1 -*
##############################################################################
#
# @SUMMARY: -- 
#             
# Compatibility : teste seulement sur PyMOL 0.99 sous windows XP - python v2.4
# @AUTHOR: M. P. Barzine
# @COPYRIGHT: M. P. Barzine (C), V. TRAN (C), 2011
# @LICENSE: Released under GPL:
# This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA 
#
# DATE  : 2011-03
# REV   : 1
# REQUIREMENTS: 
#
#############################################################################
import tkSimpleDialog
import tkMessageBox
from pymol import cmd
import pymol
import sys

def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                            'Superpose brins',
                            label = 'Superpose brins',
                            command = lambda s=self : fetchSuperposeDialog(s))
 
def superpose(struct1,struct2):

    """
        permet de superposer par pair_fitting deux structures par rapport aux diffÈrents rÈsidus sÈlectionnÈs
    """
    select=""
    #select1=""
    #select2=""
    print "on fait des tests"
   
    exec 'struct1 =pymol.d%s_sol_1[1]' %struct1 #permet de mÈmoriser le contenu de la variable "struct1" de pymol dans struct1
    exec 'struct2 =pymol.d%s_sol_1[1]' %struct2 #permet de mÈmoriser le contenu de la variable "struct2" de pymol dans struct2
    print "Struct 1 \n\n"
    cle=struct1.keys()
    print "reste plus qu'‡ finir le truc"
    select+=struct1['Ser']+", "+struct2['Ser']
   
    for c in cle: #balayage des diffÈrentes entrÈes du dictionnaire pour la structure 1
        if c in struct2:#si l'ÈlÈment existe dans la structure 2, on va les "pair_fitter"
            #if c in ['B2','B3','B4','B5','B6','B7','B8']:
            if c in ['B4','B5','B6','B7']:
                for i in range(len(struct1[c][1])):
                    select+=", "+struct1[c][1][i]
                    select+=", "+struct2[c][1][i]
    
            
            
    #‡ ce niveau l'ensemble de la premiere structure a ÈtÈ balayÈ et les diffÈrents ÈlÈments similaires ont ÈtÈ mÈmorisÈ dans le mÍme ordre
    #il ne reste plus qu'‡ fitter l'ensemble
                
    if select!="":
        #print "Printage"
        #print "pair_fit "+select
        cmd.do("pair_fit "+select)      


    
    
    
    
    
    
def fetchSuperposeDialog(app):
    s1,s2 = tkSimpleDialog.askstring('PDB Loader Service',
                                      'Please enter a 4-digit pdb code:',
                                      parent=app.root)
 
    superpose(s1,s2)
 
cmd.extend('superpose', superpose)
