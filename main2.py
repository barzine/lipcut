#!python
# -*-coding:Latin-1 -*
##############################################################################
#
# @SUMMARY: -- 
#             
# Compatibility : teste seulement sur PyMOL 0.99 sous windows XP - python v2.4
# @AUTHOR: M. P. Barzine
# @COPYRIGHT: M. P. Barzine (C), 2011
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
# DATE  : 2011-05
# REV   : 1
# REQUIREMENTS: 
#
#############################################################################

#importation des packages propres ‡ python, extÈrieurs 
#et dÈveloppÈs spÈcialement pour cette application
####packages propres ‡ python:
import os
#from glob import glob
from pymol import cmd
import string
#####packages dÈveloppÈs spÈcialement pour l'applicatio
#import utils

####

###definition des procedures et fonctions:

def alignement_sur_brins(ref,struct):
    """
        permet d'aligner les deux structures en fonctions des divers ÈlÈments utilisÈs pour 
    """

    




def compte_nb_helice_similaire(ref,struct):
    #alignement des structures selon les brins beta
    alignement_sur_brins(ref,struct)
    #print """cmd.super(string.split(os.path.basename(r),".")[0]+'&ss s',string.split(os.path.basename(f),".")[0]+'&ss s')"""
    #cmd.do('super '+ref+'&c. A&ss s, '+struct+'&c. A&ss s' )
    #x=cmd.super(ref+'&ss s',struct+'&ss s')
    x="resultat RMS apres superposition avec pair_fitting "
    print x
    

#Corps du programme
#################################################
#dÈfinition du rÈpertoire contenant les donnÈes sources :
#Structures = glob("../mini_base/*.pdb")
#References=glob("../references/*.pdb")
#################################################

#donnÈes d'entrÈes :
Structures="machin" #creation 
References="truc"

#dÈfinition du rÈpertoire des rÈsultats :
Dir_Res="../resultats_comp_helices/"
#creation du dossier resultats s'il n'existe pas
if not(os.path.isdir(Dir_Res)):
    os.makedirs(Dir_Res, mode=0777)#tous les rÈpertoires n'existant pas

#chargement des diffÈrentes structures :
#for f in Files:
#    utils.chargement_fichier_pdb2(f)
#for r in References:
#    utils.chargement_fichier_pdb2(r)
#traitement 
cpt=0 #permet de numÈroter les structures traitÈs
for s in Structures:
    #pour chaque structure :
    #incrÈmentation du numÈro de la structure
    cpt+=1
    print "\n\nnouvelle structure n∞ "+str(cpt)+" : "+s+"\n" #debug : affichage de suivi
    #chargement de la structure:
    #utils.chargement_fichier_pdb2(f)
    cpt_ref=0#mise ‡ zÈro du nombre de reference testÈs
    #parcours de la liste des references:
    for r in References:
    #while True:
        #pour chaque modele de reference: #2LIP, 1TCA, 1CRL, 1CEX, 1LPB, 1CEX
        #incrÈmentation du numÈro de la structure
        cpt_ref+=1
        print "\n\nnouvelle reference n∞ "+str(cpt_ref)+" : "+r+" ("+s+")\n" #debug : affichage de suivi
        #chargement de la structure de reference:
        #utils.chargement_fichier_pdb2(r)
        compte_nb_helice_similaire(string.split(os.path.basename(r),".")[0],string.split(os.path.basename(s),".")[0])
        #cmd.super(string.split(os.path.basename(r),".")[0]+'&ss s',string.split(os.path.basename(f),".")[0]+'&ss s')
        #avant de passer ‡ la reference suivante, la reference actuelle est supprimÈe
        #cmd.delete(string.split(os.path.basename(r),".")[0])
        #break
    #break
    #avant de passer ‡ la structure suivante, la structure actuelle est supprimÈe
    #cmd.delete(string.split(os.path.basename(f),".")[0])    
    print "done for"+s+"\n"

    

print "all done"

