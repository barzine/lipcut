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


from pymol import cmd
from pymol import stored

#import utils
import site_utils
import ser
  
####
#verifier que tous les residus selectionnes sont bien dans les feuillets dÈcrits 
#permet d'avertir si le nombre spÈcifie de residus ‡ sÈlectionner par brin est trop important par rapport ‡ la taille du brin considÈrÈ



def verif_site_catal(serine,d,x=15):
    """
        permet de verifier si la solution consideree est conforme a la topologie utilisee en reference en ce qui concerne le site catalytique
    """

    w=dict()
    #dans un premier temp, on recherche tous les carbones alpha dans un rayon de x angstrom
    #a partir de la serine utilise comme reference
    cmd.select("selection_pour_site","name ca within "+str(x)+" of (resi "+str(serine[0])+" and chain "+serine[1]+" and name ca)")
    stored.list=list()
    cmd.iterate("selection_pour_site","stored.list.append((resi,chain,resn))")
    #print "liste genere par pymol"#debug
    #print stored.list #debug
    
    
    #on recherche dans un deuxieme temps s'il existe une histidine dans cette selection
    his,w["his"]=site_utils.verif_histidine(stored.list,d)
    
    #dans un troisieme temps on recherche un aspartate ou un glutamate idealement place 
    acide,w["acide"]=site_utils.verif_acide(stored.list,d)
    
    w["dist"]=x
    
    cmd.delete("selection_pour_site")
    return [his,acide,w]

###


def estimation_assignation(listeAssign):
    """
        permet d'estimer l'assignation des brins selon la serine choisie comme reference
        permet ainsi d'aider dans le choix a traiter
        
        retourne l'ensemble un peu modifie : un troisieme champ est ajoute aux differents tuples, un dictionnaire de warnings
    """
    
    
    #print "dans la fonction estimation assignation"
    res_complete=list()
    res_incomplete=list()
    
    if len(listeAssign)<1:
        return [['none'],'Il n y a aucune solution a estimer']
    
    for sol in listeAssign:
        ser,dicoB,avert=sol
        #verification du site catalytique
        His_l,Acide_l,avert['site']=verif_site_catal(ser,dicoB)
        #print His_l
        #print Acide_l
        #autres verifications
        if (avert["antiparallelite"] and dicoB["B8"]!=["none"] and His_l!=["none"] and Acide_l!=["none"]) or (dicoB["B_AUTRE"]==['none']and His_l!=["none"] and Acide_l!=["none"]):
            res_complete.append((ser,dicoB,His_l,Acide_l,avert)) 
        else :
            res_incomplete.append((ser,dicoB,His_l,Acide_l,avert))
        print "avert"
        print avert
    #print "les solutions completes sont :"
    #utils.affichage_boucle(res_complete)
    
    #print "\n\nles solutions imcompletes sont :"
    #utils.affichage_boucle(res_incomplete)
    warn=dict()
    if len(res_complete)>0:
        warn['complete']=len(res_complete)
    if len(res_incomplete)>0:
        warn['incomplete']=len(res_incomplete)    
    
    result=res_complete+res_incomplete
    
    
    #print "on sort de la fonction estimation assignation"
    return [result,warn]
