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
# DATE  : 2011-06
# REV   : 1
# REQUIREMENTS: 
#
#############################################################################
from pymol import cmd


def construct_dico1(liste_entree=['vide']):
    """
        permet de construire un dictionnaire ‡ partir d'une liste d'ÈlÈments donnÈs
    """
    nvdico={}
    for elmt in liste_entree:
        print elmt
        nvdico[elmt[0]]=elmt[1]
    print "construction du dictionnaire terminee" 
    return nvdico


def construct_dico2(liste_entree=['vide']):
    """
        permet de construire un dictionnaire ‡ partir d'une liste d'ÈlÈments donnÈs
    """
    nvdico={}
    liste_mod=[]
    for elmt in liste_entree:
        liste_mod.append("d"+elmt+"_s1")
    for elmt in liste_mod:
        nvdico[elmt[0]]=elmt[1]
    print "construction du dictionnaire terminee" 
    return nvdico

def construct_dico(liste_entree=['vide'],dir=".",opt="defaut"):
    """
        permet de construire un dictionnaire ‡ partir d'une liste d'ÈlÈments donnÈs
    """
    
    if liste_entree==['vide']:
        print "construct_dico([liste_structure],repertoireDesSolutions(facultatif),option(facultatif))"
        print "repertoireDesSolutions = repertoire courant si non indiquer"
        print "option si non indiquee utilisation de la premiere solution utilisee (generelement solution complete ou meilleure solution le cas echeant)"
        print "       si option expert utilisee : recherche d'une solution(script) avec le nom exact saisi"
    
    #initialisation des variables
    nvdico={}
    liste_mod=[]
    prefix="d"
    suffix=""
    if opt=="defaut":
        suffix="_sol_1"
    
    #creation des nouveaux noms
    for elmt in liste_entree:
        liste_mod.append(prefix+elmt+suffix)
    
    #ouverture des fichiers necessaires et memorisation des informations
    for elmt in liste_mod:
        cmd.do("run "+dir+"/"+elmt[1:]+".py")
        nvdico[elmt[0]]=elmt[1]
    print "construction du dictionnaire terminee" 
    return nvdico


def construct_dico_int(liste_entree=['vide']):
    """
        permet de construire un dictionnaire ‡ partir d'une liste d'ÈlÈments donnÈs
    """
    nvdico={}
    liste_mod=[]
    for elmt in liste_entree:
        liste_mod.append("d"+elmt+"_s1")
    for elmt in liste_mod:
        nvdico[elmt[0]]=elmt[1]
    print "construction du dictionnaire terminee" 
    return nvdico


cmd.extend('construct_dico',construct_dico)

def removeAlt(obj="(all)", keep="A"):
        """
        from PyMOLWIKI
        removeAlt -- remove all alternate location-atoms not of altloc "keep" from object.
 
        input:
                obj -- the object(s) to remove the atoms frmo
                keep -- which type of alt loc to keep
 
        output: none -- removes atoms
 
        examples:
                removeAlt # remove all altLocations that aren't altloc A
                removeAlt pdbID, C  # remove all but C altlocations from pdbID
        """
        
        remStr = "%s and not (alt ''+%s)" % (obj, keep);
        cmd.remove(remStr);

cmd.extend('removeAlt',removeAlt)

#cmd.align("1ACJ_stride&(c;A&i;200&n;ca)","1TCA_stride&(c;A&i;105&n;ca)")
