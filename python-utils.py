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
# DATE  : 2014-04
# 
# REQUIREMENTS: 
#
#############################################################################

#libraries and personal files importation
import utils;


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
        #select & remove all non A altlocs
        remStr = "%s and not (alt ''+%s)" % (obj, keep);
        cmd.remove(remStr);
        #reset the PDB information
        #cmd.alter( obj, "alt=''")

def chargement_fichier_pdb(filename):
    """
        permet de charger le fichier pdb et de faire les pretraitements necessaires au bon fonctionnement du script dans pyMOL
    """
    cmd.delete("all") #initialise l'espace de travail dans pymol
    cmd.load(filename) #charge le fichier spÈcifiÈ par filename
    removeAlt()#retire tous jeux de coordonnÈes alternatif pour un atome donnÈ
    cmd.remove("hetatm")#permet de retirer les heteroatoms #attention : parfois certains fichiers sont mal annotÈs : un carbone alpha peu Ítre dÈcrit comme Ètant un hÈtÈro-atome
    cmd.remove("h.")#permet de retirer toutes les molÈcules d'eau
    cmd.hide("all")#permet de cacher l'ensemble de la structure
    #cmd.show("ribbon")#rÈaffiche la structure sous le format ruban
    cmd.show("cartoon")#rÈaffiche la structure sous le format cartoon
    cmd.color("gray30")#colorise l'ensemble de gris - pour faire ressortir ultÈrieurement par contraste les fragments sÈlectionnÈs


def chargement_fichier_pdb2(filename):
    """
        permet de charger le fichier pdb et de faire les pretraitements necessaires au bon fonctionnement du script dans pyMOL
    """
    #cmd.delete("all") #initialise l'espace de travail dans pymol
    cmd.load(filename) #charge le fichier spÈcifiÈ par filename
    removeAlt()#retire tous jeux de coordonnÈes alternatif pour un atome donnÈ
    cmd.remove("hetatm")#permet de retirer les heteroatoms #attention : parfois certains fichiers sont mal annotÈs : un carbone alpha peu Ítre dÈcrit comme Ètant un hÈtÈro-atome
    cmd.remove("h.")#permet de retirer toutes les molÈcules d'eau
    cmd.hide("all")#permet de cacher l'ensemble de la structure
    #cmd.show("ribbon")#rÈaffiche la structure sous le format ruban
    cmd.show("cartoon")#rÈaffiche la structure sous le format cartoon
    #cmd.color("gray30")#colorise l'ensemble de gris - pour faire ressortir ultÈrieurement par contraste les fragments sÈlectionnÈs


def identificateur_pymol(atom):
    """
        a partir d'un tuple permettant de decrire un atome,
        renvoie une chaine de caractere permettant de manipuler cet atome dans pymol
    """
    #champ[1] de atom contient l'identificateur de chaine du residu
    #champ[0] de atom contient le numÈro du rÈsidu, 
    
    return ("(c;"+atom[1]+"&i;"+atom[0]+"&n;ca)")

##PyMOL et visualisation 
def color_struct(elmt,t_elmt):
    """
        permet de colorer la structure specifier selon sa nature, permet ainsi de visualiser facilement les elements selectionnes
        pour un controle visuel de l'utilisateur
    """
    
    #"c;A&i;87&n;ca"
    #cmd.show("spheres","c;A&i;87&n;ca")
    #cmd.color("orange","c;A&i;87&n;ca")
    #cmd.color("red","c;A&i;83&n;ca")
    
    #print "dans la fonction color_struct_beta_ser"#debug
    #print elmt,t_elmt#debug
    
    if t_elmt=="SER":
        cmd.show("spheres",identificateur_pymol(elmt))
        #cmd.color("red",identificateur_pymol(elmt))
    
    elif t_elmt=="SER_CAT":
        cmd.color("chocolate",identificateur_pymol(elmt))
            
    elif t_elmt=="B5":
        for e in elmt:
            cmd.color("red",identificateur_pymol(e))
    
    elif t_elmt=="B6":
        for e in elmt:
            cmd.color("orange",identificateur_pymol(e))
    
    elif t_elmt=="B7":
        for e in elmt:
            cmd.color("tv_orange",identificateur_pymol(e))
    
    elif t_elmt=="B8":
        for e in elmt:
            cmd.color("brightorange",identificateur_pymol(e))
            
    elif t_elmt=="B4":
        for e in elmt:
            cmd.color("hotpink",identificateur_pymol(e))
    
    elif t_elmt=="B3":
        for e in elmt:
            cmd.color("magenta",identificateur_pymol(e))
    
    elif t_elmt=="B2":
        for e in elmt:
            cmd.color("lightpink",identificateur_pymol(e))
    
    #elif t_elmt=="B_AUTRE":
    #    for e in elmt:
    #        cmd.color("blue",identificateur_pymol(e))
    
    elif t_elmt=="B_AUTRE2":
        for e in elmt:
            cmd.color("skyblue",identificateur_pymol(e))
    #print "on sort de la fonction color_struct_beta_ser"#debug
    
    elif t_elmt=="ACIDE":
        cmd.show("spheres",identificateur_pymol(elmt))
        cmd.color("wheat",identificateur_pymol(elmt))
    
    elif t_elmt=="HIS":
        cmd.show("spheres",identificateur_pymol(elmt))
        cmd.color("yellow",identificateur_pymol(elmt))

