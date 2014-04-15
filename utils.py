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
# REV   : 2
# REQUIREMENTS: 
#
#############################################################################

#importation of the python (third-parties and in house) packages needed
from pymol import cmd


#Checking functions and procedures: allow to print the different structures used 
#in the application

#def affichage_boucle(Liste):
def print_loop(Liste):
    """
    Procedure printing the content of one liste, line per line
    """ 
    for elmt in Liste :
        print elmt 
          
##
#def affichage_boucle_liste_dans_list(Liste):
def print_loop_in_loop(Liste):
    """
    Procedure printing line per line the content of lists themselves included in a list 
    """ 
    il=0 #initialisation of the counter for the sub-lists
    for sublist in Liste : #for each sublist in the main list
        il+=1 #increment now since used to reference the sublist
        print "Sublist #: ",il #print the index of the sublist
        ie=0 #initialisation of the counter for the elements of the current sublist
        for elmt in sublist : 
            ie+=1 #increment the counter of the element (so it would display the correct number)
            print elmt #print the current element
        #when end of the sublist is reached, total number of elements in it is displayed
        print "Total number of elements in sublist ",il,": ",ie,"element(s).\n" 
    #when the end of the main list is reached, display of the total number of sublist it contains
    print "\nThe main list contains "+il+" sublist(s)."
##

def affichage_dico(dico):
    """
        permet d'afficher le contenu d'un dico de maniere un peu plus intelligible
    """
    
    cles=dico.keys()
    print "cle, valeur"
    for c in cles:
        print c,dico[c] 



#Fonctions de recopie 
#pour palier au pointage dynamique rÈalisÈ par Python au lieu des affectations souhaitÈes
def copy(liste_source):
    """
        Fonction dÈfinie car python ne copie pas les donnÈes,
        mais crÈe, ‡ la place, un synonyme pour la liste
    """
    #print "dans copy" #debug : affichage de contrÙle
    liste_a_retourner=[] #initialisation de la structure mÈmorisant la liste ‡ retourner
    #print "initialisation effectue" #debug : affichage de contrÙle
    for i in liste_source: #pour chaque element i de la liste_source
    #   print i #debug : affichage de contrÙle
        liste_a_retourner.append(i) #ajout de l'ÈlÈment
    #print "on sort de copy" #debug : affichage de contrÙle
    return(liste_a_retourner) #permet de crÈer une liste identique a celle reÁue en paramËtre


#fonctions travaillant avec PyMOL :

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

