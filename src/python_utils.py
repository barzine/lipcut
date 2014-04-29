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
import utils
from pymol import cmd

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

#def chargement_fichier_pdb(filename):
def load_pdb_file(filename, fresh=True):
    """
        prepare pymol environment and load the pdb file for the pre-treatment before starting the analysis
       
    """
    if (fresh): #new analysis 
        cmd.delete("all") #pyMOL is purged of all the objects
    cmd.load(filename) #load the pdb file
    removeAlt()#every alternate coordinate is removed -- won't make a too big difference anyway
    cmd.remove("hetatm")#remove the heteroatom (normaly ligand
    #### HOWEVER, some alternate carbon alpha have been spotted with the "hetatm" annotation
    #### if any funny thing appears while the analysis, check the pdb file structure
    cmd.remove("h.")#remove all the water molecules
    cmd.hide("all")#hide the structure (since displayed in ball and stick model which is not the easier for grasping the secondary structure of a molecule
    ### hard switch:
    #cmd.show("ribbon")#displays the structure as ribbon
    cmd.show("cartoon")#displays the structure as cartoon
    if (fresh): #new analysis
        cmd.color("gray30")#to improve the contrast with the different structural elements which will be displayed in color


#def chargement_fichier_pdb2(filename):
#    """
#        permet de charger le fichier pdb et de faire les pretraitements necessaires au bon fonctionnement du script dans pyMOL
#    """
#    #cmd.delete("all") #initialise l'espace de travail dans pymol
#    cmd.load(filename) #charge le fichier spÈcifiÈ par filename
#    removeAlt()#retire tous jeux de coordonnÈes alternatif pour un atome donnÈ
#    cmd.remove("hetatm")#permet de retirer les heteroatoms #attention : parfois certains fichiers sont mal annotÈs : un carbone alpha peu Ítre dÈcrit comme Ètant un hÈtÈro-atome
#    cmd.remove("h.")#permet de retirer toutes les molÈcules d'eau
#    cmd.hide("all")#permet de cacher l'ensemble de la structure
#    #cmd.show("ribbon")#rÈaffiche la structure sous le format ruban
#    cmd.show("cartoon")#rÈaffiche la structure sous le format cartoon
#    #cmd.color("gray30")#colorise l'ensemble de gris - pour faire ressortir ultÈrieurement par contraste les fragments sÈlectionnÈs

#def identificateur_pymol(atom):
def identify_pymol(atom) 
    """
        From a tuple describing one atom position 
        return a string (the id of the atom) which will allow to handle the atom inside pymol
    """
    #Note 1: named tuple weren't used since the lab was still using on more than half of the computer pyMOL with python 2.4 which doesn't handle the tuple names
    #Note 2: field[1] of the atom contains the identifier of the chain residue             
    #        field[0] contains the number of the residue
    
    #example of output:  "c;A&i;87&n;ca"
    
    return ("(c;"+atom[1]+"&i;"+atom[0]+"&n;ca)") 


## Visualisation and pymol

def color_struct(elmt,t_elmt):
    """
        color the given element in function of its nature: hence the user can easy visualy-check the selection;
        the key structures for the superposition are: SER (catalytic), B5, B4, B6, B7, B8 and B3
    """

    
    if t_elmt=="SER": #in case of serine, an additional sphere is showed to pinpoint their position
        cmd.show("spheres",identify_pymol(elmt))
    
    elif t_elmt=="SER_CAT": #for the designed serine as the catalytic one 
        cmd.color("chocolate",identify_pymol(elmt))
            
    elif t_elmt=="B5": #for the beta sheet that is referenced as the 5th on the reference model
    #it is the beta-sheet that the catalytic serine is
        for e in elmt:
            cmd.color("red",identify_pymol(e)) #all the atoms of this beta-sheet are colored in red
    
    elif t_elmt=="B6":#for the beta sheet that is referenced as the 6th on the reference model
                      #the one closest to the 5th on the C-term
        for e in elmt:
            cmd.color("orange",identify_pymol(e))
    
    elif t_elmt=="B7":#for the beta sheet that is referenced as the 7th on the reference model
                       #the one next to the 6th on the C-end - might not exist
        for e in elmt:
            cmd.color("tv_orange",identify_pymol(e))
    
    elif t_elmt=="B8":#for the beta sheet that is referenced as the 5th on the reference model
        for e in elmt:
            cmd.color("brightorange",identify_pymol(e))
            
    elif t_elmt=="B4":#for the beta sheet that is referenced as the 4th on the reference model
                      #the closest one to the 5th on the N-term
        for e in elmt:
            cmd.color("hotpink",identify_pymol(e))
    
    elif t_elmt=="B3":#for the beta sheet that is referenced as the 3rd on the reference model
                     #the closest to the 4th on the N-term - might not exist
        for e in elmt:
            cmd.color("magenta",identify_pymol(e))
    
    elif t_elmt=="B2":#for the beta sheet that is referenced as the 2nd on the reference model
                      #the closest to the 3rd on the N-term
        for e in elmt:
            cmd.color("lightpink",identify_pymol(e))
    
    elif t_elmt=="B_AUTRE2":#for the other beta-sheet in the model
        for e in elmt:
            cmd.color("skyblue",identify_pymol(e))
    #print "on sort de la fonction color_struct_beta_ser"#debug
    
    #The two other amino acids of the catalytic site are colored and stress 
    #as they give a good qualitative perception of the rate of the superposition
    elif t_elmt=="HIS":
        cmd.show("spheres",identify_pymol(elmt))
        cmd.color("yellow",identify_pymol(elmt))
    
    elif t_elmt=="ACIDE":
        cmd.show("spheres",identify_pymol(elmt))
        cmd.color("wheat",identify_pymol(elmt))
    


