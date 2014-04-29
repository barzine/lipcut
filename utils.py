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
# REV   : 3 (2014)
# REQUIREMENTS: 
#
#############################################################################

#Checking functions and procedures: allow to print the different structures used 
#in the application

#def affichage_boucle(Liste):
def print_loop(Liste):
    """
   Print the content of a liste, line per line
    """ 
    for elmt in Liste :
        print elmt 
          
##
#def affichage_boucle_liste_dans_list(Liste):
def print_loop_in_loop(Liste):
    """
    Print line per line the content of nested lists (inside a list) 
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

#def affichage_dico(dico):
def dico_print(dico):
    """
        print the content of a dictionary in a more human-lisible way
    """
    
    ring=dico.keys()
    print "key, valeur"
    for k in ring:
        print k,dico[k] 



#Copy function ... before I learned about copy and the magic deep.copy function 
#overide the dynamic allocation by a static one
def copy(sourceList,verbose=False):
    """
        created before being aware of the deep.copy [copy] function
        return a new object instead of a synonym
    """
   
    newList=[] #initialisation

    for i in sourceList: #for each i element of sourceList
        if (verbose):
            print i #for debugging purpose
        newList.append(i) #ajout de l'ÈlÈment

    return(newList) #the new list is returned to the main environment, hence creating the new object


