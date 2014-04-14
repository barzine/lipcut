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

#importation des packages propres ‡ python, extÈrieurs 
#et dÈveloppÈs spÈcialement pour cette application

#dÈfinition des procÈdures et fonctions
def calcul_distance(list1,list2):
    """
    fonction permettant de calculer la distance 
    entre deux groupes d'atomes passes en parametre sous forme de listes
    retourne le score
    plus ce dernier est petit, plus les groupes d'atomes consideres sont proches
    """
    
    #utilisation de la structure atomCA; utilisation des champs [3], [4] et [5] correspondant au x,y et z de chaque atomeCA
    
    #passage en parametres des listes d'atomes a comparer
    #elles peuvent etre de taille differente
    #elles contiennent des informations sur des atomes : elles sont donc formattees de la meme maniere
    
    #print "dans la fonction de calcul_distance" #debug

    assert len(list1)>0 and len(list2)>0 #Arret de l'application, car liste vide passee en parametre dans la fonction calculs.calcul_distance
        
    #calcul de la  somme des distances au carre entre les deux premiers de chaque liste et les deux derniers        
    score1=((float(list1[0][3])-float(list2[0][3]))*(float(list1[0][3])-float(list2[0][3]))\
                      +((float(list1[0][4])-float(list2[0][4]))*(float(list1[0][4])-float(list2[0][4]))\
                        +((float(list1[0][5])-float(list2[0][5]))*(float(list1[0][5])-float(list2[0][5])))))
    
    score1+=((float(list1[-1][3])-float(list2[-1][3]))*(float(list1[-1][3])-float(list2[-1][3]))\
                      +((float(list1[-1][4])-float(list2[-1][4]))*(float(list1[-1][4])-float(list2[-1][4]))\
                        +((float(list1[-1][5])-float(list2[-1][5]))*(float(list1[-1][5])-float(list2[-1][5])))))
    
    ##calcul de la  somme des distances au carre entre le premier de la premiere liste et le dernier de la seconde et vice et versa            
    score2=((float(list1[-1][3])-float(list2[0][3]))*(float(list1[-1][3])-float(list2[0][3]))\
                      +((float(list1[-1][4])-float(list2[0][4]))*(float(list1[-1][4])-float(list2[0][4]))\
                        +((float(list1[-1][5])-float(list2[0][5]))*(float(list1[-1][5])-float(list2[0][5])))))
            
    score2+=((float(list1[0][3])-float(list2[-1][3]))*(float(list1[0][3])-float(list2[-1][3]))\
             +((float(list1[0][4])-float(list2[-1][4]))*(float(list1[0][4])-float(list2[-1][4]))\
               +((float(list1[0][5])-float(list2[-1][5]))*(float(list1[0][5])-float(list2[-1][5])))))
    
    if score1 < score2 :
        score=score1
        pos=1
    else : 
        score=score2
        pos=-1
    
    resultat=[score,pos]
    #print "on sort de calcul_distance" #debug
    return(resultat)        
