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
#import utils
import copy

#dÈclaration des diffÈrentes fonctions et procÈdures


def feuillet_unique(liste):
    """
        Cette fonction permet de retirer tous les doublons d'une liste de brins beta donnÈe.
    """
    #En effet, la description pdb peut faire intervenir le mÍme brin plusieurs fois 
    #dans la description d'un barrilet ou pour celle d'autres structures spÈcifiques
    #
    #Un brin est considÈre comme Ètant rÈfÈrencÈ deux fois, si le rÈsidu N-terminal 
    #et le rÈsidu C-terminal du brin sont les mÍmes, c-a-d que l'on a les mÍmes numÈros, 
    #les mÍmes identificateurs de chaÓne pour les deux rÈsidus aux extrÈmitÈs du brins
    #
    #Ce traitement se fait indÈpendemment de l'identifiant des brins ou de leur orientation
    
    #VÈrification de la nÈcessitÈ de ce traitement :
    if len(liste)<2 : #la liste doit contenir au moins deux ÈlÈments
        #s'il existe moins de deux ÈlÈments dans la liste, 
        #soit c'est une liste vide, soit elle ne contient qu'un seul ÈlÈment
        #il n'existe donc pas de doublon, la liste peut Ítre retournÈe telle qu'elle
        return liste
  
    #s'il existe plus de deux ÈlÈments :  
    nv_list=[liste.pop(0)] #initialisation de la nouvelle liste ne contenant 
    #pas de doublon et qui sera retournÈe ‡ la fin vec le premier ÈlÈment de la liste d'origine 
    #note : le premier ÈlÈment de la liste est en mÍme temps retirÈ de la liste d'origine
    
    #pour chaque ÈlÈment restant de la liste d'origine, il est vÈrifiÈ qu'une mÍme description n'est pas dÈj‡ rÈfÈrencÈe
    for e in liste :
        #print "e =" #debug : affichage de contrÙle
        #print e #debug : affichage de contrÙle
        #le brin n'est recopiÈ que si le premier et le dernier rÈsidu sont diffÈrents de ceux de tous les brins dÈj‡ rÈfÈrencÈs :
        trouve=False
        i=0
        while not(trouve) and i<len(nv_list):
            if ((e[1],e[2],e[3])==(nv_list[i][1],nv_list[i][2],nv_list[i][3]) and (e[4],e[5],e[6])==(nv_list[i][4],nv_list[i][5],nv_list[i][6])): #or ((e[1],e[2],e[3])==(i[4],i[5],i[6]) and (e[4],e[5],e[6])==(i[1],i[2],i[3])) :
                trouve=True
            i+=1
        if not trouve :
            nv_list.append(e)
               
    return (nv_list)




#####

def verif_integrite(listBrin,index,listCA):
    """
        permet de verifier que pour chaque brin de la liste des brins beta, tous les atomes sont bien definis
        profite de l'occasion pour remplacer l'identifiant de la chaine par la clÈ qui lui correspond dans l'index des chaines 
        des atomes carbones alphas et ajoute un neuvieme champ ‡ la description des brins qui permet de vÈrifier de l'intÈgritÈ du brin
    """
    #initialisation des variables
    warning_brinB=""
    listR=list()
    t=tuple()
    cles=index.keys()
    #print "dans la fonction verif_integrite"#debug : affichage de contrÙle
    for b in listBrin:
        #pour chaque brin, on parcourt du carbone alpha en N-terminal, jusqu'‡ la fin o˘ jusqu'au premier atomeCA manquant
        deb=False
        fin=False
        IDdeb=""
        trouve=False
        k=0#initialisation du compteur permettant de parcourir le vecteur cle
        while (not(trouve)): #permet d'arreter le traitement des que la bonne entree dans le dictionnaire est trouvee
            if ((not deb)and b[2]==index[cles[k]][4]):
                #l'identifiant de la chaine du premier rÈsidu correspond bien ‡ celui de l'entrÈe du dictionnaire actuellement considÈrÈe
                #ainsi, il s'agit peut-Ítre de la bonne entrÈe du dictionnaire                
                if int(index[cles[k]][0])<=int(b[3])and int(index[cles[k]][2])>=int(b[3]):
                    #le numero du residu initial du brin est compris dans cette entrÈe du dictionnaire
                    trouve=True #l'entrÈe courante est bien celle contenant l'extrÈmintÈ N-terminal du brin beta qui nous interresse
                    deb=True#passe a vrai car on a trouve le debut
                    IDdeb=copy.deepcopy(cles[k])#la clÈe rÈfÈrencant le fragment concernant l'extrÈmitÈ N-terminal du brin est mÈmorisÈ
                    #IDdeb=utils.copy(cles[k])#la clÈe rÈfÈrencant le fragment concernant l'extrÈmitÈ N-terminal du brin est mÈmorisÈ
                    if index[cles[k]][4]==b[5] and int(index[cles[k]][0])<=int(b[6]) and int(index[cles[k]][2])>=int(b[6]):
                        fin =True
                    else:
                        fin=False
                #fin du if trouvant le debut
            #fin du if cherchant la bonne entree
            k+=1
        #fin du while pour les entrees du dictionnaire
        
        #crÈation du tuple pour la nouvelle description qui sera retournÈe
        if fin:#le fragment concernant l'extrÈmitÈ N-terminal et aussi celui qui concerne le c-terminal
            t=(b[0],b[1],IDdeb,b[3],b[4],IDdeb,b[6],b[7],1)
        else :#le fragment comportant l'extrÈmitÈ N-terminal ne contient pas celui du C-terminal, en d'autres termes,
            #il y a un certains nombres de rÈsidus non dÈcrits dans le fichier source
            t=(b[0],b[1],IDdeb,b[3],b[4],b[5],b[6],b[7],0)
            #avertisssement pour le rapport gÈnÈrÈ :
            warning_brinB+="\nLe brin '"+b[0]+"' est incompletement decrit dans le fichier d'origine"
            print "\nLe brin '"+b[0]+"' est incompletement decrit dans le fichier d'origine"#debug : affichage de contrÙle
        listR.append(t)
    #fin du for permettant de parcourir l'ensemble des brins beta
            
    #print "sortie de la fonction verif_integrite"#debug : affichage de contrÙle
    return [listR,warning_brinB]
