#!python
# -*-coding:Latin-1 -*
##############################################################################
#
# @SUMMARY: -- 
#             
# Compatibility : tested only on PyMOL 0.99  - python v2.4
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
# DATE  : 2011-03
# REV   : 1
# REQUIREMENTS: 
#
#############################################################################

#importation des packages propres ‡ python, extÈrieurs 
#et dÈveloppÈs spÈcialement pour cette application
import copy
import site_utils
#import utils

#dÈfinition des fonctions et procÈdures
def redefinition(ser_list,listCA):
    """
        retrieve a specific list (of possible catalitic serine) indexing all the serines
        and return that list while the index of the serines is replaced by their description
        in the listCA
    """

    #initialisation of the results list to be return
    #print "inside redefinition (serine)" #debug
    list_res=list()
    for s in ser_list:
        #print s #debug
        t=(listCA[s[0][0]],s[1],s[2])
        list_res.append(t)
     
    #print "out of redefinition (ser)" #debug
    return list_res


def serine_site_catal(ser_list,beta_list,atom_L,site_d,opt=True):
    """
        allows to fund the more likely catalytic serine
        optionnaly later: could modulate the processing based on user choices if needed 
    """
    
    #initialisation of different lists
    ser_site_p=list()
    ser_site_p0=list() #initialisation of potential catalytic site serines list -
    # research only done at the C-terminal extremity and at the end of description of the current beta 
    ser_site_p1plus=list() #initialisation of potential catalytic site serines list -
    # research only done at the (1+ C-terminal extremity) of the end of description of the current beta
    ser_site_p1moins=list() #initialisation de la liste des serines potentielles du site catalytique, on ne les recherche
    #qu'en -1 c-terminal  de la description du brin beta considere
    ser_site_p2plus=list() #initialisation de la liste des serines potentielles du site catalytique, on ne les recherche
    #qu'en +2 c-terminal  de la description du brin beta considere
    
    #print "dans serine_site_catal" #debug
    #utils.affichage_boucle(ser_list) #debug
    for e in ser_list :#e est donc un tuple renvoyant a une serine
        #print "serine=",str(e) #debug : affichage de contrÙle
        cpt=-1 #initialisation d'un compteur, permet de connaitre le tour de boucle : correspond au numero du brin concerne dans la liste des brins
        for b in beta_list :
            #print b #debug : affichage de contrÙle
            cpt+=1
            #on va traiter les differentes positions possibles de la serine en fonction du brin 
            #les positions retenues sont (en C-terminal du brin) 
            #0, +1, -1 et +2
            #les serines seront listes dans l'ordre de leur sequence mais regroupees en fonction de leur position relative
            #vis a vis du brin considere comme etant alors le b5
            #La chaine est verifiee en meme temps que le positionnement de la serine et de la sequence des feuillets
            
            #Le champ contenant le numero du residu de l'atome considere est situe dans le champ [0] dans la liste des atomes 
            #et son identifiant de chaine est dans le champ [1]
            #Les serines sont stockees selon le format [0] :index dans la liste des atomCA et [1] numero de residu dans la chaine
            #b etant la description d'un brin beta, b[6] est le numero du dernier residu concerne par le brin considere
            #et b[5] l'identifiant de la chaine (donne par le fichier pdb) de ce meme residu
            #le troisieme terme permet de calculer la borne a partir de laquelle il faudra recuperer les aa
            #de maniere ‡ ce que les residus soient DANS le feuillet et en N-terminal de la serine
            #et cela tout en excluant la serine
            if int(atom_L[int(e[0])][0])==(int(b[6])) and (atom_L[int(e[0])][1])==(b[5]).split('.')[0]: 
            #if int(atom_L[int(e[0])][0])==(int(b[6])): #si extraction par calcul sur les coordonnees des atomes CA
                t=(e,cpt,1)#enregistrement de la position de la serine dans la liste des atomes et la position du brin beta accolee dans la liste des beta
                #le troisieme terme "1" permet d'indiquer ‡ partir de quelle position il faut recuperer les aa
                #print t #debug
                ser_site_p0.append(t)
               
            elif int(atom_L[int(e[0])][0])==(int(b[6])+1)and (atom_L[int(e[0])][1])==(b[5]).split('.')[0]:
                t=(e,cpt,1)#enregistrement de la position de la serine dans la liste des atomes et la position du brin beta accolee dans la liste des beta
                #print t #debug
                ser_site_p1plus.append(t) 
               
            elif int(atom_L[int(e[0])][0])==(int(b[6])-1)and (atom_L[int(e[0])][1])==(b[5]).split('.')[0]:
                t=(e,cpt,1)#enregistrement de la position de la serine dans la liste des atomes et la position du brin beta accolee dans la liste des beta
                #print t #debug
                #le troisieme terme "1" permet d'indiquer ‡ partir de quelle position il faut recuperer les aa
                ser_site_p1moins.append(t)
            
            elif int(atom_L[int(e[0])][0])==(int(b[6])+2)and (atom_L[int(e[0])][1])==(b[5]).split('.')[0]:
                t=(e,cpt,2)#enregistrement de la position de la serine dans la liste des atomes et la position du brin beta accolee dans la liste des beta
                #print t #debug
                ser_site_p2plus.append(t)      
        #fin de la boucle for       
        
    #print "recomposition de la liste ser_site_p"#debug
    #on ne prend que les listes non vides
    if len(ser_site_p0)>0:
        ser_site_p=ser_site_p+ser_site_p0
        
    if len(ser_site_p1plus)>0:
        ser_site_p=ser_site_p+ser_site_p1plus
        
    if len(ser_site_p1moins)>0:
        ser_site_p=ser_site_p+ser_site_p1moins
        
    if len(ser_site_p2plus)>0:
        ser_site_p=ser_site_p+ser_site_p2plus
    
    
    if(len(ser_site_p)<1): #la liste des serines pouvant etre catalytiques est vide !
        list_tmp=list()
        #toutes les serines sont rÈÈtudiÈes afin d'en retirer celles qui sont le mieux placÈes par rapport aux brins beta
        for s in ser_list :#e est donc un tuple renvoyant a une serine
            cpt=-1 #initialisation d'un compteur, permet de connaitre le tour de boucle : correspond au numero du brin concerne dans la liste des brins
            for b in beta_list :
                cpt+=1
                #l'identificateur de la chaine de la serine est comparÈe ‡ celle du brin
                if (atom_L[int(s[0])][1])==(b[5]).split('.')[0] :
                    t=(int(atom_L[int(s[0])][0])-int(b[6]),cpt,s)
                    list_tmp.append(t)

        #la valeur (absolue) minimale de distance entre les serines et l'extrÈmitÈ C-terminal du brin est rÈcupÈrÈe
        if len(list_tmp)>0:
            #initialisation avec la premiere serine et le premier brin
            mini=abs(list_tmp[0][0])
            for elmt in list_tmp:
                if abs(elmt[0])<mini:
                    mini=abs(elmt[0])
        #les serines correspondant ‡ cette distance minimale sont conservÈes pour l'assignation des brins selon la topologie de rÈfÈrence
        for elmt in list_tmp:
            if abs(elmt[0])==mini:
                if elmt[0]>0:#la serine est positionnÈe aprËs le brin beta
                    ser_site_p.append((elmt[2],elmt[1],elmt[0]))
                else : #la serine est positionnÈe dans le brin beta ou avant
                    if atom_L[int(elmt[2][0])][0]>int(b[3]):
                        ser_site_p.append((elmt[2],elmt[1],1))
    
    if site_d!="" : #une description du site catalytique a ete donnee, on va croiser cette information avec la liste des serines preselectionnees
        #les serines concernees vont etre placees en tete de liste
        ser_filtre=site_utils.verif_site(ser_site_p,beta_list,site_d,atom_L)
        
        for s in ser_filtre:#on a donc s etant un tuple renvoyant a une serine
            #on va devoir placer les elements en tete de liste
            ser_site_p.remove(s)
        
        #on replace ces serines en tete de liste :
        ser_site_p=ser_filtre+copy.deepcopy(ser_site_p)
        
    if ser_site_p>0:#des serines ont ÈtÈ trouvÈes selon les diffÈrents critËres spÈcifiÈs
        #les "pointeurs" vers les serines sont remplacÈs par la description des sÈrines
        #print "toujours dans serine site catal"#debug
        #print ser_site_p #debug
        ser_site_p=redefinition(ser_site_p,atom_L)
        #print "ser site" #debug
        #utils.affichage_boucle(ser_site_p) #debug
    else: #aucune serine n'a ÈtÈ trouvÈe
        ser_site_p=["none"]
    #print "on sort de serine_site_catal" #debug
    return ser_site_p
