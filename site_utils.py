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

import utils


def verif_site(liste_ser,beta,site_desc,atom_L):
    """
        Croise les informations entre la liste des serines possiblement catalytiques et 
        la description donne par le fichier pdb du site catalytique
    """
    #print "dans verif_site" #debug : affichage de contrÙle
    list_ser_filtre=[]#permet de conserver les informations relatives aux sÈrines dÈcrites comme catalytiques - 
    #c-a-d qu'il existe une dÈfinition ("SITE") o˘ elles seraient situÈes
    site=site_desc.split("#") #permet de rÈcupÈrer les diffÈrentes lignes sous forme d'ÈlÈments dans une liste
    for s in site: #pour chaque ÈlÈment stockÈ dans la variable site
        #print s #debug : affichage de contrÙle
        for e in liste_ser : #pour chaque sÈrine considÈrÈe comme Ètant possiblement catalytique 
            #print atom_L[int(e[0][0])] #debug : affichage de contrÙle
            if atom_L[int(e[0][0])][0] in s :
                if e not in list_ser_filtre:
                    list_ser_filtre.append(e)
    
    #print "ser filtre" #debug : affichage de contrÙle
    #print list_ser_filtre #debug : affichage de contrÙle                   
    #print "on sort de verif_site" #debug : affichage de contrÙle
    return list_ser_filtre




def verif_histidine(list_aa,d):
    """
        permet de rechercher une histidine possiblement catalytique dans la liste d'aa envoyee
    """
    possibilite=list()
    avertissement=""
    #for aa in list_aa:
        #if aa[2]=='HIS':
            #preliste.append(aa)
    
    #on verifie si cette histidine est situe apres le dernier brin :
    if d["B8"]!=['none']:
        #il y a une description de brin pour le brin 8 :
        #l'histidine doit se situer en C-terminal de ce brin :
        #print aa[0] #debug
        #print d["B8"][0]#debug
        for aa in list_aa:
            if aa[2]=='HIS':
                if int(d["B8"][0][6])<int(aa[0]):
                    #avertissement+=""
                    possibilite.append(aa) 
                    
    elif d["B7"]!=['none']:
        #il y a une description de brin pour le brin 7 :
        #l'histidine doit se situer en C-terminal de ce brin, puisqu'il n'y a pas de definition pour le brin 7 :
        #print aa[0] #debug
        #print d["B7"][0]#debug
        for aa in list_aa:
            if aa[2]=='HIS':
                if int(d["B7"][0][6])<int(aa[0]):
                    avertissement="La ou les histidines donnees sont situees apres le brin B7, le brin B8 n'etant pas attribue dans la solution courante"
                    possibilite.append(aa)
                    
    elif d["B6"]!=['none']:
        #il y a une description de brin pour le brin 6 :
        #l'histidine doit se situer en C-terminal de ce brin :
        #print aa[0]#debug
        #print d["B6"][0]#debug
        for aa in list_aa:
            if aa[2]=='HIS':
                if int(d["B6"][0][6])<int(aa[0]):
                    avertissement="La ou les histidines donnees sont situes aprËs le brin B6, les brins B7 et B8 n'etant pas attribues dans la solution courante"
                    possibilite.append(aa)                                       
    
    if len(possibilite)>1:
        avertissement+="- Remarque : il y a plus d'une solution"
                    
    if len(possibilite)>0:
        #print len(possibilite)
        #print possibilite
        #on a au moins une histidine qui repond a nos criteres pour etre catalytique
        for his in possibilite:
            utils.color_struct(his,"HIS")    
        return [possibilite,avertissement]
    else :
        #on informe de la non existence d'histidine pour cette solution 
        return [['none'],"Aucune histidine n'a ete trouvee pour la solution courante"]
    
    
def verif_acide(list_aa,d):
    """
        permet de rechercher un aa acide possiblement catalytique dans la liste d'aa envoyee
    """
    possibilite=list()
    avertissement=""
    #on verifie si l'acide trouvÈ est situÈ aprËs l'avant-dernier brin :
    if d["B8"]!=['none']:
        for aa in list_aa:
            if aa[2] in ['ASP','GLU'] :
                #il y a une description de brin pour le brin 8 :
                #il y en a donc une aussi du brin 7 
                #pour Èviter les problemes de rÈsidus multiples:
                if aa[0].isdigit():
                    if int(d["B7"][0][6])<int(aa[0]):
                        possibilite.append(aa)
                else:
                    bb=(aa[0][0:-1],aa[1],aa[2])
                    if int(d["B7"][0][6])<int(bb[0]):
                        possibilite.append(bb)
    elif d["B7"]!=['none']:
        for aa in list_aa:
            if aa[2] in ['ASP','GLU'] :        
                #il y a une description de brin pour le brin 7 :
                #l'histidine doit se situer en C-terminal de ce brin, puisqu'il n'y a pas de definition pour le brin 7 :
                #print aa[0] #debug
                #print d["B7"][0]#debug
                avertissement="La ou les acides donnes sont situes apres le brin B6, le brin B8 n'etant pas attribue dans la solution courante"
                if aa[0].isdigit():
                    if int(d["B6"][0][6])<int(aa[0]):
                        possibilite.append(aa)
                else:
                    bb=(aa[0][0:-1],aa[1],aa[2])
                    if int(d["B6"][0][6])<int(bb[0]):
                        possibilite.append(bb)
    elif d["B6"]!=['none']:
        for aa in list_aa:
            if aa[2] in ['ASP','GLU'] :        
                #il y a une description de brin pour le brin 6 :
                #l'histidine doit se situer en C-terminal de ce brin :
                #print aa[0]#debug
                #print d["B6"][0]#debug
                avertissement="La ou les acides donnes sont situes apres le brin B5, le brin B7 n'etant pas attribue dans la solution courante"
                if aa[0].isdigit():
                    if int(d["B5"][0][6])<int(aa[0]):
                        possibilite.append(aa)
                else:
                    bb=(aa[0][0:-1],aa[1],aa[2])
                    if int(d["B5"][0][6])<int(bb[0]):
                        possibilite.append(bb)
    
    if len(possibilite)>1:
        avertissement+="- Remarque : il y a plus d'une solution"
    
    if len(possibilite)>0:
        #print possibilite
        #on a au moins une histidine qui repond a nos criteres pour etre catalytique
        for acide in possibilite:
            utils.color_struct(acide,"ACIDE")    
        return [possibilite,avertissement]
    else :
        #on informe de la non existence d'histidine pour cette solution 
        return [['none'],"Aucun acide amine acide n'a ete trouve pour la solution courante"]
    


