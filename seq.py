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
# DATE  : 2011-03
# REV   : 1
# REQUIREMENTS: 
#
#############################################################################


#from pymol import cmd
#from pymol import stored

import utils
import calculs
#import ser

#import extract
#import b_utils
#import site_utils
#import string
#import os
#from os import path


def parcours_liste(source,list_reference,liste_atomeCA,indexCA):
    """ 
        cette fonction prend un element de la liste des feuillets, ce qui permet de definir les bornes 
        entre lesquelles il faut parcourir la liste des des carbones alpha
        pour en retirer le groupe de residus qui serait le plus proche du groupe donne en reference (list_reference)
        le groupe resultat doit etre de la meme taille que le groupe donne en reference
        le groupe resultat (list_resultat) est retourne a la fin
    """
        
#recuperation des bornes devant etre parcouru dans la liste "liste_atomeCA"

    #dans un premier temps, on recupere la chaine ou les chaines concernees par le brin source 
    #if source[2]!=source[5]:
        #traitement a revoir
        #pass
    
    #else : #source[2]==source[5]
    
    
    #cles=indexCA.keys()
    #cles.sort()
    #utils.affichage_boucle(cles)
    #print("cle, valeur\n")
    #for c in cles:
    #    print indexCA[c]
    
    
    #print "source[3] : "+str(source[3])
    #print "source : "+str(source) 
    #print "source[2] : "+str(source[2])
    
    #print "source[6] : "+str(source[6])
    #print "source[5] : "+str(source[5])
    
    assert int(source[3])>=int(indexCA[source[2]][0])
    assert int(source[6])<=int(indexCA[source[5]][2])
    
    
    #deb=int(source[3])-1 #on a le numero du residu initial du feuillet en source[3]
    #fin=int(source[6]) ##on a le numero du residu final du feuillet en source[3]
    
    deb=int(indexCA[source[2]][1])-int(indexCA[source[2]][0])+int(source[3])
    fin=int(indexCA[source[2]][1])-int(indexCA[source[2]][0])+int(source[6])+1
    
    #print fin

    #print "dans la fonction de parcours_liste" #debug
    #initialisation des variables
    
    list_tmp=list()#permet de conserver le resultat pour chaque morceau (score et groupe) 
    #a la fin, le meilleur groupe sera renvoye
    
    try :
        #print "on va tester la taille de la liste" #debug
        assert(len(list_reference)<= int(fin)-int(deb)) #il faut que la liste envoye soit valide
        ok=True
    except:
        ok=False
    if ok:
        l=len(list_reference) #on conserve la taille de la liste
        #print l #debug
        #print "taille du truc a faire"+str(l) #debug
        #initialisation pour le premier calcul de score
        #print "on va initialiser la liste tmp" #debug
        for i in range(l):
            list_tmp.append(liste_atomeCA[int(deb)+i])
            #print liste_atomeCA[deb+i] #debug
        #print "check1"
        #print "initialisation de la liste tmp : ok" #debug
        #calcul du score   
        #print "on va calculer le score" #debug
        
        score=calculs.calcul_distance(list_reference,list_tmp)
        #print "score calcule sans prob" #debug
        #print "check2"
        
        #print "list_tmp"#debug
        #utils.affichage_boucle(list_tmp)#debug
        #score=calculs.calcul_distance(list_reference,list_tmp) #on calcule le score
        #print "score"
        #print score
        
        #python gerant lui-meme ses pointeurs :
        listp=utils.copy(list_tmp)
        #print "check3"
        #affectation dans list_resultat = > va servir pour la boucle qui suit
        #print "on va affecter la liste resultat" #debug
        list_resultat=[score,listp]
        #print "check4"
        #print "affectation reussi" #debug
        #parcours de la liste des atomes entre les bornes qui nous interessent
        #print "parcours de la liste" #debug
        #for i in range(deb+l,fin):
        #    print i
        
        #print fin #debug
        #print indexCA[source[2]][1] #debug
        #print int(indexCA[source[2]][0]) #debug
        #print int(source[6]) #debug
        #print indexCA[source[2]] #debug
        
        for i in range(deb+l,fin) :
            #print "check5."+str(i)
            #print "fonction xrange reussi" #debug
            list_tmp.pop(0) #on retire le premier element car on "avance" le long de la sequence
            #print "pop reussi" #debug
            list_tmp.append(liste_atomeCA[i])#on ajoute le carbone alpha qui suit dans la sequence
            #print "ajout en fin reussi" #debug
            #print liste_atomeCA[i] #debug
           
            #print "list_tmp" #debug
            #utils.affichage_boucle(list_tmp)#debug
            score=calculs.calcul_distance(list_reference,list_tmp) #on calcule le score
            #print "score"#debug
            #print score#debug
            
            #print score #debug
            if (float(score[0])<float(list_resultat[0][0])) and float(score[0])!=float(0):
                #print "dans le if"
                listp=utils.copy(list_tmp)
                #print "utils.copy effectue"
                list_resultat=[score,listp]
            #else :
                #print "tout va bien mais test final non passe" #debug
                         
    else : #la taille de la liste a comparer a la reference est plus petite que cette derniere
        #print "Warning : la taille de la liste a comparer est plus petite que celle de la reference"#+str(fin)+" source"+str(source)+" reference:"+str(list_reference)
        #print "liste_reference"#debug
        #print list_reference#debug
        
        #print "debut de la liste a comparer"#debug
        #print liste_atomeCA[deb]#debug
        #print "fin de la liste a comparer"#debug
        #print liste_atomeCA[fin]#debug
        
        list_tmp=[]
        list_tmp.append(liste_atomeCA[deb])
        list_tmp.append(liste_atomeCA[int(fin)-1])
        #print "list tmp" #debug
        #print list_tmp #debug
        score=calculs.calcul_distance(list_reference,list_tmp)
        list_resultat=[score,list_tmp]
        #readapter le calcul si le feuillet est moins long que le taille du segment choisi pour la coupe
    
    #print "on sort de parcours_liste" #debug
    #print list_resultat #debug
    return list_resultat



def recup_seq(seq_beta_prec,listB,listA,indexCA,place):
    """
        sous module pour recuperer le positionnement des feuillets les uns par rapport aux autres
        (description a revoir)
        listB : liste des brins beta n'etant pas encore positionnes par rapport aux autres
        seq_beta : sequence choisi sur le brin  precedent
        listA : liste des carbones alpha de l'ensemble des residus de la proteine considere
        indexCA : index explicitant la liste des atomes references dans listA (quel residu dans quelle chaine,
        a partir d'ou)
        place : permet de selectionner correctement les brins selon leur position (avant ou apres) "seq_beta_prec"
    """
    #print "dans recup_seq" #debug
    #print "seq_beta_prec" #debug
    #print seq_beta_prec #debug
    #initialisation de la structure permettant de conserver les differents resultats 
    resultat_L=[]
    #determination du beta accole suivant et n'etant pas encore positionne
    for b in listB :
        if place=="avant" and int(b[6])<int(seq_beta_prec[0][0]) or place=="apres" and int(b[3])>int(seq_beta_prec[-1][0]):
            #print "on fait le parcours de liste pour :"+str(b)
            r=(parcours_liste(b,seq_beta_prec,listA,indexCA),b) #on sauvegarde le resultat de parcours liste pour ce brin avec les informations du brin eux-memes
            resultat_L.append(r) #on ajoute l'ensemble a la liste resultat, on va ainsi pouvoir rechercher le brin ayant le meilleur score
            
    #on verifie en premier lieu qu'au moins un resultat a ete selectionne :
    if len(resultat_L)<1:
        return (["none"],["none"],listB,0)
    
    #else : des resultats possibles ont ete selectionnes, on verifie que la distance entre les deux listes de residus n'est pas extravagante :
    #for elmt in resultat_L:
    #    if not()):
    #        resultat_L.remove(elmt)
    
    #if len(resultat_L)<1:
    #    return (["none"],["none"],listB)

    #else: #en theorie, il ne doit pas y avoir plus d'une solution, mais au cas ou, on recupere la meilleure solution.
    #le plus petit score est le brin qui est le plus proche
    #on utilise une variable temporaire pour recuperer le tuple qui nous interresse
    tmp=min(resultat_L)#le score etant le premier element de chaque tuple, on peut utiliser la fonction min de python
    tmp,beta=tmp #on recupere les informations du brin dans beta et on ecrase tmp pour ne recupere seulement que le resultat de la fonction
    #parcours liste associe a ce brin c-a-d on recupere dans tmp[0] : le score et dans tmp[1] : la liste des atomes (carbones alpha) selectionnes sur ce brin pour l'alignement futur
    seq_beta=tmp[1]
    #remarque : il serait p.e judicieux de renvoyer le deuxieme terme de score c-a-d tmp[0][1] car il permet d'expliciter le "sens" dans la quelle la sequence doit etre lu.
    # ou alors faire un seq_beta=seq_beta.reverse() si tmp[0][1]==-1
    pos=tmp[0][1]
    #if tmp[0][1]==-1:
        #utils.affichage_boucle(seq_beta)  #debug
    #    seq_beta.reverse()
    
    #utils.affichage_boucle(listB) #debug
    #on retire ce brin de la liste initiale, car maintenant il est attribue :
    listB.remove(beta)
    #utils.affichage_boucle(listB) #debug
    
    #print beta      #debug
    #utils.affichage_boucle(seq_beta)  #debug
    
    #print "on sort de recup_seq" #debug
    return(beta,seq_beta,listB,pos)
