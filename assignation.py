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

#importation of libraries
import copy

#importation of local files
import b_utils
import utils
import ser
import seq
#import estimation


#former name: positionnement_feuillets(num_feuillet,ser_cat,pos,listB,lg,listA,indexCA,site):
def sheet_locus(sheet_nb,ser_cat,pos,listB,lg,listA,indexCA,site):
    """ a partir de la liste des serines et des feuillets beta, reconstruction de l'organisation des feuillets entre eux
        c'est-a-dire que l'on positionne les feuillets beta5, beta4, beta6, beta7 et beta3
    """

    #print "dans positionnement_feuillets" #debug : affichage de controle 
    
    #print "ser cat : "+str(ser_cat) #debug : affichage de controle 
    #print indexCA[ser_cat[2]] #debug  : affichage de controle 
    #print "pos",pos #debug  : affichage de controle 
    
    
    #initialisation des variables
    beta={}#dictionnaire permettant de mÈmoriser l'assignation des diffÈrents brins
    avert={}#permet d'enregistrer les diffÈrents problËmes pouvant subvenir 
    #assignation de beta5
    beta5=listB[sheet_nb]#lors de la sÈlection de la serine courante, le numÈro du brin beta ayant servi ‡ la sÈlectionner a ÈtÈ mÈmorisÈ aussi
    listB.remove(beta5)#le brin assignÈ est retirÈ de l'ensemble, seuls les brins non encore assignÈs sont conservÈs
        
    #les traitements qui suivent tirent avantage de la topologie dÈcrite pour ce genre de protÈines
    #sÈlection des n rÈsidus(utilisation des carbones alphas) du "beta5"
    #les lg rÈsidus sont rÈcupÈrÈs en C-terminal du brin considÈrÈ comme Ètant le brin numÈro b5 pour cette sÈrine
   
    seq_beta5=list()#crÈation de la liste qui va mÈmoriser les carbones alpha d'intÈrÍt
    #ces carbones vont permettre de positionner les autres brins et dans un second temps 
    #pour les alignements inter-structures

    #en utilisant la position de la serine par rapport au feuillet
    ref=int(indexCA[ser_cat[2]][1])-int(indexCA[ser_cat[2]][0])-int(pos)+int(ser_cat[0])
    #le dernier -1 est la car on ne veut pas recuperer la sÈrine elle-mÍme
    #print "ref : "+str(ref) #debug : affichage de controle 

    #print "boucle de recup pour seq_beta5 :"   #debug : affichage de controle 
    for i in range(lg) :
        #print ref-i#debug : affichage de controle 
        t=listA[ref-i]
        #print t   #debug : affichage de controle 
        seq_beta5.append(t)#le rÈsidu  est ajoutÈ ‡ la liste 
    
    #‡ ce niveau, la sÈquence est complËte

    #print seq_beta5#debug : affichage de controle          
    seq_beta5.reverse()#les rÈsidus sont repositionner dans l'ordre N-terminal vers C-terminal
    #utils.affichage_boucle(seq_beta5)
    beta["B5"]=(beta5,seq_beta5)#l'identificateur dui brin et la sÈquence associÈe est enregistrÈe dans le dictionnaire des brins assignÈs
                
    ####assignation des autres brins beta
    beta6,seq_beta6,listB,sens=seq.recup_seq(seq_beta5,listB,listA,indexCA,"apres")
    if seq_beta6!=['none']:
        beta["B6"]=(beta6,seq_beta6)
        if sens==-1:
            avert["B6"]="Remarque : le brin b6 est antiparallele au b5"
    else :
        beta["B6"]=['none']

    beta4,seq_beta4,listB,sens=seq.recup_seq(seq_beta5,listB,listA,indexCA,"avant")
    if seq_beta4!=['none']:
        beta["B4"]=(beta4,seq_beta4)
        if sens==-1:
            avert["B4"]="Remarque : le brin b4 est antiparallele au b5"
    else :
        beta["B4"]=['none']
                    
    #pour le beta7 :
    if len(listB)>0 and seq_beta6!=['none']:
        beta7,seq_beta7,listB,sens=seq.recup_seq(seq_beta6,listB,listA,indexCA,"apres")
        if seq_beta7!=['none'] :
            beta["B7"]=(beta7,seq_beta7)
            if sens==-1:
                avert["B7"]="Remarque : le brin b7 est antiparallele au b6"
        else :
            beta["B7"]=['none']
    else:
        beta["B7"]=['none']
        
    #pour le beta8 :
    if len(listB)>0 and beta["B7"]!=['none']:
        beta8,seq_beta8,listB,sens=seq.recup_seq(seq_beta7,listB,listA,indexCA,"apres")                 
        if seq_beta8!=['none'] :
            beta["B8"]=(beta8,seq_beta8)
            if sens==-1:
                avert["B8"]="Remarque : le brin b8 est antiparallele au b7"
        else:
            beta["B8"]=['none']
    else : 
        beta["B8"]=['none']
                    
    #pour le beta3 :    
    if len(listB)>0 and seq_beta4!=['none']:
        beta3,seq_beta3,listB,sens=seq.recup_seq(seq_beta4,listB,listA,indexCA,"apres")
        if seq_beta3!=['none'] :
            beta["B3"]=(beta3,seq_beta3)
            if sens==-1:
                avert["B3"]="Remarque : le brin b3 est antiparallele au b4"
        else :
            beta["B3"]=['none'] 
    else:
        beta["B3"]=['none']   
                             
    #pour le beta2 :
    if len(listB)>0 and beta["B3"]!=['none']:
        beta2,seq_beta2,listB,sens=seq.recup_seq(seq_beta3,listB,listA,indexCA,"avant")
        if seq_beta2!=['none'] :
            beta["B2"]=(beta2,seq_beta2)
            if sens ==-1:
                avert["antiparallelite"]=True
        else :
            beta["B2"]=['none']
    else : 
        beta["B2"]=['none']    
        
    
    #tous les brins restants sont assignÈs dans une entrÈe spÈciale du dictionnaire
    if len(listB)>0:
        beta["B_AUTRE"]=listB
        avert["B_AUTRE"]="Il existe "+str(len(listB))+" brins n'ayant pas ete assignes"
    else:
        beta["B_AUTRE"]=['none']
    
    if not ("antiparallelite" in avert):
        avert["antiparallelite"]=False
    #print "on sort de positionnement_feuillets" #debug : affichage de controle 
    
    if beta['B8']==['none'] or beta['B2']==['none'] :
        avert["manquant"]="Attention : Il existe des brins du modele qui n'ont aucune correspondance dans la structure analysee"
    else:
        avert["manquant"]="Tous les brins du modele ont trouvÈs une correspondance dans la structure analysee"
        
    
    return [ser_cat,beta,avert]




####
def traitement(listS,listB,lg,listA,indexCA,site):
    """
        va traiter toutes les serines potentiellement catalytique et va essayer d'organiser les brins beta en suivant
        la topologie de reference
    """
    
    #print "dans traitement"#debug : affichage de contrÙle
    
    #print "beta list :"
    #utils.affichage_boucle(listB)
    
    #initialisation des variables
    res=list()#liste contenant l'ensemble des solutions pour la structure traitÈe
    warning={}#dictionnaire des avertissements gÈnÈrÈs par le traitement

    listB=b_utils.feuillet_unique(listB)#ce traitement permet de retirer tout brin etant reference 2 fois. 
    #un brin est considere comme etant reference 2 fois, si le residu en N terminal et le residu en C terminal sont les memes
    #(c-a-d meme numero, meme chaine pour les deux residus aux extremites)
    #independemment de l'identifiant des brins ou de leur orientation
    #print "beta list  apres feuillet unique:"
    #utils.affichage_boucle(listB)
    
    
    #verification de l'intÈgritÈ des brins beta rÈfÈrencÈs et correction de l'indexation par rapport ‡ la liste des carbones alpha
    listB,warning["brinB"]=b_utils.verif_integrite(listB,indexCA,listA)
    
    ser_cat=list() #va permettre de conserver la liste des serines possiblement catalytique.
    #le traitement devra etre applique pour chaque element de cette liste
    #utils.affichage_boucle(ser_list)#debug : affichage de contrÙle
    #rÈcuperation de la sÈrine catalytique croisement de la liste de sÈrines avec celle de la liste des feuillets dans la variable l_ser_filtre
    l_ser_filtre=ser.serine_site_catal(listS,listB,listA,site)#idealement, il n'y a qu'une seule serine dans la liste l_ser_filtre
    #print "l_ser_filtre" #debug : affichage de contrÙle
    #utils.affichage_boucle(l_ser_filtre) #debug : affichage de contrÙle
    #print "nbre de serine selectionnees :"+str(len(l_ser_filtre))
    if len(l_ser_filtre)<1:
        warning["Ser"]="#Probleme de reference : il n'existe aucune serine repondant aux criteres de selection pour une serine possiblement catalytique."
        #on sort de la fonction car tout le reste du traitement depend de cette selection
        return [[""],warning]
    
    for s in l_ser_filtre:
        #print "s"+str(s)#debug : affichage de controle     
        ser_cat,num_feuillet,pos=s #on separe les differentes informations
        #print "ser_cat :" #debug : affichage de contrÙle
        #print ser_cat #debug : affichage de contrÙle
        #afin de faciliter le controle visuel par l'utilisateur, on fait ressortir cette serine
        utils.color_struct(ser_cat,"SER")
            
        #print "debut de la disposition des brins pour la serine considere %(ser)s" %{'ser':ser_cat}#debug      
        #l'assignation des brins pour cette sÈrine est conservÈe, avec l'identification de la sÈrine pour pouvoir Ítre utilisÈ ultÈrieurement
        sol=positionnement_feuillets(num_feuillet,ser_cat,pos,copy.deepcopy(listB),lg,listA,indexCA,site)
        print"sol"+str(sol)
        res.append(sol)#le triplet courant sÈrine/assignation des brins/avertissementsest ajoutÈ ‡ la liste des rÈsultats
        print"\n\n\n"
        #print sol #debug : affichage de controle     
    #fin du for permettant de parcourir toutes les sÈrines sÈlectionnÈes
  
    #print "sortie de traitement"
    utils.affichage_boucle(res)
    #print "nb de resultats"
    #print len(res)
    
    return [res,warning]
