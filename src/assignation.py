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

##libraries and personal files importation
import copy
import b_utils
import utils
import ser
import seq




#def positionnement_feuillets(num_feuillet,ser_cat,pos,listB,lg,listA,indexCA,site):
def sheets_position(sheet_num,ser_cat,pos,listB,lg,listA,indexCA, site):
    """
        from the serine list and the beta sheet list, a possible spatial organisation of the beta sheet is build,
        i.e. the beta sheet 5 (beta5), 4 (beta4), 6 (beta6), 7 (beta7) and 8 (beta8) are assigned
        sheet_num:
        ser_cat:
        pos:
        listB: list of the beta sheet defined from the pdb file
        lg: length of sequence to use for the beta sheet
        indexCA:
        site
        
    """
    #initialisation 
    beta={}#dict for the beta sheet assignation 
    avert={}#to keep record of the different errors or warnings
    
    #assignation of beta5
    beta5=listB[sheet_num]#when the current ser has been assigned, the number of the beta sheet used for the selection has been recorded as well
    listB.remove(beta5)#every assigned beta sheet is remove from listB which only contains then still-to-be assigned beta sheet
    
    #The following steps take advantage of the described topology for the lipases:
    # a) selection of the n residues of beta5 - this done by recording the alpha carbons of the main chain
    # b) lg residues at the C-term of this serine are considered to constitute the beta5

    seq_beta5=list()#list recording the Calpha of the presume beta5
    #Theses carbons are used for assigning the other beta sheets in accordance with the topology
    # and they are used in the inter-structures superpositions as well.

    
    ref=int(indexCA[ser_cat[2]][1])-int(indexCA[ser_cat[2]][0])-int(pos)+int(ser_cat[0])#to find out the beta5 we use the serint
    #the serine is just after the beta5 in the proteic sequence
    
    #the correct length of beta sheet (lg) has to be retrieved 
    for i in range(lg) :
        t=listA[ref-i] #tuple referecing the atom placed at -i from the considered serine 
        seq_beta5.append(t)#the tuple is added to the beta5 list
    #all the sequence has been recorded
    seq_beta5.reverse()#so the residues would be in the conventional order: N-term to C-term
   
    beta["B5"]=(beta5,seq_beta5)#the sequence with its name as a key are recorded in the assigned dictionary
                
    ####assignation of the other beta sheets
    beta6,seq_beta6,listB,direction=seq.retrieve_seq(seq_beta5,listB,listA,indexCA,"after")
    if seq_beta6!=['none']:
        beta["B6"]=(beta6,seq_beta6)
        if direction==-1:
            avert["B6"]="Warning: b6 is antiparallel to b5"
    else :
        beta["B6"]=['none']

    beta4,seq_beta4,listB,direction=seq.retrieve_seq(seq_beta5,listB,listA,indexCA,"before")
    if seq_beta4!=['none']:
        beta["B4"]=(beta4,seq_beta4)
        if direction==-1:
            avert["B4"]="Warning: b4 is antiparallel to b5"
    else :
        beta["B4"]=['none']
                    
    #for beta7 :
    if len(listB)>0 and seq_beta6!=['none']:
        beta7,seq_beta7,listB,direction=seq.retrieve_seq(seq_beta6,listB,listA,indexCA,"after")
        if seq_beta7!=['none'] :
            beta["B7"]=(beta7,seq_beta7)
            if direction==-1:
                avert["B7"]="Warning: b7 is antiparallel to b6"
        else :
            beta["B7"]=['none']
    else:
        beta["B7"]=['none']
        
    #for beta8 :
    if len(listB)>0 and beta["B7"]!=['none']:
        beta8,seq_beta8,listB,direction=seq.retrieve_seq(seq_beta7,listB,listA,indexCA,"after")                 
        if seq_beta8!=['none'] :
            beta["B8"]=(beta8,seq_beta8)
            if direction==-1:
                avert["B8"]="Warning: b8 is antiparallel to b7"
        else:
            beta["B8"]=['none']
    else : 
        beta["B8"]=['none']
                    
    #for beta3 :    
    if len(listB)>0 and seq_beta4!=['none']:
        beta3,seq_beta3,listB,direction=seq.retrieve_seq(seq_beta4,listB,listA,indexCA,"after")
        if seq_beta3!=['none'] :
            beta["B3"]=(beta3,seq_beta3)
            if direction==-1:
                avert["B3"]="Warning: b3 is antiparallel to b4"
        else :
            beta["B3"]=['none'] 
    else:
        beta["B3"]=['none']   
                             
    #for beta2 :
    if len(listB)>0 and beta["B3"]!=['none']:
        beta2,seq_beta2,listB,direction=seq.retrieve_seq(seq_beta3,listB,listA,indexCA,"before")
        if seq_beta2!=['none'] :
            beta["B2"]=(beta2,seq_beta2)
            if direction ==-1:
                avert["antiparallelity"]=True
        else :
            beta["B2"]=['none']
    else : 
        beta["B2"]=['none']    
        
    
    #all leftover beta sheet are regrouped in one entry of the dictionary beta
    if len(listB)>0:
        beta["B_OTHER"]=listB
        avert["B_OTHER"]=str(len(listB))+" has/have not been assigned"
    else:
        beta["B_OTHER"]=['none']
    
    if not ("antiparallelity" in avert):
        avert["antiparallelity"]=False
    
    
    if beta['B8']==['none'] or beta['B2']==['none'] :
        avert["missing"]="Warning: Some of the model beta sheets cannot been assigned in the crystal"
    else:
        avert["missing"]="Every beta sheet in the model has found a correspondance in the crystal"
        
    
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
