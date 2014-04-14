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

import os
import string
import subprocess


#dÈfinitions des structures utilisÈes:
    #format de stockage d'un atomeCA:
    #structure : AtomCA
        #premier champ [0] : numero du residu concerne
        #deuxiËme champ [1] : id de la chaine ou se situe le residu
        #troisiËme champ [2] : index dans le dictionnaire concernant ce residu (forme de l'id de la chaine plus numero) 
        #quatriËme champ [3] : abscisse de l'atome
        #cinquiËme champ [4] : ordonne de l'atome
        #sixiËme champ [5] : hauteur de l'atome
        #septiËme champ [6] : nom du residu

    #format de stockage d'une entrÈe du dictionnaire indexCA
    #structure : entrÈe d'indexCA
        #clÈ : formÈe de l'identifiant de la chaine +"."+valeur d'un compteur auto-incrÈmentÈe
        #valeur :
            #premier champ [0] : numÈro du premier rÈsidu rÈfÈrencÈ par cette clÈ
            #deuxiËme champ [1] : position (index) de ce rÈsidu dans la liste globale des atomesCA
            #troisiËme champ [2] : numÈro du dernier rÈsidu rÈfÈrencÈ par cette clÈ
            #quatriËme champ [3] : position (index) de ce rÈsidu dans la liste globale des atomesCA
            #cinquiËme champ [4] : identifiant de la chaine o˘ se situent les atomesCA concernÈes par cette entrÈe
   
    #format des brins beta :
    #structure : BrinB
            #champ[0] : identificateur du brin considere
            #champ[1] : nom du premier aa dans la sequence du brin considere
            #champ[2] : chaine sur laquelle est le premier residu du brin considere
            #champ[3] : numero dans la sequence du premier residu du brin considere
            #champ[4] : nom du dernier aa dans la sequence du brin considere
            #champ[5] : chaine sur laquelle est le dernier residu du brin considere
            #champ[6] : numero dans la sequence du dernier residu du brin considere
            #champ[7] : sens de ce brin par rapport au precedent         
    


#dÈfinition des procÈdures et des fonctions
def extraction_data_pdb_ser_atomeCA(filename) :
    """
    permet de recuperer les coordonnees des carbones alphas des residus decrits dans un fichier pdb 
    on ne recupere qu'un seul jeu de coordonnees, meme s'il existe des alternatives
    on en profite pour recupere la liste des serines aussi
    """

    #note : utilisation de la structure : AtomCA
    
    #print "extraction_data_pdb_atomeCA" #debug : affichage de contrÙle
    
    atomCA_L=[]#liste des carbones alpha des rÈsidus
    s_list=[]#va contenir la liste des sÈrines
    res=-1#compte le nombre de rÈsidus dans la sequence (il faut faire +1) - sert aussi d'index 
    
    indexCA={}#Ce dictionnaire sert d'index permettant ainsi de conserver la position d'un carbone alpha
    #dans la liste globale les rÈfÈrenÁant. Ce dictionnaire permet de pallier aux problËmes gÈnÈrÈs par l'existence
    #de plusieurs chaÓnes et par ceux provoquÈs par l'absence de dÈfinition de certains carbones alpha.
    
    #initialisation pour le dictionnaire indexCA des diffÈrentes variables permettant de construire les diffÈrentes entrÈes.
    #et pour faciliter la gÈnÈration de cet index, une premiËre entrÈe, vouÈe ‡ Ítre supprimÈe, est crÈÈe.
    chaineID="Deleted" #variable servant ‡ gÈnÈrer les clÈs d'indexCA 
    numResInit="" #initialisation du numÈro de rÈsidu intial (premier champ)
    indexResInit="" #position du premier rÈsidu dans la liste des atomesCA (deuxiËme champ)
    numResFin="" # mÈmorise le numÈro du dernier rÈsidu concernÈe par cette entrÈe (troisiËme champ)
    indexResFin=""#position de ce dernier rÈsidu dans la liste des atomesCA (quatriËme champ)
    refID="" #identifiant initial de la chaine comportant les rÈsidus prÈcÈdents
    cpt=-1 #permet de gÈnÈrer des clÈs diffÈrentes en cas de donnÈes manquantes, mÍme si les rÈsidus proviennent de la mÍme chaine 

    warning_cdc="" #permet de conserver l'ensemble des warnings gÈnÈrÈs par l'extraction, en vue d'un rapport ultÈrieuer
    warning_desc_atom="" #diffÈrent de "" s'il existe plus d'une description de coordonnÈes pour un mÍme carbone alpha
    warning_res_mult="" #diffÈrent de "" s'il existe plus d'un rÈsidu dÈcrit pour la mÍme position - c-a-d, mÍme chaine, mÍme numÈro
    #corps de la fonction :
    file=open(filename,'r')#ouverture du fichier => ne doit pas poser de probleme car deja verifie avant 
    for line in file: #le fichier est parcouru ligne par ligne
        #le premier mot de la ligne sert de parseur pour le traitement
        if line[0:4]=="ATOM" or line[0:6]=="HETATM" :   
            if line[13:15]=="CA" and line[77]=='C': #recupere que si l'atome courant est le carbone alpha
                doublon=False
                if line[16]==' ' or line[16]=='A': #permet de recuperer seulement les coordonnes une fois meme s'il existe une position alternative
                    res+=1 #index du nouveau rÈsidu dans la liste totale des carbones alpha
                    #vÈrification de "l'unicitÈ" de l'aa dÈcrit pour cette position
                    #en  effet dans certains fichiers de la PDB, pour une mÍme position on peut avoir deux aa qui soient dÈcrits
                    #la Serine est rÈcupÈrÈe de maniËre privilÈgiÈe s'il y a lieu, sinon la premiËre description enregistrÈe est conservÈe
                    if line[26]!=' ':#il existe donc plusieurs descriptions pour cette position dans la sequence
                        #la nouvelle ligne est comparÈe au dernier enregistrement
                        if (int(atomCA_L[-1][0])==int(line[22:26]) and atomCA_L[-1][1]==line[21]):
                            #il existe donc pour cette position une description dÈj‡ enregistrÈe
                            #L'existence de differents aa pour une mÍme position est conservÈe en vue d'un avertissement ultÈrieure
                            #print "Il existe plusieurs possibilites (plusieurs aa) pour une meme position dans une meme chaine"
                            warning_res_mult="Residus multiples : Il existe plusieurs possibilites (plusieurs aa)\n pour une meme position dans une meme chaine\n"
                            #a noter, il est preferable de refaire l'affectation a chaque fois au lieu de tester l'existence de cet avertissement
                            #avant de faire l'affectation (on gagne en temps de calcul)
                            res-=1 #l'index de la liste des atomes doit Ítre corrigÈ, car sur les deux descriptions, une seule sera conservÈe
                            #la Serine Ètant prÈfÈrentiellement rÈcupÈrÈe,l'enregistrement ne sera modifiÈ que si la description actuelle en rÈfËre une 
                            if line[17:20]=="SER":
                                #le dernier enregistrement est retire
                                atomCA_L.remove(atomCA_L[-1])
                            else: #ce n'est pas la description d'une serine et il existe dÈj‡ un enregistrement pour cette position
                                #cette description ne doit pas etre recupere
                                doublon=True
                        #else:
                            #on ne fait rien car aucune description n'a encore ÈtÈ enregistrÈe pour ce couple numero de residu, chaine 
                            #pass
                    
                    #creation d'un vecteur qui va permettre de conserver les index des differents residus, cela permet de palier
                    #aux problemes pouvant survenir lorsque la description pdb n'est pas complete (chaines multiples, sequences decrites incompletes)
                    #ce qui suit permet de gÈrer les discontinuitÈs dans la liste des carbonnes alpha 
                    if (not doublon and((refID!= line[21]) or (int(line[22:26])!=int(numResFin)+1))):
                        indexCA[chaineID]=(numResInit,indexResInit,numResFin,indexResFin,refID)
                        cpt+=1
                        #le residu actuel est le premier de la nouvelle chaine
                        chaineID=line[21]+"."+str(cpt)
                        #print chaineID #debug
                        #son numero et son index sont mÈmorisÈs
                        numResInit=line[22:26]
                        indexResInit=int(res)
                    #end if
                    
                    #ce qui suit est vrai pour tout carbone alpha devant Ítre rÈcupÈrÈ, c'est a dire qu'il n'a pas ÈtÈ considÈrÈ comme un doublon
                    if not(doublon):
                        refID=line[21]
                        t=(str(int(line[22:26])),refID,chaineID,str(float(line[31:38])),str(float(line[39:46])),str(float(line[47:54])),line[17:20])
                        atomCA_L.append(t)
                        numResFin=line[22:26]
                        indexResFin=int(res)
                    #pourrait etre judicieux de revoir ce point (de ne rÈcupÈrer que l'index seulement)
                    if "SER" in line : #recupere que si l'atome courant appartient a une serine
                        t=(res,line[22:26])
                        s_list.append(t)#on ajoute le numero de ce tuple dans la liste des atomeCA a la liste des serines
                    #END if
                else : #il existe d'autres positions alternatives pour certains atomes
                    warning_desc_atom="Coordonnees multiples : Il existe dans le fichier "+filename+" des atomes \npossedant plus d'un jeu de coordonnees\n"
                #END if
            #end if     
        #END if
    #END for 
    file.close() 
    #le dernier maillon est rajoutÈ ‡ l'indexCA
    indexCA[chaineID]=(numResInit,indexResInit,numResFin,indexResFin,refID)
    del indexCA["Deleted"]

    #print indexCA #debug : affichage de contrÙle
    #print "atom_list" #debug : affichage de contrÙle
    #utils.affichage_boucle(atomCA_L)#debug : affichage de contrÙle
    if len(atomCA_L)<1: #aucune descriptio, de carbone alpha n'a ÈtÈ extraite du fichier d'entrÈe
        warning_cdc+="DonnÈes introuvables : Les coordonnÈes d'aucun carbone alpha n'a pu Ítre rÈcupÈrÈs ‡ partir \ndu fichier "+filename+" .\n"
    if len(s_list)<1: #aucune description de sÈrine n'a pu Ítre extraite du fichier d'entrÈe
        warning_cdc+="SÈrine introuvable : Aucune sÈrine n'a pu Ítre extraite ‡ partir du fichier "+filename+" .\n"
    if len(indexCA)>1: #existence de plusieurs clÈs, deux cas possibles : soit plusieurs chaines dÈcrites, soit description des carbonnes alpha incomplËte
        #existence de "trous" dans la description totale
        warning_cdc+="Fragmentation de l'information : Il existe plusieurs chaÓnes dÈcrites ou bien il manque la description \nde certains carbones alpha.\n"
    #gÈnÈration d'une chaine de caracteres contenant tous les avertissements ayant ÈtÈ gÈnÈrÈs lors de l'extraction
    if warning_res_mult!="":
        warning_cdc+=warning_res_mult
    if warning_desc_atom!="":
        warning_cdc+=warning_desc_atom
    
    #print "nombre residus :" #debug
    #print res  #debug
    list_resultat=[s_list,atomCA_L,indexCA,warning_cdc]
    return(list_resultat)
 
###
def extraction_data_pdb_ser_beta_site_atomeCA(filename) :
    """
    permet de recuperer la liste des serines et celles des feuillets beta, si fichier pdb correctement formate
    A REVOIR : les feuillets beta et l'identification de la chaine sur laquelle les residus sont postionnees
    """
    
    #note : utilisation des structures BrinB, indexCA et atomeCA
    #print "extraction_data_pdb_ser_beta_site_atomeCA" #debug
    #declaration des variables locales :
    b_list=[]#va contenir la liste des brins beta
    site="" #description du site catal. dans le fichier pdb
    
    #corps de la fonction :
    file=open(filename,'r')#ouverture du fichier => ne doit pas poser de probleme car deja verifie avant 
    for line in file: #on parcourt le fichier
    #on "parse" sur le premier mot de la ligne
        if line[0:5]=="SHEET" :
            #construction de l'identifiant (3 premiers champs de la ligne apres le terme SHEET)
            #note : separation des 3 champs pour permetre de les retravailler facilement si necessaire par la suite
            #en utilisation la fonction split et en renseigant : pour separateur                   
            id=str(int(line[5:10]))+":" #ajout d'un separateur ":" apres le numero du brin dans le feuillet considere
            p=line[11:14] #identificateur du feuillet considere - utilisation d'une variable temporaire
            #pour permettre le nettoyage des espaces superflus avant d'etre ajouter a l'identifiant
            id +=p.strip()+":" #nettoyage des espaces superflus et du separateur
            id+=str(int(line[15:16])) #ajout du nombre total de brin dans le feuillet considere
            #creation d'un tuple pour l'element beta courant : identifiant du brin, residu initial, chaine du residu initial ,n∞dans la sÈquence du rÈsidu initial, 
            #residu final, chaine du residu final, position du residu final dans la sequence, sens (parallele =1 et anti-parallËle = -1) par rapport au brin precedent
            t=(id,line[17:20],line[21],str(int(line[22:26])),line[28:31],line[32],str(int(line[33:37])),str(int(line[37:40])))
            #on ajoute ce brin a la liste des brins
            b_list.append(t)
            #END if
        #END if   
        if line[0:4]=="SITE":
            site=site+"#"+line
    #les deux listes s_list et b_list
    file.close()
    #intialisation de la chaine de caractËres qui contiendra tous les avertissements dont l'utilisateur doit prendre connaissance
    warning=""
    #ajout de l'avertissements s'il y a lieu :
    if len(b_list)<1:
        warning+="Brin beta introuvable : Aucune dÈfinition de brin beta n'a pu Ítre extraite du fichier "+filename+"\n"
    
    #extraction de la liste des sÈrines, des carbones alpha, de l'index des carbones alpha. Si des avertissements sont gÈnÈrÈs, ils vont
    #Ítre rÈcupÈrÈs afin de pouvoir Ítre portÈs ‡ l'attention de l'utilisateur
    s_list,atomCA_L,indexCA,warning=extraction_data_pdb_ser_atomeCA(filename)
    
    #ajout de la remarque s'il y a lieu    
    if site=="":
        warning+="Remarque  sur le site catalytique : aucun site catalytique n'est dÈcrit dans le fichier "+filename+"\n"
     
    #affichage de contrÙle divers :
    #print "nombre residus :" #debug : affichage de contrÙle
    #print res  #debug : affichage de contrÙle
    #print "s_list :" #debug : affichage de contrÙle
    #utils.affichage_boucle(s_list)#debug : affichage de contrÙle
    #print "b_list:"#debug : affichage de contrÙle
    #utils.affichage_boucle(b_list)#debug : affichage de contrÙle
    #print "atom_list"#debug : affichage de contrÙle
    #utils.affichage_boucle(atomCA_L)#debug : affichage de contrÙle
    #print"site :" #debug : affichage de contrÙle
    #print site #debug : affichage de contrÙle
    
    #un seul element ne peut Ítre renvoyÈ avec une fonction, 
    #une liste des rÈsultats va donc Ítre crÈÈe et va Ítre retourner 
    list_resultat=[s_list,b_list,atomCA_L,indexCA,site,warning]
    return(list_resultat) #donc dans list_resultat[0] s_list et dans list_resultat[1] : b_list

###
def parcoursStride(f_res,file_source):
    """
        parcours le fichier stride pour ajouter la definition des helices et des brins beta
    """
    
    #print "dans parcoursStride" #debug :affichage de contrÙle
    try: #essai d'ouverture du fichier en mode "ajout"
        res_file=open(f_res,"a")
    except IOError:
        print "Probleme d'ouverture du fichier"+f_res
        return     

    try:
        file_stride=open(file_source,"r")
    except IOError:
        f_res.close()#fermeture du premier fichier
        print "Probleme d'ouverture du fichier"+file_source
        return
    #si le programme s'execute, c'est que l'ouverture des deux fichiers c'est rÈalisÈ sans problËme
    nbHelix=0
    nbSheet=0
    for line in file_stride:
        if line.startswith("LOC"):
            tmp=line.split()
            if tmp[1]=="AlphaHelix":
                nbHelix+=1
                res_file.write('HELIX %(nb)4d %(nb)3d %(resInit)s %(resInitChain)s %(resNumInit)4d  %(resFin)s %(resFinChain)s %(resNumFin)4d  1%(lg)36d\n' % {'nb' : int(nbHelix),'resInit':tmp[2],'resInitChain':tmp[4],'resNumInit':int(tmp[3]),'resFin':tmp[5], 'resNumFin':int(tmp[6]), 'resFinChain':tmp[7],'lg':int(tmp[6])-int(tmp[3])+1})
            if tmp[1]=="310Helix":
                nbHelix+=1
                res_file.write('HELIX %(nb)4d %(nb)3d %(resInit)s %(resInitChain)s %(resNumInit)4d  %(resFin)s %(resFinChain)s %(resNumFin)4d  5%(lg)36d\n' % {'nb' : int(nbHelix),'resInit':tmp[2],'resInitChain':tmp[4],'resNumInit':int(tmp[3]),'resFin':tmp[5], 'resNumFin':int(tmp[6]), 'resFinChain':tmp[7],'lg':int(tmp[6])-int(tmp[3])+1})
            if tmp[1]=="Strand": 
                nbSheet+=1
                res_file.write('SHEET %(nb)4d %(nb)3d 0 %(resInit)s %(resInitChain) s%(resNumInit)4d  %(resFin)s %(resFinChain)s%(resNumFin)4d  2\n' % {'nb' : int(nbSheet),'resInit':tmp[2],'resInitChain':tmp[4],'resNumInit':int(tmp[3]),'resFin':tmp[5], 'resNumFin':int(tmp[6]), 'resFinChain':tmp[7]})
        #else: #les definitions ont deje ete donnees, on passe a la boucle suivante
            #pass
    #fin du for
    res_file.close()
    file_stride.close()
    #print "sortie de parcours stride"#debug :affichage de contrÙle  
    

###
def stride2pdb(filename):    
    """
    permet de recreer un fichier pdb en utilisant l'assignation de 
    structures donnee par l'algorithme STRIDE a partir du fichier fichier.pdb et fichier.stride
    """
    
    #filename : nom du fichier pdb a modifier
 
    #creation d'un fichier temporaire permettant de sauvegarder le contenu du nouveau fichier pdb
    try:    
        file_res=open(string.split(os.path.basename(filename),".")[0]+"_stride.pdb","w")
    except IOError:
        print "Probleme d'ouverture du fichier"+file_res
        return
    
    flag=False #permet de connaitre si les descriptions des brins beta ont ÈtÈ ajoutÈs ‡ la fin du fichier ou non
    try:
        file_pdb=open(filename,"r")
    except IOError:
        print "Probleme d'ouverture du fichier"+filename
        file_res.close()
        return     
    
    #on doit recopier le debut du fichier pdb, puis introduire les definitions donnees par le fichier stride 
    #et enfin completer le fichier pdb avec le reste des donnees
    #premier =True #permet de savoir si on est deja sorti et revenu sur le fichier pdb
    credit=True #permet d'indiquer si les remarques sur la generation des structures secondaires ont ete ajoutes
    premier=True
    for line in file_pdb :
        #print line
        #on gere en premier lieu, les specificites:
        if (line[0:5]!="SHEET" and line[0:5]!="HELIX" and line[0:3]!="END"):
            if line[0:6]=="REMARK" and credit:
                file_res.write('REMARK   0 : Assignation des structures HELIX et SHEET\n')
                file_res.write('REMARK   0 : par l algorithme STRIDE\n')
                file_res.write('REMARK   0 : Generation de ce fichier %(*)s_stride.pdb par stride2pdb version 0.1\n' %{'*':string.split(os.path.basename(filename),".")[0]})
                credit=False
            file_res.write(line)
        else :#sinon on rencontre soit l'entete de ligne HELIX soit SHEET soit la fin du fichier 
            if premier:
                if line[0:3]=="END":
                    flag=True
                premier=False
                file_res.close()
                parcoursStride(string.split(os.path.basename(filename),".")[0]+"_stride.pdb",string.split(os.path.basename(filename),".")[0]+".stride")
                file_res=open(string.split(os.path.basename(filename),".")[0]+"_stride.pdb","a")
    file_pdb.close()
    
    if premier: #la description pour une raison ou une autre n'a pas ÈtÈ transcripte
        parcoursStride(string.split(os.path.basename(filename),".")[0]+"_stride.pdb",string.split(os.path.basename(filename),".")[0]+".stride")

    if not flag:
        file_res.write('END')
        
    file_res.close()
    return(string.split(os.path.basename(filename),".")[0]+"_stride.pdb")

#####
def pdb2stride2pdb(filename):
    """
    applique l'algorithme stride a un fichier pdb et recree ensuite un nouveau fichier pdb a partir
    de l'assignation de structures secondaires (helice et brin)
    retourne le nom du nouveau fichier
    """
    
    #creation du fichier .stride a partir du fichier pdb donne en parametre
    subprocess.call("Stride -o "+filename+" >"+string.split(os.path.basename(filename),".")[0]+".stride",shell=True)
    #creation du nouveau fichier pdb
    filename=stride2pdb(filename)
    #le nom du fichier nouvellement crÈer est renvoyÈ
    return filename

########




def extraction_data_pdb_ser_atomeCA2(filename) :
    """
    permet de recuperer les coordonnees des carbones alphas des residus decrits dans un fichier pdb 
    on ne recupere qu'un seul jeu de coordonnees, meme s'il existe des alternatives
    on en profite pour recupere la liste des serines aussi
    """

    #note : utilisation de la structure : AtomCA
    
    #print "extraction_data_pdb_atomeCA" #debug : affichage de contrÙle
    
    atomCA_L=[]#liste des carbones alpha des rÈsidus
    s_list=[]#va contenir la liste des sÈrines
    res=-1#compte le nombre de rÈsidus dans la sequence (il faut faire +1) - sert aussi d'index 
    
    indexCA={}#Ce dictionnaire sert d'index permettant ainsi de conserver la position d'un carbone alpha
    #dans la liste globale les rÈfÈrenÁant. Ce dictionnaire permet de pallier aux problËmes gÈnÈrÈs par l'existence
    #de plusieurs chaÓnes et par ceux provoquÈs par l'absence de dÈfinition de certains carbones alpha.
    
    #initialisation pour le dictionnaire indexCA des diffÈrentes variables permettant de construire les diffÈrentes entrÈes.
    #et pour faciliter la gÈnÈration de cet index, une premiËre entrÈe, vouÈe ‡ Ítre supprimÈe, est crÈÈe.
    chaineID="Deleted" #variable servant ‡ gÈnÈrer les clÈs d'indexCA 
    numResInit="" #initialisation du numÈro de rÈsidu intial (premier champ)
    indexResInit="" #position du premier rÈsidu dans la liste des atomesCA (deuxiËme champ)
    numResFin="" # mÈmorise le numÈro du dernier rÈsidu concernÈe par cette entrÈe (troisiËme champ)
    indexResFin=""#position de ce dernier rÈsidu dans la liste des atomesCA (quatriËme champ)
    refID="" #identifiant initial de la chaine comportant les rÈsidus prÈcÈdents
    cpt=-1 #permet de gÈnÈrer des clÈs diffÈrentes en cas de donnÈes manquantes, mÍme si les rÈsidus proviennent de la mÍme chaine 

    warning_cdc="" #permet de conserver l'ensemble des warnings gÈnÈrÈs par l'extraction, en vue d'un rapport ultÈrieuer
    warning_desc_atom="" #diffÈrent de "" s'il existe plus d'une description de coordonnÈes pour un mÍme carbone alpha
    warning_res_mult="" #diffÈrent de "" s'il existe plus d'un rÈsidu dÈcrit pour la mÍme position - c-a-d, mÍme chaine, mÍme numÈro
    #corps de la fonction :
    file=open(filename,'r')#ouverture du fichier => ne doit pas poser de probleme car deja verifie avant 
    for line in file: #le fichier est parcouru ligne par ligne
        #le premier mot de la ligne sert de parseur pour le traitement
        if line[0:4]=="ATOM" :   
            if line[13:15]=="CA" and line[77]=='C': #recupere que si l'atome courant est le carbone alpha
                doublon=False
                if line[16]==' ' or line[16]=='A': #permet de recuperer seulement les coordonnes une fois meme s'il existe une position alternative
                    res+=1 #index du nouveau rÈsidu dans la liste totale des carbones alpha
                    #vÈrification de "l'unicitÈ" de l'aa dÈcrit pour cette position
                    #en  effet dans certains fichiers de la PDB, pour une mÍme position on peut avoir deux aa qui soient dÈcrits
                    #la Serine est rÈcupÈrÈe de maniËre privilÈgiÈe s'il y a lieu, sinon la premiËre description enregistrÈe est conservÈe
                    if line[26]!=' ':#il existe donc plusieurs descriptions pour cette position dans la sequence
                        #la nouvelle ligne est comparÈe au dernier enregistrement
                        if (int(atomCA_L[-1][0])==int(line[22:26]) and atomCA_L[-1][1]==line[21]):
                            #il existe donc pour cette position une description dÈj‡ enregistrÈe
                            #L'existence de differents aa pour une mÍme position est conservÈe en vue d'un avertissement ultÈrieure
                            #print "Il existe plusieurs possibilites (plusieurs aa) pour une meme position dans une meme chaine"
                            warning_res_mult="Residus multiples : Il existe plusieurs possibilites (plusieurs aa)\n pour une meme position dans une meme chaine\n"
                            #a noter, il est preferable de refaire l'affectation a chaque fois au lieu de tester l'existence de cet avertissement
                            #avant de faire l'affectation (on gagne en temps de calcul)
                            res-=1 #l'index de la liste des atomes doit Ítre corrigÈ, car sur les deux descriptions, une seule sera conservÈe
                            #la Serine Ètant prÈfÈrentiellement rÈcupÈrÈe,l'enregistrement ne sera modifiÈ que si la description actuelle en rÈfËre une 
                            if line[17:20]=="SER":
                                #le dernier enregistrement est retire
                                atomCA_L.remove(atomCA_L[-1])
                            else: #ce n'est pas la description d'une serine et il existe dÈj‡ un enregistrement pour cette position
                                #cette description ne doit pas etre recupere
                                doublon=True
                        #else:
                            #on ne fait rien car aucune description n'a encore ÈtÈ enregistrÈe pour ce couple numero de residu, chaine 
                            #pass
                    
                    #creation d'un vecteur qui va permettre de conserver les index des differents residus, cela permet de palier
                    #aux problemes pouvant survenir lorsque la description pdb n'est pas complete (chaines multiples, sequences decrites incompletes)
                    #ce qui suit permet de gÈrer les discontinuitÈs dans la liste des carbonnes alpha 
                    if (not doublon and((refID!= line[21]) or (int(line[22:26])!=int(numResFin)+1))):
                        indexCA[chaineID]=(numResInit,indexResInit,numResFin,indexResFin,refID)
                        cpt+=1
                        #le residu actuel est le premier de la nouvelle chaine
                        chaineID=line[21]+"."+str(cpt)
                        #print chaineID #debug
                        #son numero et son index sont mÈmorisÈs
                        numResInit=line[22:26]
                        indexResInit=int(res)
                    #end if
                    
                    #ce qui suit est vrai pour tout carbone alpha devant Ítre rÈcupÈrÈ, c'est a dire qu'il n'a pas ÈtÈ considÈrÈ comme un doublon
                    if not(doublon):
                        refID=line[21]
                        t=(str(int(line[22:26])),refID,chaineID,str(float(line[31:38])),str(float(line[39:46])),str(float(line[47:54])),line[17:20])
                        atomCA_L.append(t)
                        numResFin=line[22:26]
                        indexResFin=int(res)
                    #pourrait etre judicieux de revoir ce point (de ne rÈcupÈrer que l'index seulement)
                    if "SER" in line : #recupere que si l'atome courant appartient a une serine
                        t=(res,line[22:26])
                        s_list.append(t)#on ajoute le numero de ce tuple dans la liste des atomeCA a la liste des serines
                    #END if
                else : #il existe d'autres positions alternatives pour certains atomes
                    warning_desc_atom="Coordonnees multiples : Il existe dans le fichier "+filename+" des atomes \npossedant plus d'un jeu de coordonnees\n"
                #END if
            #end if     
        #END if
    #END for 
    file.close() 
    #le dernier maillon est rajoutÈ ‡ l'indexCA
    indexCA[chaineID]=(numResInit,indexResInit,numResFin,indexResFin,refID)
    del indexCA["Deleted"]

    #print indexCA #debug : affichage de contrÙle
    #print "atom_list" #debug : affichage de contrÙle
    #utils.affichage_boucle(atomCA_L)#debug : affichage de contrÙle
    if len(atomCA_L)<1: #aucune descriptio, de carbone alpha n'a ÈtÈ extraite du fichier d'entrÈe
        warning_cdc+="DonnÈes introuvables : Les coordonnÈes d'aucun carbone alpha n'a pu Ítre rÈcupÈrÈs ‡ partir \ndu fichier "+filename+" .\n"
    if len(s_list)<1: #aucune description de sÈrine n'a pu Ítre extraite du fichier d'entrÈe
        warning_cdc+="SÈrine introuvable : Aucune sÈrine n'a pu Ítre extraite ‡ partir du fichier "+filename+" .\n"
    if len(indexCA)>1: #existence de plusieurs clÈs, deux cas possibles : soit plusieurs chaines dÈcrites, soit description des carbonnes alpha incomplËte
        #existence de "trous" dans la description totale
        warning_cdc+="Fragmentation de l'information : Il existe plusieurs chaÓnes dÈcrites ou bien il manque la description \nde certains carbones alpha.\n"
    #gÈnÈration d'une chaine de caracteres contenant tous les avertissements ayant ÈtÈ gÈnÈrÈs lors de l'extraction
    if warning_res_mult!="":
        warning_cdc+=warning_res_mult
    if warning_desc_atom!="":
        warning_cdc+=warning_desc_atom
    
    #print "nombre residus :" #debug
    #print res  #debug
    list_resultat=[s_list,atomCA_L,indexCA,warning_cdc]
    return(list_resultat)
