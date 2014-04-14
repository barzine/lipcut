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
####packages propres ‡ python:
import os
from glob import glob
import string
#####packages dÈveloppÈs spÈcialement pour l'application
import extract
import utils
import assignation
import rapport
import estimation 
import sys

####
#Corps du programme

#dÈfinition du rÈpertoire contenant les donnÈes sources :
Files = glob("../*.pdb")
#Files = glob("../base_pbm/*.pdb")
#Files = glob("../base_pbm/parcours_liste/deb/*.pdb")
#Files = glob("../base_pbm/parcours_liste/fin/*.pdb")
#Files = glob("../base_pbm/serine/*.pdb")
#Files = glob("../base_pbm/mauvais/*.pdb")
#Files = glob("../base_pbm/def_incomplete/*.pdb")
#Files = glob("../base/1HPL.pdb")
#Files=glob("../base/1CUW.pdb")
#Files = glob("../test/*.pdb")
#Files = glob("../base/2LIP.pdb")
#Files = glob("../base/1TCA.pdb")
#Files = glob("../base/1TD3.pdb")
#Files = glob("../base/1AGY.pdb")
#dÈfinition du rÈpertoire des rÈsultats :
Dir_Res="../resultats/"
#creation du dossier resultats s'il n'existe pas
if not(os.path.isdir(Dir_Res)):
    os.makedirs(Dir_Res, mode=0777)#tous les rÈpertoires n'existant pas


#intialisation des listes repertoriants les differents problemes
res1=[]
coordonnees=[]
multiple_residus=[]
trou=[]
hetatm=[]
#description avec HETATM
#probleme de definition : 1RP1
#sans serine : 1CUI et 1CUJ
#traitement 
cpt=0 #permet de numÈroter les structures traitÈs
for f in Files:
    print f
    #intialisation des dictionnaires d'avertissement
    warning={}
    w1={}
    w2={}
    w3={}
    
    #incrÈmentation du numÈro de la structure
    cpt+=1
    print "\n\nnouvelle structure n∞ "+str(cpt)+" : "+f+"\n" #debug : affichage de suivi
    
    #assignation des structures secondaires selon stride
    #filename=extract.pdb2stride2pdb(f)
    filename=string.split(os.path.basename(f),".")[0]+"_stride.pdb"
    w1['algoSS']="Stride"
    
    utils.chargement_fichier_pdb(filename)
    
    #extraction de la liste des sÈrines, des brins beta, des carbonnes alpha, de l'index des carbones, de la dÈfinition du site catalytique et
    #des avertissements ayant pu Ítre gÈnÈrÈs 
    ser_list,beta_list,atomCA_list,indexCA,site,w1['extract']=extract.extraction_data_pdb_ser_beta_site_atomeCA(filename)

    s_list,atomCA_L2,indexCA2,warning_cdc=extract.extraction_data_pdb_ser_atomeCA2(filename)
    #utils.affichage_boucle(atomCA_list)#debug : affichage de contrÙle
    #print "listB :"#debug : affichage de contrÙle
    #utils.affichage_boucle(beta_list)#debug : affichage de contrÙle
    
    resultat,w2['traitement']=assignation.traitement(ser_list,beta_list,4,atomCA_list,indexCA,site)#traitement
    resultat_commentee="Aucune solution trouvee"
    #on va calculer un score pour chaque couple serine_possiblement_catalytique et l'assignation des resultats correspondants
    #mais il est verifie qu'un tel traitement a eu lieu:
    if resultat!= ['']:
        resultat_commentee,w2['estimation']=estimation.estimation_assignation(resultat)
        
        dir_sol=rapport.creation_script_resultat(resultat_commentee,filename,Dir_Res)
        
    #if len(resultat_commentee)>0:
    #    utils.color_struct(resultat_commentee[0][0],"SER_CAT")
    #    
    #    print 
    #    print 
    #    v_color=resultat_commentee[0][1].keys()
    #    print "v_color : "+str(v_color)
    #    print resultat[0][1]
        
    #    for v in v_color:
    #        print v
    #        print resultat[0][1][v]
    #        if v!= 'B_AUTRE' and (resultat_commentee[0][1][v])!=['none']:
    #            utils.color_struct(resultat_commentee[0][1][v][1],v)
    
    
    
    
    warning=dict(w1.items()+w2.items()+w3.items())#regeneration d'un dictionnaire unique des avertissements
    print "warning :", warning["extract"], warning["algoSS"]#debug : affichage de contrÙle
    if resultat_commentee=="Aucune solution trouvee":
        print resultat_commentee
        print "Il faut analyser ce fichier : "+f+" manuelement "
        sys.exit()
        
    rapport.rapport(filename,resultat_commentee,warning,dir_sol)#gÈnÈration du rapport et des scripts permettant de visualiser les solutions
    if int(atomCA_list[0][0])!=1 :
        res1.append(filename.split('_')[0])
    
    if "Residus multiples :" in w1['extract']:
        multiple_residus.append(filename.split('_')[0])
    
    if "Coordonnees multiples" in w1["extract"]:
        coordonnees.append(filename.split('_')[0])    
    
    cles=indexCA.keys()
    cles.sort()
    num=list()
    alpha=list()
    for c in cles:
        num.append(c.split('.')[1])
        alpha.append(c.split('.')[0])
    beta=list()
    beta.append(alpha[0])
    for a in alpha:
        if a not in beta:
            beta.append(a)
    
    if len(beta)!=len(num):
        trou.append(filename.split('_')[0])
        
    if atomCA_L2!=atomCA_list:
        hetatm.append(filename.split('_')[0])  
        
    print "done for"+f+"\n"
    #partie ‡ revoir
    
    #les fichiers temporaires sont supprimer
    os.remove(filename)#nouveau fichier pdb
    os.remove(string.split(os.path.basename(f),".")[0]+".stride")#fichier stride
    
stats=open(Dir_Res+"/stats5tgl.txt","w")   

stats.write("#probleme de definition : 1RP1\n#sans serine : 1CUI et 1CUJ\n")
stats.write("\n")
stats.write("premier residu different de 1, nb :"+str(len(res1))+"\n")
stats.write(str(res1))
stats.write("\n")
stats.write("\n")
stats.write("structures possedant plus d'un jeu de coordonnees pour un meme atome, nb :"+str(len(coordonnees))+"\n")
stats.write(str(coordonnees))
stats.write("\n")
stats.write("\n")
stats.write("structures possedant plus d'un meme aa pour une meme position, nb :"+str(len(multiple_residus))+"\n")
stats.write(str(multiple_residus))
stats.write("\n")
stats.write("\n")
stats.write("structures ayant la definition de certains de leur residus en tant que heteroatomes, nb :"+str(len(hetatm))+"%\n")
stats.write(str(hetatm))
stats.write("\n")
stats.write("\n")
stats.write("structures possedant un trou dans leur definition, nb :"+str(len(trou))+"%\n")
stats.write(str(trou))



stats.close()
print "all done"
