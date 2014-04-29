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

#####packages dÈveloppÈs spÈcialement pour l'applicatio
import extract
import utils
import assignation
import rapport
import estimation 

####
#Corps du programme

#dÈfinition du rÈpertoire contenant les donnÈes sources :
Files = glob("../base/*.pdb")
#Files = glob("../base_pbm/*.pdb")
#Files = glob("../base_pbm/parcours_liste/deb/*.pdb")
#Files = glob("../base_pbm/parcours_liste/fin/*.pdb")
#Files = glob("../base_pbm/serine/*.pdb")
#Files = glob("../base_pbm/mauvais/*.pdb")
#Files = glob("../base_pbm/def_incomplete/*.pdb")
#Files = glob("../base/1HPL.pdb")
#Files=glob("../base/1CUW.pdb")
#Files = glob("../test/*.pdb")
Files = glob("../base/2LIP.pdb")
#Files = glob("../base/1TCA.pdb")
#Files = glob("../base/1TD3.pdb")
#dÈfinition du rÈpertoire des rÈsultats :
Dir_Res="../resultats/"
#creation du dossier resultats s'il n'existe pas
if not(os.path.isdir(Dir_Res)):
    os.makedirs(Dir_Res, mode=0777)#tous les rÈpertoires n'existant pas

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
    filename=extract.pdb2stride2pdb(f)
    w1['algoSS']="Stride"
    
    utils.chargement_fichier_pdb(filename)
    
    #extraction de la liste des sÈrines, des brins beta, des carbonnes alpha, de l'index des carbones, de la dÈfinition du site catalytique et
    #des avertissements ayant pu Ítre gÈnÈrÈs 
    ser_list,beta_list,atomCA_list,indexCA,site,w1['extract']=extract.extraction_data_pdb_ser_beta_site_atomeCA(filename)

    #utils.affichage_boucle(atomCA_list)#debug : affichage de contrÙle
    #print "listB :"#debug : affichage de contrÙle
    #utils.affichage_boucle(beta_list)#debug : affichage de contrÙle
    
    resultat,w2['traitement']=assignation.traitement(ser_list,beta_list,4,atomCA_list,indexCA,site)#traitement
    
    
    #on va calculer un score pour chaque couple serine_possiblement_catalytique et l'assignation des resultats correspondants
    #mais il est verifie qu'un tel traitement a eu lieu:
    if len(resultat)>0:
        resultat_commentee,w2['estimation']=estimation.estimation_assignation(resultat)
        
        rapport.creation_script_resultat(resultat_commentee,filename,Dir_Res)
        
    if len(resultat_commentee)>0:
        utils.color_struct(resultat_commentee[0][0],"SER_CAT")
    #    
        print 
        print 
        v_color=resultat_commentee[0][1].keys()
        print "v_color : "+str(v_color)
        print resultat[0][1]
        
        for v in v_color:
            print v
            print resultat[0][1][v]
            if v!= 'B_AUTRE' and (resultat_commentee[0][1][v])!=['none']:
                utils.color_struct(resultat_commentee[0][1][v][1],v)
    
    
    
    
    warning=dict(w1.items()+w2.items()+w3.items())#regeneration d'un dictionnaire unique des avertissements
    print "warning :", warning["extract"], warning["algoSS"]#debug : affichage de contrÙle
    
    #rapport.rapport(filename,resultat_commentee,warning)#gÈnÈration du rapport et des scripts permettant de visualiser les solutions
    
    print "done for"+f+"\n"
print "all done"
