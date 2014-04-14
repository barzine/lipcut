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
#import utils
#import assignation
#import rapport
#import estimation 

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
#Files = glob("../base/2LIP.pdb")
#Files = glob("../base/1TCA.pdb")
#Files = glob("../base/1TD3.pdb")
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
    
#description avec HETATM
#probleme de definition : 1RP1
#sans serine : 1CUI et 1CUJ
    

#traitement 
cpt=0 #permet de numÈroter les structures traitÈs
for f in Files:

    w1=dict()
    #incrÈmentation du numÈro de la structure
    cpt+=1
    print "\n\nnouvelle structure n∞ "+str(cpt)+" : "+f+"\n" #debug : affichage de suivi
    
    #assignation des structures secondaires selon stride
    filename=extract.pdb2stride2pdb(f)
     
    #extraction de la liste des sÈrines, des brins beta, des carbonnes alpha, de l'index des carbones, de la dÈfinition du site catalytique et
    #des avertissements ayant pu Ítre gÈnÈrÈs 
    ser_list,beta_list,atomCA_list,indexCA,site,w1['extract']=extract.extraction_data_pdb_ser_beta_site_atomeCA(filename)

    if atomCA_list[0][0]!=1 :
        res1.append(f)
    
    if "Residus multiples :" in w1['extract']:
        multiple_residus.append(filename.split('_')[0])
    
    if "Coordonnees multiples :" in w1['extract']:
        multiple_residus.append(filename.split('_')[0])    
    
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

    print "done for"+filename.split('_')[0]+"\n"
    
    
         
    print res1
    print coordonnees
    print multiple_residus
    print trou
print "all done"
