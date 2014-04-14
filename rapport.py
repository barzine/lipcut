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
import string
import datetime
import os
import utils

def formatageWarningSolution(dico):
    """
        prend en entrÈe le dictionnaire des avertissements et remarques liÈs ‡ une solution
        et renvoie l'ensemble sous forme de chaines de caracteres (plus ou moins bien) formattÈe
    """

    warning_formate=""
    #les differents warnings associees aux solutions
    #sont: 
    #brin 2 et 3 anti // ?
    warning_formate+="  -En ce qui concerne l'assignation des brins selon le modele de reference : \n"
    if not(dico['antiparallelite']):
        warning_formate+="#    L'antiparallelite des brins B2 et B3 n'a pu etre prouve\n"
        
    #existence de brins non assignes
    if "B_AUTRE" in dico:
        warning_formate+="#    "+dico["B_AUTRE"]+"\n" 
    #brin manquant par rapport a la topologie de reference
    if "manquant" in dico:
        warning_formate+="#    "+dico["manquant"]+"\n"
    if "B3" in dico :
        warning_formate+="#    "+dico["B3"]+"\n"
    if "B4" in dico :
        warning_formate+="#    "+dico["B4"]+"\n"
    if "B6" in dico :
        warning_formate+="#    "+dico["B6"]+"\n"
    if "B7" in dico :
        warning_formate+="#    "+dico["B7"]+"\n"        
    if "B8" in dico :
        warning_formate+="#    "+dico["B8"]+"\n"    
        
    #site catal :
    warning_formate+="  -En ce qui concerne le site catalytique : \n" 
    warning_formate+="#    La distance de selection autour de la serine est de "+str(dico["site"]["dist"])+" angtroms\n"
    #His:
    if dico["site"]["his"]!="":
        warning_formate+="#      -Pour l'histidine : "+dico["site"]["his"]+"\n"
    #Acide :
    if dico["site"]["acide"]!="":
        warning_formate+="#      -Pour l'acide amine acide : "+dico["site"]["acide"]+"\n"
    #comparaison avec le champ site de la description pdb
    #verification sur le site catal : His ? Acide?
    
    
    return warning_formate
    #return(str(dico))


def creation_script_resultat(solution,filename,repertoire):
    """
       crÈe les scripts nÈcessaires ‡ la visualisation des rÈsultats
    """

    dir_res_sol=repertoire#+"/"+string.split(os.path.basename(filename),".")[0] #si on veut chaque structure dans un repertoire propre
    if not(os.path.isdir(dir_res_sol)):
        os.makedirs(dir_res_sol, mode=0777)
        
    #le fichier pdb est recopie dans ce repertoire :
    file=open(filename,'r')
    filec=open(dir_res_sol+"/"+filename,"w")
    for line in file:
        filec.write(line)
    filec.close()
    file.close()
    
    #initialisation pour la date:
    jour=datetime.date.today()
      
    for cpt,sol in enumerate(solution) :
        #creation d'un nouveau fichier script pour la solution courante
        #si ce fichier existait deja, il sera recree
        #creation du nom de sortie peut etre a revoir
        file_res=open(dir_res_sol+"/"+string.split(filename,".")[0]+"_sol_"+str(cpt+1)+".py","w")
        #file_res=open(dir_res_sol+"/"+filename+"_spol_"+str(cpt+1)+".py","w")
        #creation de l'entete du fichier
        file_res.write("#!python \n")
        file_res.write("# -*-coding:Latin-1 -* \n")
        file_res.write("def removeAlt(obj='(all)', keep='A'):\n")
        file_res.write('    """\n')
        file_res.write("        from PyMOLWIKI\n")
        file_res.write("        removeAlt -- remove all alternate location-atoms not of altloc keep from object.\n")
        file_res.write("\n")
        file_res.write("        input:\n")
        file_res.write("              obj -- the object(s) to remove the atoms from\n")
        file_res.write("              keep -- which type of alt loc to keep\n")
        file_res.write("\n")
        file_res.write("        output: none -- removes atoms\n")
        file_res.write("\n")
        file_res.write("     examples:\n")
        file_res.write("              removeAlt # remove all altLocations that aren t altloc A\n")
        file_res.write("              removeAlt pdbID, C  # remove all but C altlocations from pdbID\n")
        file_res.write('    """\n')
        file_res.write("    #select & remove all non A altlocs\n")
        file_res.write('    remStr = "%s and not (alt ""+%s)" % (obj, keep);\n')
        file_res.write("    cmd.remove(remStr);\n")
    
        file_res.write("##########################################################################################\n")
        file_res.write("# Rapport gÈnÈrÈ de maniËre automatique ‡ partir de l'application RESATOLC                \n")
        file_res.write("# application dÈveloppÈe lors du stage de fin d'Ètudes de la formation M2 BioInfo(Nantes) \n")
        file_res.write("# stage effectuÈe auprËs du Professeur Vinh TRAN au sein de l'U3B, unitÈ 6204 du CNRS     \n")
        file_res.write("# auteur: M. BARZINE (barzine@gmail.com)                                                  \n")
        file_res.write("# crÈdit: M. BARZINE - V. TRAN                                                            \n")
        file_res.write("##########################################################################################\n")
        file_res.write("# RESATOLC : script de visualisation                                                \n")    
        file_res.write("# Date de crÈation (AAAA/MM/JJ): "+jour.strftime("%Y/%m/%d")+"                                               \n")
        file_res.write("##########################################################################################\n")
        file_res.write("##########################################################################################\n")
        file_res.write("##########################################################################################\n")
        file_res.write("### script pour :         ######################################################\n")
        file_res.write("# PDB_id="+string.split(os.path.basename(filename),".")[0]+"                                \n")
        file_res.write("# ID solution : "+str(cpt+1)+"\n")
        file_res.write("#importation des packages propres ‡ python, extÈrieurs \n#et dÈveloppÈs spÈcialement pour cette application\n#‡ partir de packages externes :\n")
        file_res.write("from pymol import cmd\n\n")
        file_res.write("#affichage de la solution\n")
        file_res.write('#cmd.delete("all") #initialise l espace de travail dans pymol\n')
        file_res.write('cmd.load("'+filename+'") #charge le fichier spÈcifiÈ par filename\n')
        chaineRes="d"+string.split(os.path.basename(filename),".")[0]+"_sol_"+str(cpt+1)+"=('"+string.split(os.path.basename(filename),".")[0]+"',{"#initiation du tuple ; ce dernier ne comporte que deux parties : la premiere le nom de la structure, le second la solution
        #creation du tuple est rÈalisÈe plus bas, avant l'ecriture de chaineRes dans le fichier resultat
        file_res.write("removeAlt()\n")
        file_res.write('cmd.remove("hetatm")#permet de retirer les heteroatoms #attention : parfois certains fichiers sont mal annotÈs : un carbone alpha peu Ítre dÈcrit comme Ètant un hÈtÈro-atome\n')
        file_res.write('cmd.remove("h.")#permet de retirer toutes les molÈcules d eau\n')
        file_res.write('#cmd.hide("all")#permet de cacher l ensemble de la structure \n')
        file_res.write('cmd.hide("lines","'+string.split(os.path.basename(filename),".")[0]+'")#permet de cacher l ensemble de la structure courante \n')
        file_res.write('cmd.show("cartoon","'+string.split(os.path.basename(filename),".")[0]+'")#rÈaffiche la structure sous le format cartoon\n')
        file_res.write('#cmd.color("gray30","'+string.split(os.path.basename(filename),".")[0]+'")#colorise l ensemble de gris - pour faire ressortir ultÈrieurement par contraste les fragments sÈlectionnÈs\n')
        file_res.write('#mise en valeur de la sÈrine\n')
        t=utils.identificateur_pymol(sol[0])
        tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
        chaineRes+="'Ser' : '"+tmp+"',"
        file_res.write('cmd.show("spheres","'+tmp+'")\n')
        file_res.write('cmd.color("chocolate","'+tmp+'")\n')
        file_res.write('#mise en valeur du B5\n')
        chaineRes+="'B5' : ("+str(sol[1]['B5'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b5, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
        for e in sol[1]['B5'][1]:
            t=utils.identificateur_pymol(e)
            tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
            file_res.write('cmd.color("red","'+tmp+'")\n')
            chaineRes+="'"+tmp+"'" #ajout du residu actuel
            chaineRes+="," #fin pour ce residu, la virgule prepare pour le prochain residu
        chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
        chaineRes+="])"#fin de l'entree pour ce brin
        
        if sol[1]['B6']!=['none']:
            file_res.write('#mise en valeur du B6\n')
            chaineRes+=",'B6' : ("+str(sol[1]['B6'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b6, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
            for e in sol[1]['B6'][1]:
                t=utils.identificateur_pymol(e)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.color("orange","'+tmp+'")\n')
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="])"#fin de l'entree pour ce brin
                
        if sol[1]['B7']!=['none']:
            file_res.write('#mise en valeur du B7\n')
            chaineRes+=",'B7' : ("+str(sol[1]['B7'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b7, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
            for e in sol[1]['B7'][1]:
                t=utils.identificateur_pymol(e)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.color("tv_orange","'+tmp+'")\n') 
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="])"#fin de l'entree pour ce brin
               
        if sol[1]['B8']!=['none']:
            file_res.write('#mise en valeur du B8\n')
            chaineRes+=",'B8' : ("+str(sol[1]['B8'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b8, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
            for e in sol[1]['B8'][1]:
                t=utils.identificateur_pymol(e)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.color("brightorange","'+tmp+'")\n')         
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="])"#fin de l'entree pour ce brin
            
        if sol[1]['B4']!=['none']:
            file_res.write('#mise en valeur du B4\n')
            chaineRes+=",'B4' : ("+str(sol[1]['B4'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b4, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
            for e in sol[1]['B4'][1]:
                t=utils.identificateur_pymol(e)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.color("hotpink","'+tmp+'")\n')
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="])"#fin de l'entree pour ce brin
            
        if sol[1]['B3']!=['none']:
            file_res.write('#mise en valeur du B3\n')
            chaineRes+=",'B3' : ("+str(sol[1]['B3'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b3, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
            for e in sol[1]['B3'][1]:
                t=utils.identificateur_pymol(e)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.color("magenta","'+tmp+'")\n')
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="])"#fin de l'entree pour ce brin    
        
        if sol[1]['B2']!=['none']:
            file_res.write('#mise en valeur du B2\n')
            chaineRes+=",'B2' : ("+str(sol[1]['B2'][0])+",["#mise en place du dÈbut de l'entrÈe pour le brin b3, de l'identifiant et de l'initialisation de la liste des rÈsidus sÈlectionnÈs
            for e in sol[1]['B2'][1]:
                t=utils.identificateur_pymol(e)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.color("lightpink","'+tmp+'")\n')
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="])"#fin de l'entree pour ce brin 
            
        if sol[2]!=['none']:
            file_res.write('#mise en valeur des histidines\n')
            chaineRes+=",'His' : ["#mise en place du dÈbut de l'entrÈe pour l'histidine cette entrÈe contient les diffÈrents aa repondants aux criteres de selection
            for his in sol[2]:
                t=utils.identificateur_pymol(his)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.show("spheres","'+tmp+'")\n')
                file_res.write('cmd.color("yellow","'+tmp+'")\n')
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="]"#fin de l'entree l'histidine
                
        if sol[3]!=['none']:
            file_res.write('#mise en valeur des aa acides\n')
            chaineRes+=",'Acide' : ["#mise en place du dÈbut de l'entrÈe pour l'aa acide cette entrÈe contient les diffÈrents aa repondants aux criteres de selection
            for aa in sol[3]:
                t=utils.identificateur_pymol(aa)
                tmp=string.split(os.path.basename(filename),".")[0]+"&"+t
                file_res.write('cmd.show("spheres","'+tmp+'")\n')
                file_res.write('cmd.color("wheat","'+tmp+'")\n')
                chaineRes+="'"+tmp+"'" #ajout du residu actuel
                chaineRes+=","#fin pour ce residu, la virgule prepare pour le prochain residu
            chaineRes=chaineRes[:-1] #le dernier residu etant le dernier, il faut enlever la derniere virgule
            chaineRes+="]"#fin de l'entree l'aa acide
        
        chaineRes+='})'#termine le dictionnaire et le tuple        
        file_res.write('#les lignes suivantes permettent de mÈmoriser dans une entrÈe de dictionnaire, les diffÈrentes informartions prÈcÈdentes\n')
        file_res.write("d"+string.split(os.path.basename(filename),".")[0]+"_sol_"+str(cpt+1)+"=tuple()\n") #creation du tuple
        file_res.write(chaineRes)        
                
        
        file_res.close()
    return dir_res_sol

def rapport(filename,solution,avert,repertoire):
    """
        permet d'Ècrire un rapport de maniËre automatique ‡ partir des diffÈrentes solutions 
        et permet aussi de crÈer un script par solution proposÈe afin de les visualiser 
    """
    
    #initialisation pour la date:
    jour=datetime.date.today()
    
    #creation d'un nouveau fichier rapport pour le fichier dont le nom est specifie dans filename
    #si ce fichier existait deja, il sera recree
    #creation du nom de sortie peut etre a revoir
    file_res=open(repertoire+"/"+string.split(os.path.basename(filename),".")[0]+"_rapport.txt","w")
    #creation de l'entete du fichier
    file_res.write("##########################################################################################\n")
    file_res.write("# Rapport gÈnÈrÈ de maniËre automatique ‡ partir de l'application RESATOLC                \n")
    file_res.write("# application dÈveloppÈe lors du stage de fin d'Ètudes de la formation M2 BioInfo(Nantes) \n")
    file_res.write("# stage effectuÈe auprËs du Professeur Vinh TRAN au sein de l'U3B, unitÈ 6204 du CNRS     \n")
    file_res.write("# auteur: M. BARZINE (barzine@gmail.com)                                                  \n")
    file_res.write("# crÈdit: M. BARZINE - V. TRAN                                                            \n")
    file_res.write("##########################################################################################\n")
    file_res.write("# RESATOLC : rapport unitaire de structure                                                \n")    
    file_res.write("# Date de crÈation (AAAA/MM/JJ): "+jour.strftime("%Y/%m/%d")+"                                               \n")
    file_res.write("##########################################################################################\n")
    file_res.write("##########################################################################################\n")
    file_res.write("##########################################################################################\n")
    file_res.write("### Fichier resultat pour :         ######################################################\n")
    file_res.write("PDB_id="+string.split(os.path.basename(filename),".")[0]+"                                \n")
    file_res.write("# algorithme utilisÈ pour le calcul de l'assignation des structures secondaires :"+avert['algoSS']+"\n")
    file_res.write("##########################################################################################\n")    
    file_res.write("#### Avertissements divers :         #####################################################\n")
    file_res.write("# Avertissement(s) Èventuel(s) sur le fichier d'entrÈe :                                  \n")
    if avert['extract']=="":
        file_res.write("# Aucun avertissement\n")
    file_res.write("#"+avert['extract']+"\n")
    file_res.write("# Avertissement(s) Èventuel(s) gÈnÈrÈs lors du traitement :                               \n")    
    ### cette entree contient lui-mÍme un dictionnaire
    #file_res.write("# "+str(avert['traitement'])+"\n")
    if ("Ser" in avert['traitement']):
        file_res.write("#"+avert['traitement']['Ser']+"\n")
        file_res.write("#                                                                                        #\n")
        file_res.write("# Un traitement spÈcifique (et manuelle) est grandement conseillÈ pour ce cas            #\n")
        file_res.write("##########################################################################################\n")
        return
    file_res.write("# Bilan des solutions :                                                                   \n")
    if not ("complete" in avert['estimation']):
        file_res.write("#    Aucune solution complËte n'a pu etre trouvÈe\n")
    else :
        file_res.write("#    Solutions complËtes : "+str(avert['estimation']['complete'])+"          \n")
    if ("incomplete" in avert['estimation']):   
        file_res.write("#    Solutions incomplËtes : "+str(avert['estimation']['incomplete'])+"          \n")
    file_res.write("##########################################################################################\n")
    file_res.write("### DÈtails des solutions :        #######################################################\n")
    for cpt,sol in enumerate(solution) :
        #presentation de la solution
        file_res.write("###\n")
        #ID solution (permettra au biologiste de specifier la sol.qu'il aura choisi
        file_res.write("# ID solution : "+str(cpt+1)+"\n")
        file_res.write("#\n")
        file_res.write("# sÈrine de rÈfÈrence  : Ser "+sol[0][1]+" "+sol[0][0]+"                               \n")
        file_res.write("#\n")
        #file_res.write("#c-a-d que cette serine est considÈrÈe comme faisant partie du site catalytique pour  \n")
        #file_res.write("#la solution courante                                                                 \n") 
        #traitement des avertissements ou des remarques 
        file_res.write("# Avertissements ou remarques Èventuels sur cette solution :                          \n")
        file_res.write(formatageWarningSolution(sol[4]))
        file_res.write("#\n")
        file_res.write("#\n")
        #assignation des feuillets
        file_res.write("# Assignation des brins \n")
        if sol[1]['B2']!=['none']:
            file_res.write("#  B2 : "+str(sol[1]['B2'][0])+"         \n")
        else :
            file_res.write("#  B2 : n'a pas pu etre attribue\n")
        if sol[1]['B3']!=['none']:
            file_res.write("#  B3 : "+str(sol[1]['B3'][0])+"         \n")
        else :
            file_res.write("#  B3 : n'a pas pu etre attribue\n")
        if sol[1]['B4']!=['none']:
            file_res.write("#  B4 : "+str(sol[1]['B4'][0])+"         \n")
        else :
            file_res.write("#  B4 : n'a pas pu etre attribue\n")
        file_res.write("#  B5 : "+str(sol[1]['B5'][0])+"         \n")
        if sol[1]['B6']!=['none']:
            file_res.write("#  B6 : "+str(sol[1]['B6'][0])+"         \n")
        else :
            file_res.write("#  B6 : n'a pas pu etre attribue\n")
        if sol[1]['B7']!=['none']:
            file_res.write("#  B7 : "+str(sol[1]['B7'][0])+"         \n")
        else :
            file_res.write("#  B7 : n'a pas pu etre attribue\n")
        if sol[1]['B8']!=['none']:
            file_res.write("#  B8 : "+str(sol[1]['B8'][0])+"         \n")
        else :
            file_res.write("#  B8 : n'a pas pu etre attribue\n")
        file_res.write("#\n")
        file_res.write("#\n")
        #liste des Histidines Ètant ‡ une distance infÈrieur ‡ X angstrom
        file_res.write("# Liste des acides aminÈs candidats ‡ la participation du site catalytique :\n")
        if sol[2]!=['none']:
            file_res.write("#  Pour l'histidine : \n")
            for his in sol[2]:
                file_res.write("#    His "+str(his[1])+" "+str(his[0])+"\n")
        else :
            file_res.write("#  Pour l'histidine : aucune n'a ÈtÈ trouvÈe selon les critËres spÈcifiÈs\n")
        #liste des acides aminÈs acides ‡ une distance infÈrieur ‡ X angstrom
        if sol[3]!=['none']:
            file_res.write("#  Pour l'acide aminÈ acide : \n")
            for aa in sol[3]:
                file_res.write("#    "+str(aa[2])+" "+str(aa[1])+" "+str(aa[0])+"\n")
        else :
            file_res.write("#  Pour l'acide aminÈ acide : aucun aa n'a ÈtÈ trouvÈ selon les critËres spÈcifiÈs\n")
        file_res.write("###\n")
        file_res.write("###\n")
    
    file_res.close()
