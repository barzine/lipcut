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

#libraries and packages importation 
import copy

#def feuillet_unique(liste):
def remove_duplicate(liste):
    """
        Remove all duplicates of a beta sheet record from the input list: pdb files are
        using several times the same beta sheet record to describe some specific tertionary structure
        e.g. beta barrel. 
        A beta sheet record is considered as duplicate if it has the same N and C terminal amino-acid
        i.e. these two amino-acid have the exact same identifier for the two beta sheet record.
        
        This process does not take the orientation of the beta sheet in consideration
    """

    #Check if there is a need of looking for duplicates
    if len(liste)<2 : #if there are less than 2 elements, there can not be any duplicate
        #the list is return without any modification
        return liste
  
    #if 2 or more elements in the list:  
    nv_list=[liste.pop(0)] #initialisation of the list which will contain only unique beta sheet 
    #and that would later be returned
    #reminder: pop will remove the current element from the original list
   
    #each remaining element of the list is checked 
    #that it is not a duplicate of one of the record already in the return list
    for e in liste :#element from the original list
        #the beta sheet is added in the return list only if its first and last aa-s are not already in the return list
        found=False
        i=0
        while not(found) and i<len(nv_list): #as long as no duplicates has been found
            if ((e[1],e[2],e[3])==(nv_list[i][1],nv_list[i][2],nv_list[i][3]) and (e[4],e[5],e[6])==(nv_list[i][4],nv_list[i][5],nv_list[i][6])): #or ((e[1],e[2],e[3])==(i[4],i[5],i[6]) and (e[4],e[5],e[6])==(i[1],i[2],i[3])) :
                found=True
            i+=1
        if not found : #the beta sheet is unique, hence added to return list
            nv_list.append(e)
  
    return (nv_list)




#####
#def verif_integrite(listBrin,index,listCA):
def integrity_check(listBeta,index,listCA):

    """
        This function serves sevaral purposes: 
        1/ Check that each beta sheet of a given list has all its constitutive atoms defined 
        2/ Replace the chain identifier by the corresponding key in the alpha Carbon index 
        and add a 9th field for recording the beta sheet integrity
    """
    #initialisation 
    warning_beta="" #record the warnings 
    listR=list() #
    t=tuple()
    cles=index.keys()
    
    for b in listBeta:
        #for each beta sheet, every alpha Carbon, from the N-terminal to the C-terminal or to the first missing alpha carbon
        start=False
        end=False
        IDstart=""
        found=False
        k=0#counter initialisation for the key
        while (not(found)): #stop as soon as the correct entry of the dictiaonnary is found 
            if ((not start)and b[2]==index[cles[k]][4]):
                #chain identifier is the same as the one of the current dictioannary entry
                #check if this is the entry we are looking for                
                if int(index[cles[k]][0])<=int(b[3])and int(index[cles[k]][2])>=int(b[3]):
                    #the number of the inital residue is comprised in the current record
                    found=True # we have found the correct beta sheet entry with the N-term we were looking for
                    start=True#True now that we have the start
                    IDstart=copy.deepcopy(cles[k])#Referencing the key for the fragment for N terminal extremity of the beta-sheet
                    if index[cles[k]][4]==b[5] and int(index[cles[k]][0])<=int(b[6]) and int(index[cles[k]][2])>=int(b[6]):
                        end =True
                    else:
                        end=False
            k+=1 
        
        #creation of the record with the new description
        if end:#the fragment that contains the  N-terminal part is also the one containing the c-terminal part
            t=(b[0],b[1],IDstart,b[3],b[4],IDstart,b[6],b[7],1)
        else :#the fragment containing the N-terminal part does not contain the C-terminal, i.e.
            #some residues are not described in the original file
            t=(b[0],b[1],IDstart,b[3],b[4],b[5],b[6],b[7],0)
            #warning is then created for the report:
            warning_beta+="\nThe beta sheet '"+b[0]+"' is not correctly described in the original file"
        listR.append(t)
    return [listR,warning_beta]
