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


def verif_site(liste_ser,beta,site_desc,atom_L,verbose=False):
    """
        Cross information between the list of serines that might have catalytic properties
        and the information recorded in the pdb file for the catalytic site description 
    """
    if (verbose):
        print "in verif_site function" #debug
    list_ser_filtre=[]#record information of serines described in the catalytic field of the file
    #the serines are (if) discribed in the field "SITE"
    site=site_desc.split("#") #retrieve the elements of the lignes as a list 
    for s in site: #for each elemnt in the site variable
        if (verbose):
            print s #debug
        for e in liste_ser : #for each serine which can be considered as catalytic 
            if (verbose):
                print atom_L[int(e[0][0])] #debug
            if atom_L[int(e[0][0])][0] in s :
                if e not in list_ser_filtre:
                    list_ser_filtre.append(e)
    if (verbose):
        print "ser filtred" #debug 
        print list_ser_filtre #debug                   
        print "end of verif_site function" #debug
    return list_ser_filtre




def verif_histidine(list_aa,d,verbose=False):
    """
        seek histidine which can be catalytic in an imput list of aa
    """
    possibilities=list()
    avertissement=""
    if (verbose):
        for aa in list_aa:
            if aa[2]=='HIS':
                preliste.append(aa)
    
    #Check that this histidine is located after the LAST beta sheet:
    if d["B8"]!=['none']: #if there is a description for the beta sheet refered as beta sheet 8
        #the histidine has to be located at the C-terminal of the beta sheet
        if (verbose):
            print aa[0] #debug
            print d["B8"][0]#debug
        for aa in list_aa:
            if aa[2]=='HIS':
                if int(d["B8"][0][6])<int(aa[0]):
                    possibilities.append(aa) 
           
    elif d["B7"]!=['none']: #if the last beta sheet is the one refered as beta sheet 7
        if verbose:
            print aa[0] #debug
            print d["B7"][0]#debug
        for aa in list_aa:
            if aa[2]=='HIS':
                if int(d["B7"][0][6])<int(aa[0]):
                    avertissement="The His is (are) located after the B7 since there isn't any B8 in the current solution"
                    possibilities.append(aa)
                    
    elif d["B6"]!=['none']: #same logic as in elif d["B7"]!=['none']
        if verbose:
            print aa[0]#debug
            print d["B6"][0]#debug
        for aa in list_aa:
            if aa[2]=='HIS':
                if int(d["B6"][0][6])<int(aa[0]):
                    avertissement="The His is (are) located after the B6 since there isn't any B8 or B7 in the current solution"
                    possibilities.append(aa)                                       
    
    if len(possibilities)>1:
        avertissement+="-Remark: there are more than one possible solution"
                    
    if len(possibilities)>0:
        if verbose:
            print len(possibilities)
            print possibilities
        #There is at least one HIS with all the needed criteria
        for his in possibilities:
            utils.color_struct(his,"HIS")    
        return [possibilities,avertissement]
    else :
        #No HIS has been found with the correct criteria to be considered as catalytic
        return [['none'],"The current solution does't have any valid catalytic HIS"]
    
    
def verif_acide(list_aa,d, verbose=False):
    """
        seek for an acid aa which can be catalytic in an imput list of aa
    """
    possibilities=list()
    avertissement=""
    #check if the acidic aa is after the antepenultimate beta sheet:
    if d["B8"]!=['none']:
        for aa in list_aa:
            if aa[2] in ['ASP','GLU'] :
                #B8 and B7 have BOTH a a description:
                #for avoiding problems of multiples residus:
                if aa[0].isdigit():
                    if int(d["B7"][0][6])<int(aa[0]):
                        possibilities.append(aa)
                else:
                    bb=(aa[0][0:-1],aa[1],aa[2])
                    if int(d["B7"][0][6])<int(bb[0]):
                        possibilities.append(bb)
    elif d["B7"]!=['none']:
        for aa in list_aa:
            if aa[2] in ['ASP','GLU'] :
                #there is a description for the beta 7
                #HIS should be at the C-term of this beta 7 since there isn't any definition for Beta 8
                if verbose:
                    print aa[0] #debug
                    print d["B7"][0]#debug
                avertissement="The acidic aa is (are) located after B6, since there isn't any B8 and HIS has to be after B7"
                if aa[0].isdigit():
                    if int(d["B6"][0][6])<int(aa[0]):
                        possibilities.append(aa)
                else:
                    bb=(aa[0][0:-1],aa[1],aa[2])
                    if int(d["B6"][0][6])<int(bb[0]):
                        possibilities.append(bb)
    elif d["B6"]!=['none']:
        for aa in list_aa:
            if aa[2] in ['ASP','GLU'] :        
                #same logic as before:
                if verbose:
                    print aa[0]#debug
                    print d["B6"][0]#debug
                avertissement="The acidic aa is (are) located after B5, since there isn't any B7 and HIS has to be after B6 in the current solution"                if aa[0].isdigit():
                    if int(d["B5"][0][6])<int(aa[0]):
                        possibilities.append(aa)
                else:
                    bb=(aa[0][0:-1],aa[1],aa[2])
                    if int(d["B5"][0][6])<int(bb[0]):
                        possibilities.append(bb)
    
    if len(possibilities)>1:
        avertissement+="- Remark : there are more than one solution"
    
    if len(possibilities)>0:
        if verbose:
            print possibilities
        #There is at least one acidic aa which corresponds to the needed criteria to be considered as possibly catalytic
        for acide in possibilities:
            utils.color_struct(acide,"ACIDIC")    
        return [possibilities,avertissement]
    else :
        #There wasn't acidic aa found which can be catalytic
        return [['none'],"The current solution doesn't have any valid catalytic acidic aa"]
    


