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


source="C:/Documents and Settings/parissa/Mes documents/tests/d2LIP_stride/uploads/JIDU35909L.poc"
file=open(source,'r')
res=open("C:/Documents and Settings/parissa/Mes documents/tests/d2LIP_stride/uploads/rap.txt","w")
resn=0
num="a"
for line in file :
    print line
    if line[68:71]=="33 ":
        res.write(line)
        if line[24:26]!=num:
            num=line[24:26]
            resn+=1

res.write("nombre de residus :"+str(resn))
res.close()
file.close()
        


