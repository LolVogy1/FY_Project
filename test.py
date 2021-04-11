# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 17:29:33 2021

@author: Alex
"""
from pykml import parser 

root = parser.fromstring(open('map.kml', 'r').read())
print(root.Document.Placemark.LineString.coordinates)