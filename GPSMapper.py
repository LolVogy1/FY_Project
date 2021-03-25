# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:14:01 2021

@author: Alex
"""
#Data Processing
from xml.etree import ElementTree as et
import os

#GUI
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox as mb

#Matplotlib modules
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  NavigationToolbar2Tk)
from matplotlib import style
style.use('ggplot')
import numpy as np
import scipy.spatial
import geopy.distance

#TODO
#panning wihtout toolbar
#add plot on top of current one
#A* only works in a certain direction (make graph undirected)

class MapGUI:
    
    def __init__(self, window):
        self.window = window
        #button to open a file
        self.buttonFile = tk.Button(window, text = "Open File", command = lambda : self.openfile())
        self.buttonFile.pack()
        #variable to track whether marking is enabled
        self.markOn = tk.IntVar()
        self.markOn.set(0)
        self.markOn2 = tk.IntVar()
        self.markOn2.set(0)
        #tracks co-ordinates of start and end points of shortest path
        #there is no variable for tuples within tkinter 
        self.startx = tk.DoubleVar()
        self.starty = tk.DoubleVar()
        self.endx = tk.DoubleVar()
        self.endy = tk.DoubleVar()
    
    def openfile(self):
        #open file dialog window
        filepath = filedialog.askopenfilename(filetypes=(("Map files", "*.kml"),("All Files", "*.")))
        #get only the file name
        filepath = os.path.split(filepath)[1]
        #parse the file
        self.openmap(filepath)
            
    def openmap(self,filename):
        #set up an empty list
        coordList = list()
        #parses the XML data
        tree = et.parse(filename)
        #gets the first element in the XML/KML file
        root = tree.getroot()
        #use .iter to find every element in the sub tree where the tag matches below
        for child in root.iter('{http://www.google.com/kml/ext/2.2}coord'):
                #get the text from the element
                x = child.text
                #text is "1 2 3" etc
                #Coordinates stored as longitude, latitude, altitude
                #split the text into an array and convert from string to float
                #either works but map function is marginally quicker
                """coord = [float(y) for y in x.split()]"""
                coord = list(map(float, x.split()))
                #add to list
                coordList.append(coord)
        #draw the map
        self.drawmap(coordList)
                
    def drawmap(self,coordList):
        #set up graph
        fig = plt.figure(figsize = (8,6), dpi= 100)
        plot1 = fig.add_subplot(111,picker =1)
        
        #make a 2D list of coordinates
        coordList2D = list()
        for i in coordList:
            coordList2D.append((i[0],i[1]))
        
        #add button to mark clicks
        self.buttonMark = tk.Button(window, text = "Mark Start", command = lambda : toggleOn())
        self.buttonMark.pack()
        self.buttonMark2 = tk.Button(window, text = "Mark End", command = lambda : toggleOn2())
        self.buttonMark2.pack()
        self.buttonCalc = tk.Button(window, text = "Find Path", command = lambda : findPath(coordList2D, fig))
        self.buttonCalc.pack()
        
        #split coordinates into x and y
        long = [x[0] for x in coordList]
        lat = [y[1] for y in coordList]
        
        #creates a ckdtree for identifying closest data point to a click
        points = np.column_stack([long,lat])
        ckdtree = scipy.spatial.cKDTree(points)
         
        canvas = FigureCanvasTkAgg(fig, window)   
        canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
        
        # creating the Matplotlib toolbar 
        toolbar = NavigationToolbar2Tk(canvas, window) 
        toolbar.update() 
        
        #plot coordinates
        plot1.plot(long,lat)
        plot1.axis('off')
        #scale = 1.5
        #f = zoom_factory(plot1,base_scale = scale)
        
        def toggleOn():
            markOn = self.markOn.get()
            if markOn == 0:
                #so you cant mark the start and the end point at the same time
                self.markOn.set(1)
                self.markOn2.set(0)
            else:
                self.markOn.set(0)
            
        def toggleOn2():
            markOn = self.markOn2.get()
            if markOn == 0:
                self.markOn.set(0)
                self.markOn2.set(1)
            else:
                self.markOn2.set(0)
                
        def findPath(coords, fig):
            #if you haven't set a start and end point
            testStart = self.startx.get()
            testEnd = self.endx.get()
            startPoint = (testStart, self.starty.get())
            endPoint = (testEnd, self.endy.get())
            if testStart == 0.0 or testEnd == 0.0:
                mb.showerror("Missing Points!","You haven't marked a start/end point")
            else:
                plot2 = fig.add_subplot(111)
                mapGraph = graphConvert(coords)
                path = aStar(startPoint, endPoint, mapGraph)
                long = [x[0] for x in path]
                lat = [y[0] for y in path]
                plot2.plot(long, lat)
                fig.canvas.draw()
                print(path)

                
                
        """converts a list of coordinates into a graph for pathfinding """    
        def graphConvert(coords):
            #create an empty list for the graph
            graph = list()
            for i in range(0, len(coords)-1):
                #if the start of the edge is not the end node
                if i != len(coords)-1:
                    #an edge is a tuple of the coordinates of the start node, end node and the distance between them
                    dist =  geopy.distance.distance(coords[i],coords[i+1]).km
                    edge = (coords[i],coords[i+1], dist)
                    graph.append(edge)
                    #add edge in opposite direction (graph is undirected)
                    #edge = (coords[i+1],coords[i], dist)
                    #graph.append(edge)
                else:
                    break
            return graph
        
        """A* algorithm"""
        def aStar(start, endp, graph):
            #empty list of visited nodes
            visited = list()
            #add the starting point to the list of visited nodes
            visited.append(start)
            current = start
            end = endp
            print("Start:",current)
            print("End:",end)
            while current != end:
                #empty list of neighbours
                neighbours = list()
                f = float('inf')
                
                """getting all of the neighbouring nodes"""
                for edge in graph:
                    #find edges from the current node
                    if edge[0] == current:
                        #estimate the distance from the neighbour to the end
                        h = geopy.distance.distance(edge[0], end).km
                        #the neighbour is a tuple of the  next point, distance to the point and the estimated distance to the end
                        neighbour = (edge[1], edge[2], h)
                        neighbours.append(neighbour)
                        
                """find the best neighbour"""
                for node in neighbours:
                    #if the end node is found
                    if node[0] == end:
                        visited.append(node[0])
                        current = end
                        print("End reached")
                        break;
                    else:
                        #distance of edge + distance to end
                        tempF = node[1] + node[2]
                        #if the total distance is lower than the current neighbour
                        if tempF < f:
                            f = tempF
                            #select this node as the best one
                            current = node[0]
                #add the neighbour to the visited nodes
                visited.append(current)
                #set the current node to this one
                print("next node:",current)
            return visited
         
        
        """display info on click"""
        def on_click(event):
            if event.inaxes is not None:
                #get the closest co-ordinates of the click
                nearestpoint = closest_point_coords(ckdtree, event.xdata, event.ydata).tolist()
                #plots a circle on click (Start)
                if self.markOn.get() == 1:
                    #mark on map
                    plot1.plot(nearestpoint[0],nearestpoint[1], 'o')
                    fig.canvas.draw()
                    #record coordinates
                    self.startx.set(nearestpoint[0])
                    self.starty.set(nearestpoint[1])
                    #hide the button to mark the start point
                    self.markOn.set(0)
                    self.buttonMark.pack_forget()
                #plots a square on click (End)
                elif self.markOn2.get() == 1:
                    plot1.plot(nearestpoint[0],nearestpoint[1], 's')
                    self.endx.set(nearestpoint[0])
                    self.endy.set(nearestpoint[1])
                    fig.canvas.draw()
                    self.markOn2.set(0)
                    self.buttonMark2.pack_forget()
                
            else:
                #if the click occurred outside the axes
                print("clicked outside")
        
                
        #connect on click event to plot
        fig.canvas.callbacks.connect('button_press_event',on_click)
        
        


"""for getting the closest point to a click"""

def closest_point_distance(ckdtree, x, y):
    #returns distance to closest point
    return ckdtree.query([x, y])[0]

def closest_point_id(ckdtree, x, y):
    #returns index of closest point
    return ckdtree.query([x, y])[1]

def closest_point_coords(ckdtree, x, y):
    # returns coordinates of closest point
    return ckdtree.data[closest_point_id(ckdtree, x, y)]
    # ckdtree.data is the same as points
'''
#zooming with mouse wheel
def zoom_factory(ax,base_scale = 2.):
    def zoom_fun(event):
        # get the current x and y limits
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
        cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1/base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            print(event.button)
        # set new limits
        ax.set_xlim([xdata - cur_xrange*scale_factor,
                     xdata + cur_xrange*scale_factor])
        ax.set_ylim([ydata - cur_yrange*scale_factor,
                     ydata + cur_yrange*scale_factor])
        plt.draw() # force re-draw

    fig = ax.get_figure() # get the figure of interest
    # attach the call back
    fig.canvas.mpl_connect('scroll_event',zoom_fun)
    #return the function
    return zoom_fun
'''
    

#set up GUI
window = tk.Tk()
window.geometry("1000x700")
window.title("GPS Mapper")
#initialise GUI 
mapgui = MapGUI(window)
window.mainloop()

