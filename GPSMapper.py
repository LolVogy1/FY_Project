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
from collections import defaultdict

#TODO
#panning wihtout toolbar
#remove existing lines

"""only used if using A star algorithm"""
class Node:
    def __init__(self, parent = None, position = None):
        self.parent = parent
        self.position = position
        self.g = 0
        self.h = 0
        self.f = 0
    def __eq__(self, other):
        return self.position == other.position

"""graph class for use in Dijkstra's algorithm"""
class Graph:
    def __init__(self):
        #will be a dictionary of all connected nodes e.g {(X,Y):[(X:Y),(X,Y)]}
        self.edges = defaultdict(list)
        #a set of weights between every set of two connected nodes. 
        #The tuple of nodes is the key e.g {((X,Y),(X.Y)): 1.2,((X,Y),(X.Y)): 1.4 }
        self.weights = {}
        
    def add_edge(self, from_node, to_node, weight):
        #adds the edge and the weight of the edge
        #assumes edge is bi-directional
        self.edges[from_node].append(to_node)
        self.edges[to_node].append(from_node)
        self.weights[(from_node, to_node)] = weight
        self.weights[(to_node,from_node)] = weight
        
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
        #add button to mark clicks
        self.buttonMark = tk.Button(window, text = "Mark Start", command = lambda : self.toggleOn())
        self.buttonMark2 = tk.Button(window, text = "Mark End", command = lambda : self.toggleOn2())
        self.labelClick = tk.Label(window, text = "Last Click:")
        self.labelX = tk.Label(window)
        self.labelY = tk.Label(window)

        
    def toggleOn(self):
            markOn = self.markOn.get()
            if markOn == 0:
                #so you cant mark the start and the end point at the same time
                self.markOn.set(1)
                self.markOn2.set(0)
            else:
                self.markOn.set(0)
            
    def toggleOn2(self):
        markOn = self.markOn2.get()
        if markOn == 0:
            self.markOn.set(0)
            self.markOn2.set(1)
        else:
            self.markOn2.set(0)
    
    def openfile(self):
        #open file dialog window
        filepath = filedialog.askopenfilename(filetypes=(("Map files", "*.kml"),("All Files", "*.")))
        #get only the file name
        #filepath = os.path.split(filepath)[1]
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
        #attempt for different format
        if not coordList:
            linestrings = tree.findall('.//{http://www.opengis.net/kml/2.2}LineString')
            for attributes in linestrings:
                for subAttribute in attributes:
                    if subAttribute.tag == '{http://www.opengis.net/kml/2.2}coordinates':
                        x = subAttribute.text
                        #coordinates are stored in one long string
                        xlist = x.split("\n")
                        #first and last elements are '' so need to be removed
                        xlist.pop(0)
                        xlist.pop()
                        for item in xlist:
                            coord = list(map(float,item.split(",")))
                            coordList.append(coord)
            
        #draw the map
        self.drawmap(coordList)
                
    def drawmap(self,coordList):
     
        #make a 2D list of coordinates
        coordList2D = list()
        for i in coordList:
            coordList2D.append((i[0],i[1]))
            
        """converts a list of coordinates into a graph for pathfinding """    
        def graphConvert(coords):
            #create an empty list for the graph
            graph = []
            for i in range(0, len(coords)-1):
                #if the start of the edge is not the end node
                if i != len(coords)-1:
                    #an edge is a tuple of the coordinates of the start node, end node and the distance between them
                    dist =  geopy.distance.distance(coords[i],coords[i+1]).km
                    edge = (coords[i],coords[i+1], dist)
                    graph.append(edge)
                    #add edge in opposite direction (graph is undirected)
                    #dist =  geopy.distance.distance(coords[i+1],coords[i]).km
                    #edge = (coords[i+1],coords[i], dist)
                    #graph.append(edge)
                else:
                    break
            #make a graph object
            dGraph = Graph()
            #add the edges to the graph
            for edge in graph:
                dGraph.add_edge(*edge)
            return dGraph
            
        mGraph = graphConvert(coordList2D)
        
        self.buttonCalc = tk.Button(window, text = "Find Path", command = lambda : findPath(coordList2D, mGraph))
        self.buttonClear = tk.Button(window, text = "Clear Markers", command = lambda : clearWindow(coordList2D))
        self.buttonCalc.place(x=800,y=633)
        self.buttonClear.place(x=900,y=633)
        self.buttonMark.place(x=800,y=600)
        self.buttonMark2.place(x=900,y=600)
        
        #split coordinates into x and y
        long = [x[0] for x in coordList]
        lat = [y[1] for y in coordList]
        
        
        
        #creates a ckdtree for identifying closest data point to a click
        points = np.column_stack([long,lat])
        #has to be global
        global ckdtree
        ckdtree = scipy.spatial.cKDTree(points)
         
        #plot coordinates
        ax1.clear()
        ax1.plot(long,lat)
        #rescale axes
        ax1.relim()
        ax1.autoscale_view()
        ax1.axis('off')
        fig1.canvas.draw()
        scale = 1.5
        zp = ZoomPan()
        fzoom = zp.zoom_factory(ax1,base_scale = scale)
        fpan = zp.pan_factory(ax1)
        
        """display info on click"""
        def on_click(event):
            if event.inaxes is not None:
                #get the closest co-ordinates of the click
                nearestpoint = closest_point_coords(ckdtree, event.xdata, event.ydata).tolist()
                self.labelX.configure(text = nearestpoint[0])
                self.labelY.configure(text = nearestpoint[1])
                self.labelClick.place(x = 800, y = 500)
                self.labelX.place(x = 800, y = 530)
                self.labelY.place(x = 800, y = 560)
                #plots a circle on click (Start)
                if self.markOn.get() == 1:
                    #mark on map
                    start, = ax1.plot(nearestpoint[0],nearestpoint[1], 'o')
                    fig1.canvas.draw()
                    #record coordinates
                    self.startx.set(nearestpoint[0])
                    self.starty.set(nearestpoint[1])
                    #disable the button and stops you from adding another point
                    self.markOn.set(0)
                    self.buttonMark['state']= 'disabled'
                #plots a square on click (End)
                elif self.markOn2.get() == 1:
                    ax1.plot(nearestpoint[0],nearestpoint[1], 's')
                    self.endx.set(nearestpoint[0])
                    self.endy.set(nearestpoint[1])
                    fig1.canvas.draw()
                    self.markOn2.set(0)
                    self.buttonMark2['state']= 'disabled'
                
            else:
                #if the click occurred outside the axes
                print("clicked outside")
            
        #connect on click event to the figure
        fig1.canvas.callbacks.connect('button_press_event',on_click)
        
        """clears the additional markers and path"""
        def clearWindow(coords):
            long = [x[0] for x in coords]
            lat = [y[1] for y in coords]
            #clear the plot
            ax1.clear()
            #redraw
            ax1.plot(long,lat)
            ax1.axis('off')
            #ax1.axis('off')
            fig1.canvas.draw()
            #restore button functionality
            self.buttonCalc['state'] = 'normal'
            self.buttonMark['state'] = 'normal'
            self.buttonMark2['state'] = 'normal'
        
        def findPath(coords, mapGraph):
            #if you haven't set a start and end point
            testStart = self.startx.get()
            testEnd = self.endx.get()
            startPoint = (testStart, self.starty.get())
            endPoint = (testEnd, self.endy.get())
            if testStart == 0.0 or testEnd == 0.0:
                mb.showerror("Missing Points!","You haven't marked a start/end point")
            else:
                #calculate the path
                path = dijkstra(startPoint, endPoint, mapGraph)
                long = [x[0] for x in path]
                lat = [y[1] for y in path]
                #plot
                ax1.plot(long,lat,color = "blue")
                fig1.canvas.draw()
                #disable button
                self.buttonCalc['state'] = 'disabled'

        def dijkstra(start, end, dGraph):
                        
            shortestPaths = {start: (None,0)}
            current = start
            visited = set()
            
            while current != end:
                visited.add(current)
                destinations = dGraph.edges[current]
                weightToCurrent = shortestPaths[current][1]
                
                for nextNode in destinations:
                    weight = dGraph.weights[(current,nextNode)] + weightToCurrent
                    if nextNode not in shortestPaths:
                        shortestPaths[nextNode] = (current, weight)
                    else:
                        currentShortestWeight = shortestPaths[nextNode][1]
                        if currentShortestWeight > weight:
                            shortestPaths[nextNode] = (current, weight)
                            
                nextDestinations = {node: shortestPaths[node] for node in shortestPaths 
                                    if node not in visited}
                if not nextDestinations:
                    return "Route not Possible"
                
                current = min(nextDestinations, key = lambda k: nextDestinations[k][1])
                
            path = []
            while current is not None:
                path.append(current)
                nextNode = shortestPaths[current][0]
                current = nextNode
            return path[::1]
                
                
        """A* algorithm"""
        """"not used"""
    
        def aStar(start,end, graph):
            
            openList = list()
            closedList = list()
            startNode = Node(None, start)
            startNode.g = startNode.h = startNode.f = 0
            endNode = Node(None, end)
            endNode.g = endNode.h = endNode.f = 0
            openList.append(startNode)
            
            while len(openList) > 0:
                currentNode = openList[0]
                currentIndex = 0
                
                for index, item in enumerate(openList):
                    if item.f < currentNode.f:
                        currentNode = item
                        currentIndex = index
                openList.pop(currentIndex)
                closedList.append(currentNode)
                print(currentNode.position)
                
                if currentNode == endNode:
                    print("reached the end")
                    path = []
                    current = currentNode
                    while current is not None:
                        path.append(current.position)
                        current = current.parent
                    return path[::-1]
                children = []
                distanceList = []
                
                for edge in graph:
                    if edge[0] == currentNode.position:
                        newNode = Node(currentNode, edge[1])
                        children.append(newNode)
                        distanceList.append(edge[2])
                        
                for index, child in enumerate(children):
                    #if the child has been visited already
                    #ignore it
                    for closedChild in closedList:
                        if child == closedChild:
                            continue
                    child.g = currentNode.g + distanceList[index]
                    child.h = geopy.distance.distance(child.position,endNode.position).km
                    child.f = child.g + child.h
                    for openNode in openList:
                        if child == openNode and child.g > openNode.g:
                            continue
                    openList.append(child)

        


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
'''

class ZoomPan:
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None


    def zoom_factory(self, ax, base_scale = 2.):
        def zoom(event):
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
        fig.canvas.mpl_connect('scroll_event',zoom)
        #return the function
        return zoom

    def pan_factory(self, ax):
        def onPress(event):
            if event.inaxes != ax: return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press

        def onRelease(event):
            self.press = None
            ax.figure.canvas.draw()

        def onMotion(event):
            if self.press is None: return
            if event.inaxes != ax: return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax.set_xlim(self.cur_xlim)
            ax.set_ylim(self.cur_ylim)

            ax.figure.canvas.draw()

        fig = ax.get_figure() # get the figure of interest

        # attach the call back
        fig.canvas.mpl_connect('button_press_event',onPress)
        fig.canvas.mpl_connect('button_release_event',onRelease)
        fig.canvas.mpl_connect('motion_notify_event',onMotion)

        #return the function
        return onMotion



#set up graph
fig1, ax1 = plt.subplots()
fig1.set_size_inches(8, 6)
ax1.axis('off')

#set up GUI
window = tk.Tk()
window.geometry("1000x700")
window.title("GPS Mapper")
#initialise GUI 
mapgui = MapGUI(window)
canvas = FigureCanvasTkAgg(fig1, window)   
canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
#toolbar = NavigationToolbar2Tk(canvas, window) 
#toolbar.update() 
window.mainloop()
    

