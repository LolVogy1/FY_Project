# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 16:14:01 2021

@author: Alex Vong
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

#Mathematical modules
import numpy as np
import scipy.spatial
import geopy.distance
from collections import defaultdict


"""graph class for use in Dijkstra's algorithm"""
class Graph:
    def __init__(self):
        #edges are a dictionary of all possible next nodes
        #e.g {A: [B, D, E], B:[...]}
        self.edges = defaultdict(list)
        #weights contains the weights between each pair of nodes.
        #The pair of nodes forms a tuple which acts as the key
        #e.g {(X,Y):12}
        self.weights = {}
        
    def add_edge(self, from_node, to_node, weight):
        #adds the edge and the weight of the edge
        #assumes edge is bi-directional
        self.edges[from_node].append(to_node)
        self.edges[to_node].append(from_node)
        self.weights[(from_node, to_node)] = weight
        self.weights[(to_node,from_node)] = weight
        
"""enables use of panning and zooming in matplotlib"""
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

    """zooming"""
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
    
    """panning"""
    def pan_factory(self, ax):
        #get data of mouse press
        def onPress(event):
            
            if event.inaxes != ax: return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press
        
        #re-draw the figure when click released
        def onRelease(event):
            
            self.press = None
            ax.figure.canvas.draw()
            
        #adjusts the limits of the axes to simulate movement
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

"""GUI class"""   
class MapGUI:
    def __init__(self, window):
        self.window = window
        #button to open a file
        self.buttonFile = tk.Button(window, text = "Open File", command = lambda : self.openfile())
        self.buttonFile.pack()
        #variable to track whether marking is enabled
        self.markStart = tk.IntVar()
        self.markStart.set(0)
        self.markEnd = tk.IntVar()
        self.markEnd.set(0)
        #tracks co-ordinates of start and end points of shortest path
        #needs to be globally acessable
        self.startx = tk.DoubleVar()
        self.starty = tk.DoubleVar()
        self.endx = tk.DoubleVar()
        self.endy = tk.DoubleVar()
        #add button to mark clicks
        self.buttonMarkStart = tk.Button(window, text = "Mark Start", command = lambda : self.toggleStart())
        self.buttonMarkEnd = tk.Button(window, text = "Mark End", command = lambda : self.toggleEnd())
        #labels to display co-ordinates of last click
        self.labelClick = tk.Label(window, text = "Last Click:")
        self.labelX = tk.Label(window)
        self.labelY = tk.Label(window)

    #toggles marking the start points
    def toggleStart(self):
        self.markStart.set(1)
        self.markEnd.set(0)
    
    #toggles marking the end point
    def toggleEnd(self):
        self.markStart.set(0)
        self.markEnd.set(1)

    
    #open a file
    def openfile(self):
        #open file dialog window
        filepath = filedialog.askopenfilename(filetypes=(("Map files", "*.kml"),("All Files", "*.")))
        filetype = filepath.split(".")[1]
        if filetype == "kml":
            #parse the file
            self.openmap(filepath)
        else:
            mb.showerror("Invalid file", "File not KML")
    
    #gets data from the file
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
        
    """Draw the map and enable interactive functions"""
    def drawmap(self,coordList):
        #make a 2D list of coordinates (data is 3D)
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
                else:
                    break
            #make a graph object
            dGraph = Graph()
            #add the edges to the graph
            for edge in graph:
                dGraph.add_edge(*edge)
            return dGraph
        
        #create a graph from the co-ordinates    
        mGraph = graphConvert(coordList2D)
        print(mGraph.edges)
        
        #buttons to calculate path and mark points appear once the map is drawn
        self.buttonCalc = tk.Button(window, text = "Find Path", command = lambda : findPath(coordList2D, mGraph))
        self.buttonClear = tk.Button(window, text = "Clear Markers", command = lambda : clearWindow(coordList2D))
        self.buttonCalc.place(x=800,y=633)
        self.buttonClear.place(x=900,y=633)
        self.buttonMarkStart.place(x=800,y=600)
        self.buttonMarkEnd.place(x=900,y=600)
        #make sure start/end points are unmarked
        self.startx.set(0.0)
        self.starty.set(0.0)
        self.endx.set(0.0)
        self.endy.set(0.0)
        
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
        #enable panning and zooming
        scale = 1.5
        zp = ZoomPan()
        fzoom = zp.zoom_factory(ax1,base_scale = scale)
        fpan = zp.pan_factory(ax1)
        
        """display info on click"""
        def on_click(event):
            if event.inaxes is not None:
                #get the closest co-ordinates of the click
                nearestpoint = closest_point_coords(ckdtree, event.xdata, event.ydata).tolist()
                #labels displaying data appear on click
                self.labelX.configure(text = nearestpoint[0])
                self.labelY.configure(text = nearestpoint[1])
                self.labelClick.place(x = 800, y = 500)
                self.labelX.place(x = 800, y = 530)
                self.labelY.place(x = 800, y = 560)
                #plots a circle on click (Start)
                if self.markStart.get() == 1:
                    #mark on map
                    start, = ax1.plot(nearestpoint[0],nearestpoint[1], 'o')
                    fig1.canvas.draw()
                    #record coordinates
                    self.startx.set(nearestpoint[0])
                    self.starty.set(nearestpoint[1])
                    #disable the button and stops you from adding another point
                    self.markStart.set(0)
                    self.buttonMarkStart['state']= 'disabled'
                #plots a square on click (End)
                elif self.markEnd.get() == 1:
                    ax1.plot(nearestpoint[0],nearestpoint[1], 's')
                    self.endx.set(nearestpoint[0])
                    self.endy.set(nearestpoint[1])
                    fig1.canvas.draw()
                    self.markEnd.set(0)
                    self.buttonMarkEnd['state']= 'disabled'
                
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
            fig1.canvas.draw()
            #restore button functionality
            self.buttonCalc['state'] = 'normal'
            self.buttonMarkStart['state'] = 'normal'
            self.buttonMarkEnd['state'] = 'normal'
            #start/end points are no longer marked
            self.startx.set(0.0)
            self.starty.set(0.0)
            self.endx.set(0.0)
            self.endy.set(0.0)
        
        """calculate the shortest path and plot it"""
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

        """Dijkstra's algorithm for finding the shortest path"""
        def dijkstra(start, end, dGraph):
            #shortest paths is a dict of nodes
            #whose value is a tuple of (previous node, weight)
            shortestPaths = {start: (None,0)}
            current = start
            visited = set()
            #while we haven't reached the end
            while current != end:
                #add the current node to the set of visited nodes
                visited.add(current)
                #get all of the neighbouring nodes
                destinations = dGraph.edges[current]
                #get the weight to the current node
                weightToCurrent = shortestPaths[current][1]
                
                #for each neighbour
                for nextNode in destinations:
                    #weight is distance from current to next + weight from start to current
                    weight = dGraph.weights[(current,nextNode)] + weightToCurrent
                    #if we haven't explored this path
                    if nextNode not in shortestPaths:
                        #add the path
                        shortestPaths[nextNode] = (current, weight)
                    else:
                        currentShortestWeight = shortestPaths[nextNode][1]
                        #if there is a shorter path
                        if currentShortestWeight > weight:
                            shortestPaths[nextNode] = (current, weight)
                
                #next node is the destination with the lowest weight
                nextDestinations = {node: shortestPaths[node] for node in shortestPaths 
                                    if node not in visited}
                if not nextDestinations:
                    return "Route not Possible"
                
                current = min(nextDestinations, key = lambda k: nextDestinations[k][1])
            
            #Work back through destinations in shortest path    
            path = []
            while current is not None:
                path.append(current)
                nextNode = shortestPaths[current][0]
                current = nextNode
            #return the reversed path (start to end)
            return path[::1]
                

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
    


