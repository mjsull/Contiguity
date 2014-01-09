#!/usr/bin/env python
# Contiguity   Written by: Mitchell Sullivan   mjsull@gmail.com
# Supervisor: Dr. Scott Beatson
# Version 0.9 08.01.2014
# License: GPLv3

from Tkinter import *
import math
import tkFileDialog
import tkSimpleDialog
import tkFont
import tkMessageBox
import gzip
import os
import subprocess
import string
import threading
import Queue

transtab = string.maketrans('atcgATCG', 'tagcTAGC')


class contig:
    def __init__(self, name, sequence):
        self.name = name
        self.forseq = sequence.lower()
        tempseq = self.forseq[::-1]
        self.revseq = tempseq.translate(transtab)
        self.length = len(sequence)
        self.xlength = None
        if self.length >= 1000000000:
            self.strlen = str(self.length / 1000000000) + 'Gb'
        elif self.length >= 1000000:
            self.strlen = str(self.length / 1000000) + 'Mb'
        elif self.length >= 1000:
            self.strlen = str(self.length / 1000) + 'Kb'
        else:
            self.strlen = str(self.length) + 'bp'
        self.visible = False
        self.to = []
        self.fr = []
        self.xpos = None
        self.ypos = None
        self.orient = [True]
        self.short = False
        self.dupnumber = 1
        self.startanchor = 0
        self.stopanchor = 0



class App:
    def __init__(self, master):
        self.menubar = Menu(master)
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Create CSAG file", command=self.create_edges)
        self.filemenu.add_command(label="Create comparison", command=self.create_comp)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Load assembly", command=self.load_assembly)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=root.quit)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.viewmenu = Menu(self.menubar, tearoff=0)
        self.viewmenu.add_command(label="View graph", command=self.view_options)
        self.viewmenu.add_command(label="Self comparison", command=self.self_compare)
        self.viewmenu.add_command(label="Add contig", command=self.add_contig_dialogue)
        self.viewmenu.add_separator()
        self.viewmenu.add_command(label="Find contig", command=self.find_contig)
        self.viewmenu.add_separator()
        self.viewmenu.add_command(label="Zoom in", command=self.zoominmenu)
        self.viewmenu.add_command(label="Zoom out", command=self.zoomoutmenu)
        self.viewmenu.add_separator()
        self.viewmenu.add_command(label="Shrink", command=self.shrink)
        self.viewmenu.add_command(label="Stretch", command=self.stretch)
        self.menubar.add_cascade(label="View", menu=self.viewmenu)
        self.graphmenu = Menu(self.menubar, tearoff=0)
        self.graphmenu.add_command(label="Find paths", command=self.findPaths)
        self.graphmenu.add_command(label="Write fasta", command=self.writeFasta)
        self.graphmenu.add_command(label="Write multifasta", command=self.writeMultiFasta)
        self.menubar.add_cascade(label="Selected", menu=self.graphmenu)
        self.currxscroll = 1000000
        self.curryscroll = 1000000
        self.fontsize=12
        self.customFont = tkFont.Font(family="Helvetica", size=self.fontsize)
        master.config(menu=self.menubar)
        self.panes = PanedWindow(root, orient=VERTICAL, sashrelief=SUNKEN, sashpad=5)
        self.panes.pack(fill=BOTH, expand=1)
        self.mainframe = Frame(self.panes)
        self.mainframe.grid_rowconfigure(0, weight=1)
        self.mainframe.grid_columnconfigure(0, weight=1)
        xscrollbar = Scrollbar(self.mainframe, orient=HORIZONTAL)
        xscrollbar.grid(row=1, column=0, sticky=E+W)
        yscrollbar = Scrollbar(self.mainframe)
        yscrollbar.grid(row=0, column=1, sticky=N+S)
        self.canvas = Canvas(self.mainframe, bd=0, bg='#FFFAF0', scrollregion=(0, 0, self.currxscroll, self.curryscroll),
                        xscrollcommand=xscrollbar.set,
                        yscrollcommand=yscrollbar.set)
        self.canvas.grid(row=0, column=0, sticky=N+S+E+W)
        xscrollbar.config(command=self.canvas.xview)
        yscrollbar.config(command=self.canvas.yview)
        self.panes.add(self.mainframe, sticky=N+S+E+W, minsize=30, width=500, height=400)
        self.utility = PanedWindow(self.panes)
        self.panes.add(self.utility, minsize=50)
        self.consoletext = StringVar(value='')
        self.console = Label(self.utility, bg='#FFFF99', relief=SUNKEN, textvariable=self.consoletext, anchor=SW, justify=LEFT)
        self.utility.add(self.console, sticky=NSEW, height=30, width=90, minsize=100)
        self.contiglist = Frame(self.utility)
        self.utility.add(self.contiglist, sticky=NSEW, minsize=120)
        self.dirframe = Frame(self.utility)
        self.utility.add(self.dirframe, sticky=NSEW, minsize=70)
        self.lengthframe = Frame(self.utility)
        self.utility.add(self.lengthframe, sticky=NSEW, minsize=70)
        self.contiglist.grid_rowconfigure(1, weight=1, minsize=30)
        self.contiglist.grid_columnconfigure(0, weight=1)
        self.dirframe.grid_rowconfigure(1, weight=1, minsize=30)
        self.dirframe.grid_columnconfigure(0, weight=1)
        self.lengthframe.grid_rowconfigure(1, weight=1, minsize=30)
        self.lengthframe.grid_columnconfigure(0, weight=1)
        self.clscroll = Scrollbar(self.lengthframe, orient=VERTICAL)
        self.namelabel = Label(self.contiglist, text='Contig name', bg='#FFFFFF',  relief=SUNKEN)
        self.namelabel.grid(row=0, column=0, sticky=EW)
        self.namelist = Listbox(self.contiglist, yscrollcommand=self.clscroll.set, exportselection=0, width=30, height=5)
        self.namelist.bind('<Button-1>', self.setselectedcontig)
        self.namelist.bind('<Double-Button-1>', self.doublecontig)
        self.namelist.bind('<MouseWheel>', self.onmousewheel)
        self.namelist.bind('<Button-4>', self.onmousewheel)
        self.namelist.bind('<Button-5>', self.onmousewheel)
        self.namelist.grid(row=1, column=0, sticky=NSEW)
        self.dirlabel = Label(self.dirframe, bg='#FFFFFF', text='Strand', relief=SUNKEN)
        self.dirlabel.grid(row=0, column=0, sticky=NSEW)
        self.dirlist = Listbox(self.dirframe, yscrollcommand=self.clscroll.set, exportselection=0, width=7, height=5)
        self.dirlist.bind('<Button-1>', self.setselectedcontig)
        self.dirlist.bind('<Double-Button-1>', self.doublecontig)
        self.dirlist.bind('<MouseWheel>', self.onmousewheel)
        self.dirlist.bind('<Button-4>', self.onmousewheel)
        self.dirlist.bind('<Button-5>', self.onmousewheel)
        self.dirlist.grid(row=1, column=0, sticky=NSEW)
        self.lengthlabel = Label(self.lengthframe, text='Length', bg='#FFFFFF',  relief=SUNKEN)
        self.lengthlabel.grid(row=0, column=0, sticky=EW)
        self.lengthlist = Listbox(self.lengthframe, yscrollcommand=self.clscroll.set, exportselection=0, width=7, height=5)
        self.lengthlist.bind('<Button-1>', self.setselectedcontig)
        self.lengthlist.bind('<Double-Button-1>', self.doublecontig)
        self.lengthlist.bind('<MouseWheel>', self.onmousewheel)
        self.lengthlist.bind('<Button-4>', self.onmousewheel)
        self.lengthlist.bind('<Button-5>', self.onmousewheel)
        self.lengthlist.grid(row=1, column=0, sticky=NSEW)
        self.clscroll = Scrollbar(self.lengthframe, orient=VERTICAL)
        self.clscroll.config(command=self.contiglistview)
        self.clscroll.grid(row=1, column=1, sticky=NS)
        master.geometry('+10+20')
        root.minsize(600, 400)
        self.canvas.tag_bind('map', '<B1-Motion>', self.dragMove)
        self.canvas.tag_bind('map', '<Button-1>', self.recordMark)
        self.canvas.tag_bind('map', '<Double-Button-1>', self.addtolist)
        self.rcmenu = Menu(root, tearoff=0)
        self.rcmenu.add_command(label="Reverse", command=self.reverse_contig)
        self.rcmenu.add_command(label="Remove", command=self.remove_contig)
        self.rcmenu.add_command(label="Duplicate", command=self.duplicate_contig)
        self.rcmenu.add_separator()
        self.rcmenu.add_command(label="Show to", command=self.show_to)
        self.rcmenu.add_command(label="Show from", command=self.show_from)
        self.canvas.tag_bind('map', '<Button-3>', self.rightClick)
        self.canvas.bind('<Button-2>', self.beginDrag)
        self.canvas.bind('<B2-Motion>', self.dragCanvas)
        root.bind('w', self.zoomin)
        root.bind('s', self.zoomout)
        root.bind('a', self.shrink)
        root.bind('d', self.stretch)
        root.bind('<Left>', self.scrollleft)
        root.bind('<Right>', self.scrollright)
        root.bind('<Up>', self.scrollup)
        root.bind('<Down>', self.scrolldown)
        self.canvas.bind('<Button-5>', self.zoomout)
        self.canvas.bind('<Button-4>', self.zoomin)
        self.canvas.bind('<MouseWheel>', self.zoomcanvas)
        self.canvas.bind('<Button-1>', self.removerc)
        self.csagfile = StringVar(value='')
        self.selfcomparefile = StringVar(value='')
        self.outfile = StringVar(value='')
        self.buffer = StringVar(value='nnnnnnnnnn')
        self.reffile = StringVar(value='/home/mitch/plas.fa')
        self.blastfile = StringVar(value='/home/mitch/temp/query_tempdb.out')
        self.contigfile = StringVar(value='/home/mitch/curr_proj/coif/UTI89-2_51_Contigs.fasta')
        self.csag = StringVar(value='')
        self.readfile = StringVar(value='/home/mitch/curr_proj/coif/UTI89-2_trimmed_ill.fastq')
        self.workingDir = StringVar(value='temp')
        try:
            os.makedirs(self.workingDir.get())
        except:
            pass
        self.minlength = IntVar(value=100)
        self.minlengthblast = IntVar(value=100)
        self.intra = IntVar(value=0)
        self.minlengthratio = DoubleVar(value=0.1)
        self.minident = DoubleVar(value=95.0)
        self.minbitscore = DoubleVar(value=0)
        self.maxtrim = IntVar(0)
        self.maxevalue = DoubleVar(value=0.05)
        self.getOverlap = IntVar(value=1)
        self.minoverlap = IntVar(value=30)
        self.maxmm = IntVar(value=2)
        self.getdb = IntVar(value=1)
        self.nmersize = IntVar(value=31)
        self.maxdist = IntVar(value=300)
        self.cutauto = IntVar(value=1)
        self.nmercut = IntVar(value=4)
        self.nmerave = IntVar(value=10)
        self.getPaired = IntVar(value=1)
        self.insertsize = IntVar(value=600)
        self.orientation = StringVar(value='--> <--')
        self.minrl = IntVar(value=75)
        self.minpairedge = IntVar(value=5)
        self.view = StringVar(value='BLAST')
        self.viewlist = StringVar(value='')
        self.viewref = IntVar(value=1)
        self.shorten = StringVar(value='No')
        self.getLong = IntVar(value=0)
        self.minlong = IntVar(value=2)
        self.minlongover = IntVar(value=100)
        self.minlongident = IntVar(value=90)
        self.visible = set()
        self.contigheight = 25
        self.scaledown = IntVar(value=100)
        self.maxbp = IntVar(value=5000)
        self.maxnode = IntVar(value=10)
        self.maxpath = IntVar(value=1000000)
        self.onlyshort = IntVar(value=1)
        self.originalxscroll = self.currxscroll
        self.originalyscroll = self.curryscroll
        self.refline = 50
        self.contigline = self.refline + self.contigheight + 200
        self.contigDict = {}
        self.hitlist = None
        self.leftmost = None
        self.rightmost = None
        self.selected = []
        self.newscaledown = self.scaledown.get()

    def addtolist(self, event):
        thetag = self.canvas.gettags(CURRENT)[0]
        thecontig = self.canvas.find_withtag(thetag)[0]
        thecolour = self.canvas.itemcget(thecontig, "fill")
        if thecolour == "#009989":
            self.canvas.itemconfig(thecontig, fill='#A3A948')
            contig = thetag[1:]
            if '_dup' in contig:
                dupno = int(contig.split('_dup')[1])
                contig = contig.split('_dup')[0]
            else:
                dupno = 0
            self.namelist.insert(END, contig)
            if self.contigDict[contig].orient[dupno]:
                self.dirlist.insert(END, '+')
            else:
                self.dirlist.insert(END, '-')
            self.lengthlist.insert(END, str(self.contigDict[contig].length))
            self.selected.append(thetag)
        else:
            self.canvas.itemconfig(thecontig, fill='#009989')
            for i in range(len(self.selected)):
                if self.selected[i] == thetag:
                    self.namelist.delete(i)
                    self.dirlist.delete(i)
                    self.lengthlist.delete(i)
                    self.selected.pop(i)
                    break



    def contiglistview(self, *args):
        apply(self.namelist.yview, args)
        apply(self.dirlist.yview, args)
        apply(self.lengthlist.yview, args)

    def onmousewheel(self, event):
        pass
    
    def setselectedcontig(self, event):
        pass
    
    def doublecontig(self, event):
        x = self.namelist.nearest(event.y)
        self.goto(self.namelist.get(x))


    def beginDrag(self, event):
        self.canvas.scan_mark(event.x, event.y)  # initial middle-mouse click
         
    def dragCanvas(self, event):
        self.canvas.scan_dragto(event.x, event.y, 1) # capture dragging

    def rightClick(self, event):
        self.rcmenu.post(event.x_root, event.y_root)
        self.rctag = self.canvas.gettags(CURRENT)[0]

    def removerc(self, event):
        self.rcmenu.unpost()

    def remove_contig(self):
        thetag = self.rctag
        starttag = thetag + 's'
        endtag = thetag + 'e'
        btag = thetag + 'b'
        self.visible.remove(thetag[1:])
        self.canvas.delete(thetag)
        self.canvas.delete(starttag)
        self.canvas.delete(endtag)
        self.canvas.delete(btag)

    def reverse_contig(self):
        thetag = self.rctag
        if '_dup' in thetag:
            contig = thetag.split('_dup')[0][1:]
            dupnum = int(thetag.split('_dup')[1])
        else:
            contig = thetag[1:]
            dupnum = 0
        if contig in self.contigDict:
            self.contigDict[contig].orient[dupnum] = not self.contigDict[contig].orient[dupnum]
        starttag = thetag + 's'
        bb = self.canvas.coords(thetag)
        start, end = bb[0], bb[2]
        arcs= self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            newx = start + end - x[0]
            self.canvas.coords(i, newx, x[1], (newx + x[4]) / 2, abs(newx - x[4]) / 4 + (x[1] + x[5]) / 2, x[4], x[5])
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            newx = start + end - x[4]
            self.canvas.coords(i, x[0], x[1], (x[0] + newx) / 2, abs(x[0] - newx) / 4 + (x[1] + x[5]) / 2, newx, x[5])
        htag = thetag + 'h'
        selfhits = self.canvas.find_withtag(htag)
        for i in selfhits:
            x = self.canvas.coords(i)
            newx1 = start + end - x[2]
            newx2 = start + end - x[0]
            self.canvas.coords(i, newx1, x[1], newx2, x[3])
            thecolour = self.canvas.itemcget(i, "fill")
            if thecolour == "#F85931":
                thecolour = "#EDB92E"
            else:
                thecolour = "#F85931"
            self.canvas.itemconfig(self.canvas.gettags(i)[-1], fill=thecolour)
        endtag = thetag + 'b'
        bhits = self.canvas.find_withtag(endtag)
        if endtag[0] == 'r':
            for i in bhits:
                x = self.canvas.coords(i)
                self.canvas.coords(i, end - (x[0] - start), x[1], end - (x[2] - start), x[3], x[4], x[5], x[6], x[7])
                thecolour = self.canvas.itemcget(i, "fill")
                if thecolour == "#F85931":
                    thecolour = "#EDB92E"
                else:
                    thecolour = "#F85931"
                self.canvas.itemconfig(i, fill=thecolour)
        else:
            for i in bhits:
                x = self.canvas.coords(i)
                self.canvas.coords(i, x[0], x[1], x[2], x[3], end - (x[4] - start), x[5], end - (x[6] - start), x[7])
                thecolour = self.canvas.itemcget(i, "fill")
                if thecolour == "#F85931":
                    thecolour = "#EDB92E"
                else:
                    thecolour = "#F85931"
                self.canvas.itemconfig(i, fill=thecolour)
        texttag = thetag + 't'
        try:
            textitem = self.canvas.find_withtag(texttag)[0]
            thetext = self.canvas.itemcget(textitem, "text")
            namelen = len(contig)
            if len(thetext) >= namelen + 2:
                start = thetext[:namelen + 1]
                strand = thetext[namelen + 1]
                end = thetext[namelen + 2:]
                if strand == '+':
                    strand = '-'
                else:
                    strand = '+'
                self.canvas.itemconfig(textitem, text=start + strand + end)
        except IndexError:
            pass




    def duplicate_contig(self):
        thetag = self.rctag
        contig = self.canvas.find_withtag(thetag)
        if self.rctag[0] == 'r':
            return
        if '_dup' in thetag:
            contigtag = thetag.split('_dup')[0][1:]
            dupnum = int(thetag.split('_dup')[1])
        else:
            contigtag = thetag[1:]
            dupnum = 0
        bb = self.canvas.coords(contig[0])
        newtag = 'c' + contigtag + '_dup' + str(self.contigDict[contigtag].dupnumber)
        self.visible.add(newtag[1:])
        self.contigDict[contigtag].dupnumber += 1
        self.contigDict[contigtag].orient.append(self.contigDict[contigtag].orient[dupnum])
        mod = int(1.5*self.contigheight)
        btag = thetag + 'b'
        bhits = self.canvas.find_withtag(btag)
        if btag[0] == 'r':
            pass
        else:
            for i in bhits:
                x = self.canvas.coords(i)
                colour = self.canvas.itemcget(i, "fill")
                self.canvas.create_polygon(x[0], x[1], x[2], x[3], x[4] + mod, x[5] + mod, x[6] + mod, x[7] + mod, \
                                           fill=colour, outline="black", tags=(newtag + 'b', self.canvas.gettags(i)[1], 'blast'))
        self.canvas.create_rectangle(bb[0] + mod, bb[1] + mod, bb[2] + mod, bb[3] + mod, fill='#009989', tags=(newtag, 'contig', 'map'))
        htag = thetag + 'h'
        selfhits = self.canvas.find_withtag(htag)
        for i in selfhits:
            colour = self.canvas.itemcget(i, "fill")
            x = self.canvas.coords(i)
            self.canvas.create_rectangle(x[0] + mod, x[1] + mod, x[2] + mod, x[3] + mod, fill=colour, tags=(newtag, 'sblast', 'map', newtag + 'h'))
        starttag = thetag + 's'
        arcs= self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.create_line(x[0] + mod, x[1] + mod, (x[0] + mod + x[4]) / 2, abs(x[0] + mod - x[4]) / 4 + (x[1] + mod + x[5]) / 2, x[4], x[5],\
                                     smooth=True, width=3, tags=(newtag + 's', self.canvas.gettags(i)[1] , 'arc'))
        endtag = thetag + 'e'
        arcs= self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.create_line(x[0], x[1], (x[0] + x[4] + mod) / 2, abs(x[0] - x[4] - mod) / 4 + (x[1] + x[5] + mod) / 2, x[4] + mod, x[5] + mod,\
                                    smooth=True, width=3, tags=(self.canvas.gettags(i)[0], newtag + 'e', 'arc'))
        texttag = thetag + 't'
        try:
            textitem = self.canvas.find_withtag(texttag)[0]
            thetext = self.canvas.itemcget(textitem, "text")
            x = self.canvas.coords(textitem)
            self.canvas.create_text(x[0] + mod, x[1] + mod, fill='white', font=self.customFont,
                                    anchor=W, text=thetext, tags=(newtag, 'map', 'text', newtag + 't'))
        except IndexError:
            pass

    def add_contig(self, i):
        if not i in self.visible:
            self.visible.add(i)
            placed = False
            x1 = self.contigDict[i].xpos
            y1 = self.contigDict[i].ypos
            x2 = self.contigDict[i].xpos + self.contigDict[i].xlength
            y2 = self.contigDict[i].ypos + self.contigheight
            while not placed:
                placed = True
                for j in self.canvas.find_overlapping(x1, y1, x2, y2):
                    if 'contig' in self.canvas.gettags(j):
                        placed = False
                        y1 += 30
                        y2 += 30
                        if y2 >= self.curryscroll:
                            y1 = self.contigDict[i].ypos + 10
                            y2 = self.contigDict[i].ypos + self.contigheight + 10
                            placed = True
                        break
            self.canvas.create_rectangle(x1, y1, x2, y2, fill="#009989", tags=('c' + i, 'contig', 'map'))
            if self.contigDict[i].orient[0]:
                dir = '+'
            else:
                dir = '-'
            thetext = i + ' ' + dir + ' ' + self.contigDict[i].strlen
            text = self.canvas.create_text(x1 + 2, y1 + self.contigheight/2, fill='white', font=self.customFont,
                                           anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
            if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                self.canvas.delete(text)
                thetext = i + ' ' + dir
                text = self.canvas.create_text(x1 + 2, y1 + self.contigheight/2, fill='white', font=self.customFont,
                                               anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                    self.canvas.delete(text)
                    thetext = i
                    text = self.canvas.create_text(x1 + 2, y1 + self.contigheight/2, fill='white', font=self.customFont,
                                                   anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                    if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                        self.canvas.delete(text)
            for j in self.contigDict[i].to:
                if j[0] in self.visible:
                    tocoords = self.canvas.coords('c' + j[0])
                    startx = x2
                    starty = y1 + self.contigheight /2
                    endy = (tocoords[1] + tocoords[3]) /2
                    if (j[1] and self.contigDict[j[0]].orient[0]) or (not j[1] and not self.contigDict[j[0]].orient[0]):
                        endx = tocoords[0]
                    else:
                        endx = tocoords[2]
                    self.canvas.create_line(startx, starty, (startx + endx) / 2, abs(startx - endx) /4 + (starty + endy) / 2, endx, endy,\
                                            smooth=True, width=3, tags=('c' + i + 's', 'c' + j[0] + 'e', 'arc'))
                    if self.contigDict[j[0]].dupnumber != 1:
                        for k in range(1, int(self.contigDict[j[0]].dupnumber)):
                            newtag = j[0] + '_dup' + str(k)
                            if newtag in self.visible:
                                tocoords = self.canvas.coords('c' + newtag)
                                endy = (tocoords[1] + tocoords[3]) /2
                                if (j[1] and self.contigDict[j[0]].orient[k]) or (not j[1] and not self.contigDict[j[0]].orient[k]):
                                    endx = tocoords[0]
                                else:
                                    endx = tocoords[2]
                                self.canvas.create_line(startx, starty, (startx + endx) / 2, abs(startx - endx) /4 + (starty + endy) / 2, endx, endy,
                                                        smooth=True, width=3, tags=('c' + i + 's', 'c' + newtag + 'e', 'arc'))
            for j in self.contigDict[i].fr:
                if j[0] in self.visible:
                    tocoords = self.canvas.coords('c' + j[0])
                    startx = x1
                    starty = y1 + self.contigheight /2
                    endy = (tocoords[1] + tocoords[3]) /2
                    if (j[1] and self.contigDict[j[0]].orient[0]) or (not j[1] and not self.contigDict[j[0]].orient[0]):
                        endx = tocoords[0]
                    else:
                        endx = tocoords[2]
                    self.canvas.create_line(startx, starty, (startx + endx) / 2, abs(startx - endx) /4 + (starty + endy) / 2, endx, endy,\
                                            smooth=True, width=3, tags=('c' + i + 's', 'c' + j[0] + 'e', 'arc'))
                    if self.contigDict[j[0]].dupnumber != 1:
                        for k in range(1, int(self.contigDict[j[0]].dupnumber)):
                            newtag = j[0] + '_dup' + str(k)
                            if newtag in self.visible:
                                tocoords = self.canvas.coords('c' + newtag)
                                endy = (tocoords[1] + tocoords[3]) /2
                                if (j[1] and self.contigDict[j[0]].orient[k]) or (not j[1] and not self.contigDict[j[0]].orient[k]):
                                    endx = tocoords[0]
                                else:
                                    endx = tocoords[2]
                                self.canvas.create_line(startx, starty, (startx + endx) / 2, abs(startx - endx) /4 + (starty + endy) / 2, endx, endy,
                                                        smooth=True, width=3, tags=('c' + i + 's', 'c' + newtag + 'e', 'arc'))

    def add_contig_dialogue(self):
        contig = tkSimpleDialog.askstring('Add Contig', 'Contig name')
        if contig is None:
            return
        if contig in self.visible:
            self.goto(contig)
        else:
            self.contigDict[contig].visible = True
            self.contigDict[contig].xpos = 50
            self.contigDict[contig].ypos = 275
            self.contigDict[contig].xlength = self.contigDict[contig].length / self.newscaledown
            self.add_contig(contig)
            self.canvas.xview_moveto(0)
            self.canvas.yview_moveto(0)


    def show_to(self):
        thetag = self.rctag
        if '_dup' in thetag:
            contig = thetag.split('_dup')[0][1:]
            dupnum = int(thetag.split('_dup')[1])
        else:
            contig = thetag[1:]
            dupnum = 0
        coords = self.canvas.coords(thetag)
        startx = coords[2] + 10
        starty = coords[1]
        if self.contigDict[contig].orient[dupnum]:
            for i in self.contigDict[contig].to:
                if not i[0] in self.visible:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                    starty += 35
        else:
            for i in self.contigDict[contig].fr:
                if not i[0] in self.visible:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                    starty += 35



    def show_from(self):
        thetag = self.rctag
        if '_dup' in thetag:
            contig = thetag.split('_dup')[0][1:]
            dupnum = int(thetag.split('_dup')[1])
        else:
            contig = thetag[1:]
            dupnum = 0
        coords = self.canvas.coords(thetag)
        startx = coords[0]
        starty = coords[1]
        if not self.contigDict[contig].orient[dupnum]:
            for i in self.contigDict[contig].to:
                if not i[0] in self.visible:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.contigDict[i[0]].xpos = startx - 10 - self.contigDict[i[0]].xlength
                    self.add_contig(contig)
                    starty += 35
        else:
            for i in self.contigDict[contig].fr:
                if not i[0] in self.visible:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.contigDict[i[0]].xpos = startx - 10 - self.contigDict[i[0]].xlength
                    self.add_contig(i[0])
                    starty += 35
    
    def find_contig(self):
        contig = tkSimpleDialog.askstring('Find Contig', 'Contig name')
        if contig != None:
            self.goto(contig)

    def goto(self, contig):
        try:
            x = self.canvas.coords('c' + contig)
            self.canvas.xview_moveto(max([0, x[0] - 200]) / self.currxscroll)
            self.canvas.yview_moveto(max([0, x[1] - 200]) / self.curryscroll)
        except:
            pass

    def zoominmenu(self):
        self.canvas.scale(ALL, 0, 0, 1.1, 1.1)
        self.currxscroll *= 1.1
        self.curryscroll *= 1.1
        self.contigheight *= 1.1
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.fontsize *= 1.1
        self.customFont.configure(size=int(round(self.fontsize)))


    def zoomoutmenu(self):
        if self.fontsize < 1:
            pass
        else:
            self.canvas.scale(ALL, 0, 0, 0.909090909, 0.909090909)
            self.currxscroll *= 0.909090909
            self.curryscroll *= 0.909090909
            self.contigheight *= 0.909090909
            self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
            self.fontsize *= 0.909090909
            self.customFont.configure(size=int(round(self.fontsize)))

    def zoomin(self, event):
        self.canvas.scale(ALL, 0, 0, 1.1, 1.1)
        self.currxscroll *= 1.1
        self.curryscroll *= 1.1
        self.contigheight *= 1.1
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.canvas.xview_moveto((self.canvas.canvasx(event.x) * 1.1 - event.x) / self.currxscroll)
        self.canvas.yview_moveto((self.canvas.canvasy(event.y) * 1.1 - event.y) / self.curryscroll)
        self.fontsize *= 1.1
        self.customFont.configure(size=int(round(self.fontsize)))

    def zoomout(self, event):
        if self.fontsize < 1:
            pass
        else:
            self.canvas.scale(ALL, 0, 0, 0.909090909, 0.909090909)
            self.currxscroll *= 0.909090909
            self.curryscroll *= 0.909090909
            self.contigheight *= 0.909090909
            self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
            self.canvas.xview_moveto((self.canvas.canvasx(event.x) * 0.909090909 - event.x) / self.currxscroll)
            self.canvas.yview_moveto((self.canvas.canvasy(event.y) * 0.909090909 - event.y) / self.curryscroll)
            self.fontsize *= 0.909090909
            self.customFont.configure(size=int(round(self.fontsize)))

    def shrink(self, event=None):
        if self.newscaledown > 10 * self.scaledown.get():
            pass
        else:
            self.canvas.scale(ALL, 0, 0, 0.909090909, 1)
            self.currxscroll *= 0.909090909
            self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
            self.newscaledown /= 0.909090909
            textitems = self.canvas.find_withtag('text')
            for z in textitems:
                contigtag = self.canvas.gettags(z)[0]
                if '_dup' in contigtag:
                    i = contigtag.split('_dup')[0][1:]
                    dupnum = int(contigtag.split('_dup')[1])
                else:
                    i = contigtag[1:]
                    dupnum = 0
                if self.contigDict[i].orient[dupnum]:
                    dir = '+'
                else:
                    dir = '-'
                textcoords = self.canvas.coords(z)
                textend = self.canvas.bbox(z)[2]
                contigend = self.canvas.coords('c' + i)[2]
                if textend >= contigend - 2:
                    self.canvas.delete(z)
                    thetext = i + ' ' + dir
                    text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                   anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                    textend = self.canvas.bbox(text)[2]
                    if textend >= contigend - 2:
                        self.canvas.delete(text)
                        thetext = i
                        text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                       anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                        textend = self.canvas.bbox(text)[2]
                        if textend >= contigend:
                            self.canvas.delete(text)


    def stretch(self, event=None):
        if self.newscaledown < 0.1 * self.scaledown.get():
            pass
        else:
            self.canvas.scale(ALL, 0, 0, 1.1, 1)
            self.currxscroll *= 1.1
            self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
            self.newscaledown /= 1.1
            textitems = self.canvas.find_withtag('contig')
            self.canvas.delete('text')
            for z in textitems:
                contigtag = self.canvas.gettags(z)[0]
                if '_dup' in contigtag:
                    i = contigtag.split('_dup')[0][1:]
                    dupnum = int(contigtag.split('_dup')[1])
                else:
                    i = contigtag[1:]
                    dupnum = 0
                if self.contigDict[i].orient[dupnum]:
                    dir = '+'
                else:
                    dir = '-'
                length = self.contigDict[i].strlen
                textcoords = (self.canvas.coords(z)[0] + 2, (self.canvas.coords(z)[1] + self.canvas.coords(z)[3]) /2)
                contigend = self.canvas.coords('c' + i)[2]
                thetext = i + ' ' + dir + ' ' + length
                text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                textend = self.canvas.bbox(text)[2]
                if textend >= contigend - 2:
                    self.canvas.delete(text)
                    thetext = i + ' ' + dir
                    text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                   anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                    textend = self.canvas.bbox(text)[2]
                    if textend >= contigend - 2:
                        self.canvas.delete(text)
                        thetext = i
                        text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                       anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                        textend = self.canvas.bbox(text)[2]
                        if textend >= contigend - 2:
                            self.canvas.delete(text)


    def zoomcanvas(self, event):
        pass

    def scrollleft(self, event):
        self.canvas.xview_scroll(-2, 'units')
    
    def scrollright(self, event):
        self.canvas.xview_scroll(2, 'units')
    
    def scrollup(self, event):
        self.canvas.yview_scroll(-2, 'units')
    
    def scrolldown(self, event):
        self.canvas.yview_scroll(2, 'units')

    def recordMark(self, event):
        self.oldx = self.canvas.canvasx(event.x)
        self.oldy = self.canvas.canvasy(event.y)

    def dragMove(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        diffx = x - self.oldx
        diffy = y - self.oldy
        self.oldx = x
        self.oldy = y
        thetag = self.canvas.gettags(CURRENT)[0]
        self.canvas.move(thetag, diffx, diffy)
        starttag = thetag + 's'
        arcs= self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0] + diffx, x[1] + diffy, x[2] + diffx/2, abs(x[0] + diffx - x[4]) / 4 + (x[1] + diffy + x[5]) / 2, x[4], x[5])
        endtag = thetag + 'e'
        arcs= self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0], x[1], x[2] + diffx/2, abs(x[0] - x[4] + diffx) / 4 + (x[1] + x[5] + diffy) / 2, x[4] + diffx, x[5] + diffy)
        endtag = thetag + 'b'
        arcs= self.canvas.find_withtag(endtag)
        if endtag[0] == 'r':
            for i in arcs:
                x = self.canvas.coords(i)
                self.canvas.coords(i, x[0] + diffx, x[1] + diffy, x[2] + diffx, x[3] + diffy, x[4], x[5], x[6], x[7])
        else:
            for i in arcs:
                x = self.canvas.coords(i)
                self.canvas.coords(i, x[0], x[1], x[2], x[3], x[4] + diffx, x[5] + diffy, x[6] + diffx, x[7] + diffy)

    def update_console(self, text):
        self.consoletext.set(self.consoletext.get() + '\n' + text)

    def clear_all(self):
        self.canvas.delete(ALL)
        self.hitlist = []
        self.edgelist = []
        self.contigDict = {}
        self.visible = set()
        self.currxscroll = self.originalxscroll
        self.curryscroll = self.originalyscroll
        self.contigheight = 25

    def load_assembly(self):
        filename = tkFileDialog.askopenfilename()
        if filename == '':
            return
        self.clear_all()
        self.csagfile.set(filename)
        whatsthis = open(self.csagfile.get())
        what = whatsthis.readline()
        whatsthis.close()
        self.contigDict = {}
        self.edgelist = []
        if what[0] == '>':
            self.load_fasta()
            self.update_console('FASTA loaded.')
        elif what.split()[0] == 'NODE':
            self.load_csag()
            self.update_console('CSAG loaded.')
        elif len(what.split()) == '4' and what.split()[0].isdigit() and what.split()[1].isdigit()\
            and what.split()[2].isdigit() and what.split()[3].isdigit():
            self.load_lastgraph()
            self.update_console('LastGraph loaded.')
        else:
            tkMessageBox.showerror('Invalid format', 'Contiguity cannot recognise file type.')
        self.writeWorkCont()

    def load_fasta(self):
        fastafile = open(self.csagfile.get())
        namelist = []
        maxlen = 0
        for line in fastafile:
            if line.startswith('>'):
                namelist.append(line.rstrip())
                if maxlen < len(namelist[-1]):
                    maxlen = len(namelist[-1])
        for i in range(maxlen):
            startletter = namelist[0][i]
            same = True
            for j in namelist[1:]:
                if j[i] != startletter:
                    same = False
                    break
            if not same:
                cuta = i
                break
            else:
                cuta = 0
        for i in range(1, maxlen):
            startletter = namelist[0][-i]
            same = True
            for j in namelist[1:]:
                if j[-i] != startletter:
                    same = False
                    break
            if not same:
                cutb = -i + 1
                break
            else:
                cutb = None
        fastafile.close()
        fastafile = open(self.csagfile.get())
        first = True
        for line in fastafile:
            if line.startswith('>'):
                if first:
                    first = False
                else:
                    aninstance = contig(name, seq)
                    self.contigDict[name] = aninstance
                name = line.rstrip()[cuta:cutb]
                seq = ''
            else:
                seq += line.rstrip()
        if name[:5] == 'NODE_':
            name = name.split('_')[1]
        aninstance = contig(name, seq)
        self.contigDict[name] = aninstance
        fastafile.close()

    def load_csag(self):
        csag = open(self.csagfile.get())
        for line in csag:
            if line.split()[0] == 'NODE':
                title, name, seq = line.split()
                aninstance = contig(name, seq)
                self.contigDict[name] = aninstance
            elif line.split()[0] == 'EDGE':
                title, n1, d1, n2, d2, overlap = line.split()
                if d1 == 'True':
                    d1 = True
                else:
                    d1 = False
                if d2 == 'True':
                    d2 = True
                else:
                    d2 = False
                if overlap.isdigit():
                    overlap = int(overlap)
                self.edgelist.append((n1, d1, n2, d2, overlap))
        for i in self.edgelist:
            contiga, dira, contigb, dirb, overlap = i
            if dira and dirb:
                self.contigDict[contiga].to.append((contigb, True, overlap))
                self.contigDict[contigb].fr.append((contiga, False, overlap))
            elif dira and not dirb:
                self.contigDict[contiga].to.append((contigb, False, overlap))
                self.contigDict[contigb].to.append((contiga, False, overlap))
            elif not dira and dirb:
                self.contigDict[contiga].fr.append((contigb, True, overlap))
                self.contigDict[contigb].fr.append((contiga, True, overlap))
            else:
                self.contigDict[contiga].fr.append((contigb, False, overlap))
                self.contigDict[contigb].to.append((contiga, True, overlap))


    def loadref(self):
        filename = tkFileDialog.askopenfilename(parent=self.blast_options)
        if filename == '':
            return
        self.reffile.set(filename)

    def loadselfcomp(self):
        filename = tkFileDialog.askopenfilename(parent=self.blast_options)
        if filename == '':
            return
        self.selfcomparefile.set(filename)

    def loadblast(self):
        filename = tkFileDialog.askopenfilename(parent=self.blast_options)
        if filename == '':
            return
        self.blastfile.set(filename)
        
    def loadcontig(self):
        filename = tkFileDialog.askopenfilename(parent=self.create_edges_top)
        if filename == '':
            return
        self.contigfile.set(filename)
        
    def loadview(self):
        filename = tkFileDialog.askopenfilename(parent=self.view_options)
        if filename == '':
            return
        self.viewlist.set(filename)

    def loadread(self):
        filename = tkFileDialog.askopenfilename(parent=self.self.create_edges_top)
        if filename == '':
            return
        self.readfile.set(filename)

    def loadwork(self):
        filename = tkFileDialog.askdirectory(parent=self.load_files)
        if filename == '':
            root.quit()
        self.workingDir.set(filename)

    def loadoutfile(self):
        filename = tkFileDialog.asksaveasfilename(parent=self.write_fasta)
        if filename == '':
            return
        else:
            self.outfile.set(filename)

    def create_edges(self):
        try:
            self.create_edges_top.destroy()
        except:
            pass
        self.create_edges_top = Toplevel()
        self.create_edges_top.geometry('+20+30')
        self.create_edges_top.title('Create CSAG')
        self.create_edges = Frame(self.create_edges_top, padx=5, pady=5)
        self.contigfilelabel = Label(self.create_edges, text='Contig file:')
        self.contigfilelabel.grid(row=0, column=0)
        self.contigfileentry = Entry(self.create_edges, textvariable=self.contigfile)
        self.contigfileentry.grid(row=0, column=1, sticky=E)
        self.contigfileentrybutton = Button(self.create_edges, text='...', command=self.loadcontig)
        self.contigfileentrybutton.grid(row=0, column=2)
        self.readfilelabel = Label(self.create_edges, text='Read file:')
        self.readfilelabel.grid(row=1, column=0)
        self.readfileentry = Entry(self.create_edges, textvariable=self.readfile)
        self.readfileentry.grid(row=1, column=1, sticky=E)
        self.readfileentrybutton = Button(self.create_edges, text='...', command=self.loadread)
        self.readfileentrybutton.grid(row=1, column=2)
        self.create_edges.grid(row=0, column=0)
        self.overlapEdges = Frame(self.create_edges, relief=SUNKEN, bd=2)
        self.overlapEdges.grid(row=2, column=0, padx=5, pady=3, sticky=EW, columnspan=3)
        self.getOverlapLabel = Label(self.overlapEdges, text='Get overlapping edges:')
        self.getOverlapLabel.grid(row=0, column=0)
        self.getOverlapCheck = Checkbutton(self.overlapEdges, variable=self.getOverlap)
        self.getOverlapCheck.grid(row=0, column=1, sticky=W)
        self.minoverlaplabel = Label(self.overlapEdges, text='Minimum overlap length:')
        self.minoverlaplabel.grid(row=1, column=0)
        self.minoverlapentry = Entry(self.overlapEdges, textvariable=self.minoverlap)
        self.minoverlapentry.grid(row=1, column=1, sticky=W)
        self.maxmmlabel = Label(self.overlapEdges, text='Maximum mismatches in overlap:')
        self.maxmmlabel.grid(row=2, column=0)
        self.maxmmentry = Entry(self.overlapEdges, textvariable=self.maxmm)
        self.maxmmentry.grid(row=2, column=1, sticky=W)
        self.dbEdges = Frame(self.create_edges, relief=SUNKEN, bd=2)
        self.dbEdges.grid(row=3, column=0, padx=5, pady=3, sticky=EW, columnspan=3)
        self.getdblabel = Label(self.dbEdges, text='Get de bruijn edges:')
        self.getdblabel.grid(row=0, column=0)
        self.getdbentry = Checkbutton(self.dbEdges, variable=self.getdb)
        self.getdbentry.grid(row=0, column=1, sticky=W)
        self.nmersizelabel = Label(self.dbEdges, text='Nmer size:')
        self.nmersizelabel.grid(row=1, column=0)
        self.nmersizeentry = Entry(self.dbEdges, textvariable=self.nmersize)
        self.nmersizeentry.grid(row=1, column=1, sticky=EW)
        self.maxdistlabel = Label(self.dbEdges, text='Max. distance:')
        self.maxdistlabel.grid(row=2, column=0)
        self.maxdistentry = Entry(self.dbEdges, textvariable=self.maxdist)
        self.maxdistentry.grid(row=2, column=1, sticky=W)
        self.cutautolabel = Label(self.dbEdges, text='Auto detect cutoffs:')
        self.cutautolabel.grid(row=3, column=0)
        self.cutautoentry = Checkbutton(self.dbEdges, variable=self.cutauto)
        self.cutautoentry.grid(row=3, column=1, sticky=W)
        self.nmercutlabel = Label(self.dbEdges, text='Nmer cutoff:')
        self.nmercutlabel.grid(row=4, column=0)
        self.nmercutentry = Entry(self.dbEdges, textvariable=self.nmercut)
        self.nmercutentry.grid(row=4, column=1, sticky=W)
        self.nmeravelabel = Label(self.dbEdges, text='Nmer average:')
        self.nmeravelabel.grid(row=5, column=0)
        self.nmeraveentry = Entry(self.dbEdges, textvariable = self.nmerave)
        self.nmeraveentry.grid(row=5, column=1, sticky=W)
        self.pairedEdges = Frame(self.create_edges, relief=SUNKEN, bd=2)
        self.pairedEdges.grid(row=4, column=0, padx=5, pady=3, sticky=NSEW, columnspan=3)
        self.getPairedlabel = Label(self.pairedEdges, text='Get edges using paired reads:')
        self.getPairedlabel.grid(row=0, column=0)
        self.getPairedentry = Checkbutton(self.pairedEdges, variable=self.getPaired)
        self.getPairedentry.grid(row=0, column=1, sticky=W)
        self.insertsizelabel = Label(self.pairedEdges, text='Max. insert size:')
        self.insertsizelabel.grid(row=1, column=0)
        self.insertsizeentry = Entry(self.pairedEdges, textvariable=self.insertsize)
        self.insertsizeentry.grid(row=1, column=1, sticky=W)
        self.orientationlabel = Label(self.pairedEdges, text='Read orientation:')
        self.orientationlabel.grid(row=2, column=0)
        self.orientationentry = OptionMenu(self.pairedEdges, self.orientation, '--> <--', '<-- -->', '--> -->')
        self.orientationentry.grid(row=2, column=1, sticky=W)
        self.minrllabel = Label(self.pairedEdges, text='Minimum read length:')
        self.minrllabel.grid(row=3, column=0)
        self.minrlentry = Entry(self.pairedEdges, textvariable=self.minrl)
        self.minrlentry.grid(row=3, column=1, sticky=W)
        self.minpairedgelabel = Label(self.pairedEdges, text='Min. reads for edge:')
        self.minpairedgelabel.grid(row=4, column=0)
        self.minpairedgeentry = Entry(self.pairedEdges, textvariable=self.minpairedge)
        self.minpairedgeentry.grid(row=4, column=1, sticky=W)
        self.longEdges = Frame(self.create_edges, relief=SUNKEN, bd=2)
        self.longEdges.grid(row=7, column=0, padx=5, pady=3, sticky=NSEW, columnspan=3)
        self.getLonglabel = Label(self.longEdges, text='Get edges using long reads:')
        self.getLonglabel.grid(row=0, column=0)
        self.getLongentry = Checkbutton(self.longEdges, variable=self.getLong)
        self.getLongentry.grid(row=0, column=1, sticky=W)
        self.minlonglabel = Label(self.longEdges, text='Min. read number:')
        self.minlonglabel.grid(row=1, column=0)
        self.minlongentry = Entry(self.longEdges, textvariable=self.minlong)
        self.minlongentry.grid(row=1, column=1, sticky=W)
        self.minlongoverlabel = Label(self.longEdges, text='Min. overlap:')
        self.minlongoverlabel.grid(row=2, column=0)
        self.minlongoverentry = Entry(self.longEdges, textvariable=self.minlongover)
        self.minlongoverentry.grid(row=2, column=1, sticky=W)
        self.minlongidentlabel = Label(self.longEdges, text='Min. Identity:')
        self.minlongidentlabel.grid(row=3, column=0)
        self.minlongidententry = Entry(self.longEdges, textvariable=self.minlongident)
        self.minlongidententry.grid(row=3, column=1, sticky=W)
        self.okedges = Button(self.create_edges, text='Ok', command=self.ok_edges)
        self.okedges.grid(row=8, column=0, sticky=E)

    def ok_edges(self):
        if (not os.path.exists(self.readfile.get()) and (self.getdb.get() == 1 or self.getPaired.get() == 1 or self.getLong.get() == 1))\
          or not os.path.exists(self.contigfile.get()):
            tkMessageBox.showerror('File missing', 'Please add a contig and read file.')
            return
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        self.create_edges_top.destroy()
        self.edgeconsole = []
        self.lasteclen = 0
        self.queue = Queue.Queue()
        self.thethread = threading.Thread(target=self.edge_thread)
        self.thethread.start()
        self.update_edges()

    def update_edges(self):
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text != 'Creating CSAG finished.':
                    root.after(1000, self.update_edges)
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            self.update_console('Edge creation failed, please check console output.')
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_edges)
            return


    def edge_thread(self):
        self.edgelist = []
        self.contigDict = {}
        self.populateContigDict()
        self.queue.put('Creating edges now..')
        self.writeWorkCont()
        if self.getOverlap.get() == 1:
            count = len(self.edgelist)
            self.get_overlap_edges()
            count = len(self.edgelist) - count
            self.queue.put('Overlaps found ' + str(count) + ' edges.')
        if self.getdb.get() == 1:
            count = len(self.edgelist)
            self.get_nmer_freq()
            self.read_nmer_freq()
            if self.cutauto == 1:
                self.getnmercut()
            self.get_db_edges()
            count = len(self.edgelist) - count
            self.queue.put('De Bruijn found ' + str(count) + ' edges.')
        if self.getPaired.get() == 1:
            count = len(self.edgelist)
            self.get_paired_edge()
            count = len(self.edgelist) - count
            self.queue.put('Paired end found ' + str(count) + ' edges.')
        if self.getLong.get() == 1:
            count = len(self.edgelist)
            self.get_long_edge()
            count = len(self.edgelist) - count
            self.queue.put('Long reads found '+ str(count) + 'edges.')
        self.removeDuplicates()
        self.untangleEdges()
        self.queue.put(str(len(self.edgelist)) + ' edges found.')
        for i in self.edgelist:
            contiga, dira, contigb, dirb, overlap = i
            if dira and dirb:
                self.contigDict[contiga].to.append((contigb, True, overlap))
                self.contigDict[contigb].fr.append((contiga, False, overlap))
            elif dira and not dirb:
                self.contigDict[contiga].to.append((contigb, False, overlap))
                self.contigDict[contigb].to.append((contiga, False, overlap))
            elif not dira and dirb:
                self.contigDict[contiga].fr.append((contigb, True, overlap))
                self.contigDict[contigb].fr.append((contiga, True, overlap))
            else:
                self.contigDict[contiga].fr.append((contigb, False, overlap))
                self.contigDict[contigb].to.append((contiga, True, overlap))
        self.writeCSAG()
        self.csag.set(self.workingDir.get() + '/CSAG.txt')
        self.queue.put('Creating CSAG finished.')

    def populateContigDict(self):
        fastafile = open(self.contigfile.get())
        namelist = []
        maxlen = 0
        getvel = False
        for line in fastafile:
            if line.startswith('>'):
                namelist.append(line.split()[0])
                if maxlen < len(namelist[-1]):
                    maxlen = len(namelist[-1])
        splitlist = namelist[0].split('_')
        if len(splitlist) > 2 and splitlist[0] == '>NODE' and splitlist[1].isdigit():
            getvel = True
        else:
            for i in range(maxlen):
                startletter = namelist[0][i]
                same = True
                for j in namelist[1:]:
                    if j[i] != startletter:
                        same = False
                        break
                if not same:
                    cuta = i
                    break
                else:
                    cuta = 0
            for i in range(1, maxlen):
                startletter = namelist[0][-i]
                same = True
                for j in namelist[1:]:
                    if j[-i] != startletter:
                        same = False
                        break
                if not same:
                    cutb = -i + 1
                    if cutb == 0:
                        cutb = None
                    break
                else:
                    cutb = None
        fastafile.close()
        fastafile = open(self.contigfile.get())
        first = True
        for line in fastafile:
            if line.startswith('>'):
                if first:
                    first = False
                else:
                    aninstance = contig(name, seq)
                    self.contigDict[name] = aninstance
                if getvel:
                    name = line.split('_')[1]
                else:
                    name = line.split()[0][cuta:cutb]
                seq = ''
            else:
                seq += line.rstrip()
        if name[:5] == 'NODE_':
            name = name.split('_')[1]
        aninstance = contig(name, seq)
        self.contigDict[name] = aninstance

    def get_overlap_edges(self):
        minoverlap, maxmm, maxtrim = self.minoverlap.get(), self.maxmm.get(), self.maxtrim.get()
        stderr = open(self.workingDir.get() + '/bwaerr.txt', 'wa')
        subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/contigdb -in ' +
                         self.workingDir.get() + '/contigs.fa', shell=True, stdout=stderr).wait()
        subprocess.Popen('blastn -db ' + self.workingDir.get() + '/contigdb -outfmt 6 -query ' +
                         self.workingDir.get() + '/contigs.fa -out ' + self.workingDir.get() + '/contigscontigs.out', shell=True).wait()
        stderr.close()
        blastfile = open(self.workingDir.get() + '/contigscontigs.out')
        for line in blastfile:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm, indel = map(int, [qstart, qstop, rstart, rstop, length, mm, indel])
            mm = mm + indel
            if length <= minoverlap + 10:
                maxmmrev = maxmm
            else:
                maxmmrev = int(0.15 * length)
            if qstart - 1 <= maxmmrev and rstop - 1 <= maxmmrev and length + qstart - 2 + rstop >= minoverlap and mm + qstart + rstop - 2 <= maxmmrev \
              and self.contigDict[query].length > qstop + 20 and self.contigDict[subject].length > rstart + 20:
                self.edgelist.append((query, False, subject, True, rstart + qstart - 1))
            if qstart - 1 <= maxmmrev and rstop >= self.contigDict[subject].length - maxmmrev and length + qstart + self.contigDict[subject].length- rstop - 1 >= minoverlap \
              and mm + qstart + self.contigDict[subject].length- rstop - 1 <= maxmmrev \
              and self.contigDict[query].length > qstop + 20 and rstart > 20:
                self.edgelist.append((query, False, subject, False, self.contigDict[subject].length- rstart + qstart))
            if qstop >= self.contigDict[query].length - maxmmrev and rstart - 1 <= maxmmrev and length + qstop - self.contigDict[query].length + rstart - 1 >= minoverlap \
              and mm + rstart + qstop - self.contigDict[query].length - 1 <= maxmmrev \
              and qstart > 20 and self.contigDict[subject].length > rstop + 20:
                self.edgelist.append((query, True, subject, True, rstop + qstop - self.contigDict[query].length))
            if qstop >= self.contigDict[query].length - maxmmrev and rstart >= self.contigDict[subject].length - maxmmrev and length + qstop - self.contigDict[query].length \
              + rstart - self.contigDict[subject].length >= minoverlap and mm + qstop - self.contigDict[query].length + rstart - self.contigDict[subject].length <= maxmmrev \
              and qstart > 20 and rstop > 20:
                self.edgelist.append((query, True, subject, False, self.contigDict[subject].length - rstop + qstop - self.contigDict[query].length + 1))

    def dbpath(self, currpath, endnmer):
        nmersize, nmercut, maxdist, nmerave = self.nmersize.get(), self.nmercut.get(), self.maxdist.get(), self.nmerave.get()
        todo = [currpath]
        paths = []
        temptodo = todo[:]
        tempmaxdist = maxdist
        count = 0
        while len(todo) > 0:
            if len(todo) > 20000:
                paths = []
                todo = temptodo[:]
                tempmaxdist -= 10
                count += 1
            thepath = todo.pop()
            search = True
            while search:
                if len(thepath) > tempmaxdist + 20:
                    search = False
                if len(thepath) > nmersize + 20 and thepath[-nmersize - 20:] in endnmer:
                    paths.append(thepath)
                if search:
                    newn = thepath[-nmersize + 1:] + 'a'
                    if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        acount = self.nmerdict[newn]
                    else:
                        acount = 0
                    newn = thepath[-nmersize + 1:] + 't'
                    if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        tcount = self.nmerdict[newn]
                    else:
                        tcount = 0
                    newn = thepath[-nmersize + 1:] + 'c'
                    if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        ccount = self.nmerdict[newn]
                    else:
                        ccount = 0
                    newn = thepath[-nmersize + 1:] + 'g'
                    if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        gcount = self.nmerdict[newn]
                    else:
                        gcount = 0
                    if max([acount, tcount, ccount, gcount]) < nmercut:
                        search = False
                    else:
                        maxnucl = max([(acount, 'a'), (tcount, 't'), (ccount, 'c'), (gcount, 'g')])[1]
                        bases = 'atcg'.replace(maxnucl, '')
                        for i in bases:
                            newn = thepath[-nmersize + 1:] + i
                            if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                                newn = newn[::-1]
                                newn = newn.translate(transtab)
                            if newn in self.nmerdict and self.nmerdict[newn] > nmerave:
                                todo.append(thepath + i)
                        thepath += maxnucl
        countdict = {}
        for i in paths:
            if i[-nmersize - 20:] in countdict:
                acount = 0
                total = 0
                for j in range(0, len(i) - nmersize):
                    newn = i[j:j+nmersize]
                    if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    acount += 1
                    try:
                        total += self.nmerdict[newn]
                    except:
                        pass
                if total * 1.0 / acount > countdict[i[-nmersize - 20:]][1]:
                    countdict[i[-nmersize - 20:]] = (i, total * 1.0 / acount)
            else:
                acount = 0
                total = 0
                for j in range(0, len(i) - nmersize):
                    newn = i[j:j+nmersize]
                    if newn[nmersize/2] == 'a' or newn[nmersize/2] == 'c':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    acount += 1
                    try:
                        total += self.nmerdict[newn]
                    except:
                        pass
                countdict[i[-nmersize - 20:]] = (i, total * 1.0 / acount)
        paths = []
        for i in countdict:
            paths.append(countdict[i][0])
        return paths

    def getPerms(self, nmer):
        outlist = set()
        for i in range(len(nmer)):
            outlist.add(nmer[:i] + 'a' + nmer[i+1:])
            outlist.add(nmer[:i] + 't' + nmer[i+1:])
            outlist.add(nmer[:i] + 'c' + nmer[i+1:])
            outlist.add(nmer[:i] + 'g' + nmer[i+1:])
        return(list(outlist))

    def get_db_edges(self):
        nmersize = self.nmersize.get()
        ends = {}
        for i in self.contigDict:
            name = i
            seq = self.contigDict[i].forseq
            nmer = seq[:nmersize + 20]
            nmerlist = self.getPerms(nmer)
            for nmer in nmerlist:
                if nmer in ends:
                    ends[nmer] += ((name, True),)
                else:
                    ends[nmer] = ((name, True),)
            ends[nmer] = ((name, True),)
            nmer = seq[-nmersize - 20:]
            nmer = nmer[::-1]
            nmer = nmer.translate(transtab)
            nmerlist = self.getPerms(nmer)
            for nmer in nmerlist:
                if nmer in ends:
                    ends[nmer] += ((name, False),)
                else:
                    ends[nmer] = ((name, False),)
        for i in self.contigDict:
            forpaths = self.dbpath(self.contigDict[i].forseq[-nmersize:], ends)
            for j in forpaths:
                nmer = j[-nmersize - 20:]
                for k in ends[nmer]:
                    if not (j in self.contigDict[k[0]].forseq or j in self.contigDict[k[0]].revseq):
                        if len(j) < nmersize * 2 + 20:
                            appendpath = nmersize * 2 - len(j) + 20
                        else:
                            appendpath = j[nmersize:-nmersize-20]
                        if k[1]:
                            self.edgelist.append((i, True, k[0], True, appendpath))
                        else:
                            self.edgelist.append((i, True, k[0], False, appendpath))
            revpaths = self.dbpath(self.contigDict[i].revseq[-nmersize:], ends)
            for j in revpaths:
                nmer = j[-nmersize - 20:]
                for k in ends[nmer]:
                    if not (j in self.contigDict[k[0]].forseq or j in self.contigDict[k[0]].revseq):
                        if len(j) < nmersize * 2 + 20:
                            appendpath = nmersize * 2 - len(j) + 20
                        else:
                            appendpath = j[nmersize:-nmersize-20]
                        if k[1]:
                            self.edgelist.append((i, False, k[0], True, appendpath))
                        else:
                            self.edgelist.append((i, False, k[0], False, appendpath))    
        del self.nmerdict

    def getflag(self, n, count=11):
        return "".join([str((int(n) >> y) & 1) for y in range(count-1, -1, -1)])
    
    def get_paired_edge(self):
        minoverlap = self.minoverlap.get()
        tempread = open(self.readfile.get())
        firstline = tempread.readline()
        if firstline[0] == '@':
            readType = 'fq'
        elif firstline[0] == '>':
            readType = 'fa'
        else:
            self.edgeconsole.append('Read file not valid type, skipping paired edge creation.')
            tempread.close()
            return
        tempread.close()
        insertsize, minpairedge, minrl = self.insertsize.get(), self.minpairedge.get(), self.minrl.get()
        thestderr = open(self.workingDir.get() + '/bwaerr.txt', 'wa')
        subprocess.Popen('bowtie2-build -f ' + self.workingDir.get() + '/contigs.fa ' + self.workingDir.get() + '/fullcontigs', shell=True, stderr=thestderr, stdout=thestderr).wait()
        if readType == 'fa':
            subprocess.Popen('bowtie2 --local -x ' + self.workingDir.get() + '/fullcontigs -U ' + self.readfile.get() + ' -S ' + self.workingDir.get() + '/aln2.sam -f -a', shell=True, stderr=thestderr, stdout=thestderr).wait()
        elif readType == 'fq':
            subprocess.Popen('bowtie2 --local -x ' + self.workingDir.get() + '/fullcontigs -U ' + self.readfile.get() + ' -S ' + self.workingDir.get() + '/aln2.sam -q -a', shell=True, stderr=thestderr, stdout=thestderr).wait()
        thestderr.close()
        samfile = open(self.workingDir.get() + '/aln2.sam')
        zereadpair = False
        count = 0
        edgedict = {}
        for line in samfile:
            if not line.startswith('@'):
                count += 1
                query, flag, ref, pos, mapq, cigar = line.split()[:6]
                theflag = self.getflag(flag)
                if theflag[-9] == '0':
                    zereadpair = not zereadpair
                    if zereadpair:
                        firststart = []
                        firstend = []
                pos = int(pos)
                num = ''
                length = 0
                if cigar != '*':
                    for q in cigar:
                        if not q.isdigit():
                            if q != 'I' and q != 'H' and q != 'S':
                                length += int(num)
                            num = ''
                        else:
                            num += q
                if zereadpair:
                    if pos > 0 and length >= minrl:
                        if pos + length < insertsize and theflag[-5] != '0':
                            firststart.append((ref, pos + length))
                        if pos > self.contigDict[ref].length - insertsize and theflag[-5] == '0':
                            firstend.append((ref, self.contigDict[ref].length - pos))
                else:
                    secstart = False
                    secend = False
                    if pos > 0 and length >= minrl:
                        if pos + length < insertsize and theflag[-5] != '0':
                            j = (ref, pos + length)
                            secstart = True
                        if pos > self.contigDict[ref].length - insertsize and theflag[-5] == '0':
                            j = (ref, self.contigDict[ref].length - pos)
                            secend = True
                        if secstart:
                            for i in firststart:
                                if i[0] < j[0]:
                                    if '-' + i[0] in edgedict:
                                        if '-' + j[0] in edgedict['-' + i[0]]:
                                            edgedict['-' + i[0]]['-' + j[0]].append(i[1] + j[1])
                                        else:
                                            edgedict['-' + i[0]]['-' + j[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict['-' + i[0]] = {}
                                        edgedict['-' + i[0]]['-' + j[0]] = [i[1] + j[1]]
                                else:
                                    if '-' + j[0] in edgedict:
                                        if '-' + i[0] in edgedict['-' + j[0]]:
                                            edgedict['-' + j[0]]['-' + i[0]].append(i[1] + j[1])
                                        else:
                                            edgedict['-' + j[0]]['-' + i[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict['-' + j[0]] = {}
                                        edgedict['-' + j[0]]['-' + i[0]] = [i[1] + j[1]]
                            for i in firstend:
                                if i[0] < j[0]:
                                    if i[0] in edgedict:
                                        if '-' + j[0] in edgedict[i[0]]:
                                            edgedict[i[0]]['-' + j[0]].append(i[1] + j[1])
                                        else:
                                            edgedict[i[0]]['-' + j[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict[i[0]] = {}
                                        edgedict[i[0]]['-' + j[0]] = [i[1] + j[1]]
                                else:
                                    if '-' + j[0] in edgedict:
                                        if i[0] in edgedict['-' + j[0]]:
                                            edgedict['-' + j[0]][i[0]].append(i[1] + j[1])
                                        else:
                                            edgedict['-' + j[0]][i[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict['-' + j[0]] = {}
                                        edgedict['-' + j[0]][i[0]] = [i[1] + j[1]]
                        elif secend:
                            for i in firststart:
                                if i[0] < j[0]:
                                    if '-' + i[0] in edgedict:
                                        if j[0] in edgedict['-' + i[0]]:
                                            edgedict['-' + i[0]][j[0]].append(i[1] + j[1])
                                        else:
                                            edgedict['-' + i[0]][j[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict['-' + i[0]] = {}
                                        edgedict['-' + i[0]][j[0]] = [i[1] + j[1]]
                                else:
                                    if j[0] in edgedict:
                                        if '-' + i[0] in edgedict[j[0]]:
                                            edgedict[j[0]]['-' + i[0]].append(i[1] + j[1])
                                        else:
                                            edgedict[j[0]]['-' + i[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict[j[0]] = {}
                                        edgedict[j[0]]['-' + i[0]] = [i[1] + j[1]]
                            for i in firstend:
                                if i[0] < j[0]:
                                    if i[0] in edgedict:
                                        if j[0] in edgedict[i[0]]:
                                            edgedict[i[0]][j[0]].append(i[1] + j[1])
                                        else:
                                            edgedict[i[0]][j[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict[i[0]] = {}
                                        edgedict[i[0]][j[0]] = [i[1] + j[1]]
                                else:
                                    if j[0] in edgedict:
                                        if i[0] in edgedict[j[0]]:
                                            edgedict[j[0]][i[0]].append(i[1] + j[1])
                                        else:
                                            edgedict[j[0]][i[0]] = [i[1] + j[1]]
                                    else:
                                        edgedict[j[0]] = {}
                                        edgedict[j[0]][i[0]] = [i[1] + j[1]]
        newedges = []
        for i in edgedict:
            for j in edgedict[i]:
                if len(edgedict[i][j]) >= minpairedge:
                    if i[0] == '-':
                        startit = (i[1:], False)
                    else:
                        startit = (i, True)
                    if j[0] == '-':
                        endit = (j[1:], True)
                    else:
                        endit = (j, False)
                    gotit1, gotit2, gotit = True, True, True
                    for qq in self.edgelist:
                        if (qq[0] == startit[0] and qq[1] == startit[1]) or (qq[2] == startit[0] and qq[3] != startit[1]):
                            gotit1 = False
                        if (qq[2] == endit[0] and qq[3] == endit[1]) or (qq[0] == endit[0] and qq[1] != endit[1]):
                            gotit2 = False
                        if qq[0] == startit[0] and qq[1] == startit[1] and qq[2] == endit[0] and qq[3] == endit[1]:
                            gotit = False
                        if qq[0] == endit and qq[1] != endit[1] and qq[2] == startit[0] and qq[3] != startit[1]:
                            gotit = False
                    overlap = 'nnnnnnnnn'
                    if startit[1]:
                        seqa = self.contigDict[startit[0]].forseq
                    else:
                        seqa = self.contigDict[startit[0]].revseq
                    if endit[1]:
                        seqb = self.contigDict[endit[0]].forseq
                    else:
                        seqb = self.contigDict[endit[0]].revseq
                    if sum(edgedict[i][j])/len(edgedict[i][j]) < insertsize:
                        for q in range(5, minoverlap):
                            if seqa[-q:] == seqb[:q]:
                                overlap = q
                    if overlap == 'nnnnnnnnn':
                        if gotit1 or gotit2:
                            newedges.append((startit[0], startit[1], endit[0], endit[1], overlap))
                    else:
                        if gotit:
                            newedges.append((startit[0], startit[1], endit[0], endit[1], overlap))
        self.edgelist = self.edgelist + newedges

    def get_long_edge(self):
        readfile = open(self.readfile.get())
        longestread = 0
        first = True
        out = open(self.workingDir.get() + '/reads.fa', 'w')
        readlendict = {}
        for line in readfile:
            if line.startswith('>'):
                if first:
                    first = False
                    getfa = True
                else:
                    readlendict[name] = len(seq)
                    if len(seq) > longestread:
                        longestread = len(seq)
                name = line.rstrip()[1:]
                seq = ''
                out.write(line)
            elif getfa:
                seq += line.rstrip()
                out.write(line)
        out.close()
        readfile.close()
        readlendict[name] = len(seq)
        reffile = open(self.workingDir.get() + '/longref.fa', 'w')
        for i in self.contigDict:
            reffile.write('>s' + i + '\n' + self.contigDict[i].forseq + '\n')
            reffile.write('>e' + i + '\n' + self.contigDict[i].revseq + '\n')
        reffile.close()
        subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/longref -in ' + self.workingDir.get() + '/longref.fa', shell=True, stdout=subprocess.PIPE).wait()
        subprocess.Popen('blastn -db ' + self.workingDir.get() + '/longref -outfmt 6 -num_threads 8 -query ' + self.workingDir.get() + '/reads.fa -out ' + self.workingDir.get() + '/long.out', shell=True).wait()
        blastfile = open(self.workingDir.get() + '/long.out')
        edgedict = {}
        lastq = ''
        aroundzero = 20
        for line in blastfile:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
            ident = float(ident)
            if query != lastq:
                if lastq != '':
                    for i in readstart:
                        for j in readend:
                            if i < j:
                                if i in edgedict:
                                    if j in edgedict[i]:
                                        edgedict[i][j] += 1
                                    else:
                                        edgedict[i][j] = 1
                                else:
                                    edgedict[i] = {j:1}
                            else:
                                if j in edgedict:
                                    if i in edgedict[j]:
                                        edgedict[j][i] += 1
                                    else:
                                        edgedict[j][i] = 1
                                else:
                                    edgedict[j] = {i:1}
                readstart = []
                readend = []
                lastq = query
            if ident > self.minlongident.get() and length > self.minlongover.get():
                if qstart <= aroundzero and rstop <= aroundzero:
                    readstart.append(subject)
                if qstop  >= readlendict[query] - aroundzero and rstart <= aroundzero:
                    readend.append(subject)
        for i in edgedict:
            for j in edgedict[i]:
                if edgedict[i][j] >= self.minlong:
                    if i[0] == 's':
                        if j[0] == 's':
                            self.edgelist.append((i[1:], False, j[1:], True, 'nnnnnnnnn'))
                        else:
                            self.edgelist.append((i[1:], False, j[1:], False, 'nnnnnnnnn'))
                    else:
                        if j[0] == 's':
                            self.edgelist.append((i[1:], True, j[1:], True, 'nnnnnnnn'))
                        else:
                            self.edgelist.append((i[1:], True, j[1:], False, 'nnnnnnnnn'))
            

    def get_nmer_freq(self):
        nmersize, nmercut, reads  = self.nmersize.get(), self.nmercut.get(), self.readfile.get()
        nucl = set('atcg')
        self.nmerdict = {}
        if reads[-3:] == '.gz':
            readfile = gzip.open(reads)
        else:
            readfile = open(reads)
        getseq = False
        getfaseq = False
        for line in readfile:
            if line.startswith('@'):
                header = line[1:].rstrip()
                getseq = True
            elif line.startswith('>'):
                if getfaseq:
                    for i in range(0, len(seq) - nmersize):
                        nmer = seq[i:i+nmersize]
                        if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                            nmer = nmer[::-1]
                            nmer = nmer.translate(transtab)
                        if nmer in self.nmerdict:
                            self.nmerdict[nmer] += 1
                        else:
                            if set(nmer).issubset(nucl):
                                self.nmerdict[nmer] = 1
                getfaseq = True
                header = line[1:].rstrip()
                seq = ''
            elif getfaseq:
                seq += line.rstrip().lower()
            elif getseq:
                getseq = False
                seq = line.rstrip().lower()
                for i in range(0, len(seq) - nmersize + 1):
                    nmer = seq[i:i+nmersize]
                    if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                        nmer = nmer[::-1]
                        nmer = nmer.translate(transtab)
                    if nmer in self.nmerdict:
                        self.nmerdict[nmer] += 1
                    else:
                        if set(nmer).issubset(nucl):
                            self.nmerdict[nmer] = 1
        if getfaseq:
            for i in range(0, len(seq) - nmersize + 1):
                nmer = seq[i:i+nmersize]
                if nmer[nmersize/2] == 'a' or nmer[nmersize/2] == 'c':
                    nmer = nmer[::-1]
                    nmer = nmer.translate(transtab)
                if nmer in self.nmerdict:
                    self.nmerdict[nmer] += 1
                else:
                    if set(nmer).issubset(nucl):
                        self.nmerdict[nmer] = 1
        out = open(self.workingDir.get() + '/nmers', 'w')
        for i in self.nmerdict:
            out.write(i + '\t' + str(self.nmerdict[i]) + '\n')
        out.close()
        out = open(self.workingDir.get() + '/nmerfreq.csv', 'w')
        freqDict = {}
        for i in self.nmerdict:
            if self.nmerdict[i] in freqDict:
                freqDict[self.nmerdict[i]] += 1
            else:
                freqDict[self.nmerdict[i]] = 1
        for i in range(1, max(freqDict) + 1):
            if i in freqDict:
                out.write(str(i) + '\t' + str(freqDict[i]) + '\n')
            else:
                out.write(str(i) + '\t0\n')
        out.close()
        self.nmerfile = self.workingDir.get() + '/nmers'

    def read_nmer_freq(self):
        self.nmerdict = {}
        nmercut = self.nmercut.get()
        nmers = open(self.nmerfile)
        for line in nmers:
            nmer, freq = line.split()
            freq = int(freq)
            if freq > nmercut:
                self.nmerdict[nmer] = freq

    def getnmercut(self):
        nmerfile = open(self.nmerfile)
        freqdict = {}
        for line in nmerfile:
            freq = int(line.split()[1])
            if freq in freqdict:
                freqdict[freq] += 1
            else:
                freqdict[freq] = 1
        nmerfile.close()
        maxfreq = 0
        trfreq = float('inf')
        mf = None
        tf = None
        rising = True
        for freq in range(0, max(freqdict)):
            try:
                if freqdict[freq] > maxfreq and not rising:
                    mf = freq
                    maxfreq = freqdict[freq]
                if rising:
                    if freqdict[freq] < trfreq:
                        tf = freq
                        trfreq = freqdict[freq]
                    else:
                        rising = False
            except:
                pass
        nmerave = int(mf)/2
        nmercut = tf

    def untangleEdges(self):
        edgeDict = {}
        for i in self.edgelist:
            if not i[0] in edgeDict:
                edgeDict[i[0]] = ([], [])
            if not i[2] in edgeDict:
                edgeDict[i[2]] = ([], [])
            if i[1]:
                edgeDict[i[0]][0].append((i[2], i[3], i[4]))
            else:
                edgeDict[i[0]][1].append((i[2], i[3], i[4]))
            if type(i[4]) == int:
                overlap = i[4]
            else:
                overlap = i[4][::-1].translate(transtab)
            if i[3]:
                edgeDict[i[2]][1].append((i[0], not i[1], overlap))
            else:
                edgeDict[i[2]][0].append((i[0], not i[1], overlap))
        paths = []
        tosearch = set()
        for i in self.edgelist:
            tosearch.add(i[:2])
        for i in range(3, 5):
            for j in tosearch:
                todo = [[j]]
                while len(todo) > 0:
                    currpath = todo.pop()
                    if len(currpath) >= i:
                        paths.append(currpath)
                    else:
                        if currpath[-1][1]:
                            for k in edgeDict[currpath[-1][0]][0]:
                                todo.append(currpath + [k])
                        else:
                            for k in edgeDict[currpath[-1][0]][1]:
                                todo.append(currpath + [k])                         
        pathseqlist = []
        temppaths = []
        for i in paths:
            if i[0][1]:
                pathseq = self.contigDict[i[0][0]].forseq
            else:
                pathseq = self.contigDict[i[0][0]].revseq
            for q in range(1, len(i)):
                overlap = i[q][2]
                if i[q][1]:
                    if type(overlap) == int:
                        pathseq += self.contigDict[i[q][0]].forseq[overlap:]
                    else:
                        pathseq += overlap
                        pathseq += self.contigDict[i[q][0]].forseq
                else:
                    if type(overlap) == int:
                        pathseq += self.contigDict[i[q][0]].revseq[overlap:]
                    else:
                        pathseq += overlap
                        pathseq += self.contigDict[i[q][0]].revseq
            if - self.contigDict[i[-1][0]].length + 100 < 0:
                pathseq = pathseq[max([0, self.contigDict[i[0][0]].length - 100]): - self.contigDict[i[-1][0]].length + 100]
            else:
                pathseq = pathseq[max([0, self.contigDict[i[0][0]].length - 100]):]
            if len(pathseq) <= self.maxdist.get() + 200:
                temppaths.append(i)
                pathseqlist.append(pathseq)
        paths = temppaths
        edgeseqlist = []
        for i in self.edgelist:
            if i[1]:
                pathseq = self.contigDict[i[0]].forseq
            else:
                pathseq = self.contigDict[i[0]].revseq
            overlap = i[4]
            if i[3]:
                if type(overlap) == int:
                    pathseq += self.contigDict[i[2]].forseq[overlap:]
                else:
                    pathseq += overlap
                    pathseq += self.contigDict[i[2]].forseq
            else:
                if type(overlap) == int:
                    pathseq += self.contigDict[i[2]].revseq[overlap:]
                else:
                    pathseq += overlap
                    pathseq += self.contigDict[i[2]].revseq
            if - self.contigDict[i[2]].length + 100 < 0:
                pathseq = pathseq[max([0, self.contigDict[i[0]].length - 100]): - self.contigDict[i[2]].length + 100]
            else:
                pathseq = pathseq[max([0, self.contigDict[i[0]].length - 100]):]
            edgeseqlist.append(pathseq)
        tempout = open(self.workingDir.get() + '/multiedge.fa', 'w')
        count = 0
        for i in pathseqlist:
            if len(i) > 100:
                tempout.write('>' + str(count) + '\n' + i + '\n')
            count += 1
        tempout.close()
        tempout = open(self.workingDir.get() + '/singedge.fa', 'w')
        count = 0
        for i in edgeseqlist:
            if len(i) >= 100:
                tempout.write('>' + str(count) + '\n' + i + '\n')
            count += 1
        tempout.close()
        subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/multidb -in ' + self.workingDir.get() + '/multiedge.fa', shell=True, stdout=subprocess.PIPE).wait()
        subprocess.Popen('blastn -db ' + self.workingDir.get() + '/multidb -outfmt 6 -num_threads 8 -query ' + self.workingDir.get() + '/singedge.fa -out ' + self.workingDir.get() + '/singlemult.out', shell=True).wait()
        singlemalt = open(self.workingDir.get() + '/singlemult.out')
        minident = 99
        minlength = 0.99
        vetoset = set()
        for line in singlemalt:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
            ident = float(ident)
            if length > len(edgeseqlist[int(query)]) * minlength and length > len(pathseqlist[int(subject)]) * minlength and ident > minident \
              and self.edgelist[int(query)][:2] == paths[int(subject)][0][:2] and self.edgelist[int(query)][2:4] == paths[int(subject)][-1][:2]:
                vetoset.add(int(query))
        outlist = []
        for i in range(len(self.edgelist)):
            if not i in vetoset:
                outlist.append(self.edgelist[i])
        self.edgelist = outlist


    def removeDuplicates(self):
        newedgelist = {}
        edgeCountDict = {}
        for i in self.edgelist:
            if i[:4] in edgeCountDict:
                edgeCountDict[i[:4]] += 1
            else:
                edgeCountDict[i[:4]] = 1
            if i[:4] in newedgelist:
                if type(i[4]) == int:
                    lena = i[4]
                else:
                    if i[4] == 'nnnnnnnnn':
                        lena = - float('inf')
                    else:
                        lena = - len(i[4])
                if type(newedgelist[i[:4]]) == int:
                    lenb = newedgelist[i[:4]]
                else:
                    if newedgelist[i[:4]] == 'nnnnnnnnn':
                        lenb = - float('inf')
                    else:
                        lenb = - len(newedgelist[i[:4]])
                if lena > lenb:
                    newedgelist[i[:4]] = i[4]
            elif (i[2], not i[3], i[0], not i[1]) in newedgelist:
                if type(i[4]) == int:
                    lena = i[4]
                else:
                    if i[4] == 'nnnnnnnnn':
                        lena = - float('inf')
                    else:
                        lena = - len(i[4])
                if type(newedgelist[(i[2], not i[3], i[0], not i[1])]) == int:
                    lenb = newedgelist[(i[2], not i[3], i[0], not i[1])]
                else:
                    if newedgelist[(i[2], not i[3], i[0], not i[1])] == 'nnnnnnnnn':
                        lenb = - float('inf')
                    else:
                        lenb = - len(newedgelist[(i[2], not i[3], i[0], not i[1])])
                if lena > lenb:
                    newedgelist[(i[2], not i[3], i[0], not i[1])] = i[4]            
            else:
                newedgelist[i[:4]] = i[4]
        self.edgelist = []
        for i in newedgelist:
            self.edgelist.append(i + (newedgelist[i],))

    def writeCSAG(self):
        out = open(self.workingDir.get() + '/CSAG.txt', 'w')
        for i in self.contigDict:
            out.write('NODE\t' + i + '\t' + self.contigDict[i].forseq + '\n')
        for i in self.edgelist:
            out.write('EDGE\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\n')
        out.close()

    def create_comp(self):
        try:
            self.blast_options.destroy()
        except:
            pass
        self.blast_options = Toplevel()
        self.frame2 = Frame(self.blast_options)
        self.blast_options.geometry('+20+30')
        self.blast_options.title('Create comparison')
        self.reffilelabel = Label(self.frame2, text='Reference file')
        self.reffilelabel.grid(column=0, row=1)
        self.reffileentry = Entry(self.frame2, textvariable=self.reffile)
        self.reffileentry.grid(column=1, row=1)
        self.reffilebutton = Button(self.frame2, text='...', command=self.loadref)
        self.reffilebutton.grid(column=2, row=1)
        self.blastfilelabel = Label(self.frame2, text='Blast file')
        self.blastfilelabel.grid(column=0, row=2)
        self.blastfileentry = Entry(self.frame2, textvariable=self.blastfile)
        self.blastfileentry.grid(column=1, row=2)
        self.blastfilebutton = Button(self.frame2, text='...', command=self.loadblast)
        self.blastfilebutton.grid(column=2, row=2)
        self.minlengthblastlabel = Label(self.frame2, text='Minimum length (bp)')
        self.minlengthblastlabel.grid(column=0, row=3)
        self.minlengthblastentry = Entry(self.frame2, textvariable=self.minlengthblast)
        self.minlengthblastentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.minlengthratiolabel = Label(self.frame2, text='Minimum length (ratio)')
        self.minlengthratiolabel.grid(column=0, row=4)
        self.minlengthratioentry = Entry(self.frame2, textvariable=self.minlengthratio)
        self.minlengthratioentry.grid(column=1, row=4, columnspan=2, sticky=EW)
        self.minidentlabel = Label(self.frame2, text='Minimum identity')
        self.minidentlabel.grid(column=0, row=5)
        self.minidententry = Entry(self.frame2, textvariable=self.minident)
        self.minidententry.grid(column=1, row=5, columnspan=2, sticky=EW)
        self.minbitscorelabel = Label(self.frame2, text='Minimum bitscore')
        self.minbitscorelabel.grid(column=0, row=6)
        self.minbitscoreentry = Entry(self.frame2, textvariable=self.minbitscore)
        self.minbitscoreentry.grid(column=1, row=6, columnspan=2, sticky=EW)
        self.maxevaluelabel = Label(self.frame2, text='Minimum evalue')
        self.maxevaluelabel.grid(column=0, row=7)
        self.maxevalueentry = Entry(self.frame2, textvariable=self.maxevalue)
        self.maxevalueentry.grid(column=1, row=7, columnspan=2, sticky=EW)
        self.closeblastbut = Button(self.frame2, text='Ok', command=self.okblast)
        self.closeblastbut.grid(column=2, row=8, sticky=E, pady=5)
        self.frame2.grid(padx=20, pady=20)

    def okblast(self):
        try:
            self.blast_options.destroy()
        except:
            pass
        ref = open(self.reffile.get())
        first = True
        self.refDict = {}
        self.reforder = []
        for line in ref:
            if line.startswith('>'):
                if first:
                    first = False
                else:
                    self.refDict[name] = seq
                    self.reforder.append(name)
                name = line.split()[0][1:]
                seq = ''
            else:
                seq += line.rstrip()
        self.refDict[name] = seq
        self.reforder.append(name)
        if self.blastfile.get() == '' and self.contigfile.get() != '' and self.reffile.get() != '':
            blastit = tkMessageBox.askquestion('No blast files.', 'Create BLAST output?')
            if blastit == 'yes':
                subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/tempdb -in ' + self.reffile.get(), stdout=subprocess.PIPE, shell=True).wait()
                subprocess.Popen('blastn -task blastn -db ' + self.workingDir.get() + '/tempdb -outfmt 6 -query ' + self.workingDir.get() + '/contigs.fa -out ' + self.workingDir.get() + '/query_tempdb.out', shell=True).wait()
                self.blastfile.set(self.workingDir.get() + '/query_tempdb.out')
                self.update_console('BLAST file created.')
        if self.blastfile.get() != '':
            blastfile = open(self.blastfile.get())
            self.hitlist = []
            for line in blastfile:
                query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
                qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
                ident, eval, bitscore = map(float, [ident, eval, bitscore])
                if ident >= self.minident.get() and length >= self.minlengthblast.get() and mm <= self.maxmm.get() and eval <= self.maxevalue.get() and bitscore >= self.minbitscore.get():
                    self.hitlist.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
        self.update_console('Comparison created.')


    def writeWorkCont(self):
        contout = open(self.workingDir.get() + '/contigs.fa', 'w')
        for i in self.contigDict:
            contout.write('>' + i + '\n')
            for j in range(0, len(self.contigDict[i].forseq), 60):
                contout.write(self.contigDict[i].forseq[j:j+60] + '\n')
        contout.close()

    def self_compare(self):
        try:
            self.blast_options.destroy()
        except:
            pass
        self.blast_options = Toplevel()
        self.frame2 = Frame(self.blast_options)
        self.blast_options.geometry('+20+30')
        self.blast_options.title('Self comparison')
        self.reffilelabel = Label(self.frame2, text='Comparison file')
        self.reffilelabel.grid(column=0, row=1)
        self.reffileentry = Entry(self.frame2, textvariable=self.selfcomparefile)
        self.reffileentry.grid(column=1, row=1)
        self.reffilebutton = Button(self.frame2, text='...', command=self.loadselfcomp)
        self.reffilebutton.grid(column=2, row=1)
        self.intralabel = Label(self.frame2, text='Show intra contig hits')
        self.intralabel.grid(column=0, row=2)
        self.intraentry = Checkbutton(self.frame2, variable=self.intra)
        self.intraentry.grid(column=1, row=2, columnspan=2, sticky=EW)
        self.minlengthblastlabel = Label(self.frame2, text='Minimum length (bp)')
        self.minlengthblastlabel.grid(column=0, row=3)
        self.minlengthblastentry = Entry(self.frame2, textvariable=self.minlengthblast)
        self.minlengthblastentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.minlengthratiolabel = Label(self.frame2, text='Minimum length (ratio)')
        self.minlengthratiolabel.grid(column=0, row=4)
        self.minlengthratioentry = Entry(self.frame2, textvariable=self.minlengthratio)
        self.minlengthratioentry.grid(column=1, row=4, columnspan=2, sticky=EW)
        self.minidentlabel = Label(self.frame2, text='Minimum identity')
        self.minidentlabel.grid(column=0, row=5)
        self.minidententry = Entry(self.frame2, textvariable=self.minident)
        self.minidententry.grid(column=1, row=5, columnspan=2, sticky=EW)
        self.minbitscorelabel = Label(self.frame2, text='Minimum bitscore')
        self.minbitscorelabel.grid(column=0, row=6)
        self.minbitscoreentry = Entry(self.frame2, textvariable=self.minbitscore)
        self.minbitscoreentry.grid(column=1, row=6, columnspan=2, sticky=EW)
        self.maxevaluelabel = Label(self.frame2, text='Minimum evalue')
        self.maxevaluelabel.grid(column=0, row=7)
        self.maxevalueentry = Entry(self.frame2, textvariable=self.maxevalue)
        self.maxevalueentry.grid(column=1, row=7, columnspan=2, sticky=EW)
        self.closeblastbut = Button(self.frame2, text='Ok', command=self.ok_self_compare)
        self.closeblastbut.grid(column=2, row=8, sticky=E, pady=5)
        self.frame2.grid(padx=20, pady=20)


    def ok_self_compare(self):
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        self.blast_options.destroy()
        self.queue = Queue.Queue()
        self.thethread = threading.Thread(target=self.self_hits_thread)
        self.thethread.start()
        self.update_self_hits()


    def update_self_hits(self):
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text != 'Hits found.':
                    root.after(1000, self.update_self_hits)
                else:
                    self.draw_self_hits()
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            self.update_console('Self comparison failed, please check console output.')
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_self_hits)
            return

    def self_hits_thread(self):
        if self.selfcomparefile.get() == '':
            stderr = open(self.workingDir.get() + '/bwaerr.txt', 'wa')
            subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/contigdb -in ' +
                             self.workingDir.get() + '/contigs.fa', shell=True, stdout=stderr).wait()
            subprocess.Popen('blastn -db ' + self.workingDir.get() + '/contigdb -outfmt 6 -num_threads 8 -query ' +
                             self.workingDir.get() + '/contigs.fa -out ' + self.workingDir.get() + '/contigscontigs.out', shell=True).wait()
            stderr.close()
            blastfile = open(self.workingDir.get() + '/contigscontigs.out')
        else:
            blastfile = open(self.selfcomparefile.get())
        self.queue.put('Comparison created.')
        self.selfhit = []
        for line in blastfile:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
            ident, eval, bitscore = map(float, [ident, eval, bitscore])
            if ident >= self.minident.get() and length >= self.minlengthblast.get() and mm <= self.maxmm.get() and eval <= self.maxevalue.get() and bitscore >= self.minbitscore.get()\
              and (length >= self.minlengthratio.get() * self.contigDict[query].length or length >= self.minlengthratio.get() * self.contigDict[subject].length):
                if query != subject:
                    if query < subject:
                        self.selfhit.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
                elif self.intra.get() == 1:
                    if qstart < min([rstart, rstop]):
                        self.selfhit.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
        blastfile.close()
        self.queue.put('Hits found.')


    def draw_self_hits(self):
        # need to add in hits for duplicates
        self.canvas.delete('sblast')
        hitnum = 0
        for i in self.selfhit:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = i
            if query in self.visible and subject in self.visible:
                hitnum += 1
                if rstart < rstop and self.contigDict[query].orient[0] == self.contigDict[subject].orient[0]:
                    colour = '#F85931'
                elif rstart < rstop:
                    colour = '#EDB92E'
                elif rstop < rstart and self.contigDict[query].orient[0] == self.contigDict[subject].orient[0]:
                    colour = '#EDB92E'
                else:
                    colour = '#F85931'
                qcoords = self.canvas.coords('c' + query)
                rcoords = self.canvas.coords('c' + subject)
                if self.contigDict[query].orient[0]:
                    self.canvas.create_rectangle(qcoords[0] + qstart / self.contigDict[query].scaledown, qcoords[1],
                                                 qcoords[0] + qstop / self.contigDict[query].scaledown, qcoords[3],
                                                  fill=colour, tags=('c' + query, 'sblast', 'map', 'c' + query + 'h', 'self' + str(hitnum)))
                    startx1 = qcoords[0] + qstart / self.contigDict[query].scaledown
                    startx2 = qcoords[0] + qstop / self.contigDict[query].scaledown
                else:
                    self.canvas.create_rectangle(qcoords[2] - qstart / self.contigDict[query].scaledown, qcoords[1],
                                                 qcoords[2] - qstop / self.contigDict[query].scaledown, qcoords[3],
                                                  fill=colour, tags=('c' + query, 'sblast', 'map', 'c' + query + 'h', 'self' + str(hitnum)))

                    startx1 = qcoords[2] - qstart / self.contigDict[query].scaledown
                    startx2 = qcoords[2] - qstop / self.contigDict[query].scaledown
                if self.contigDict[subject].orient[0]:
                    self.canvas.create_rectangle(rcoords[0] + rstart / self.contigDict[subject].scaledown, rcoords[1],
                                                 rcoords[0] + rstop / self.contigDict[subject].scaledown, rcoords[3],
                                                  fill=colour, tags=('c' + subject, 'sblast', 'map', 'c' + subject + 'h', 'self' + str(hitnum)))
                    endx1 = rcoords[0] + rstart / self.contigDict[subject].scaledown
                    endx2 = rcoords[0] + rstop / self.contigDict[subject].scaledown
                else:
                    self.canvas.create_rectangle(rcoords[2] - rstart / self.contigDict[subject].scaledown, rcoords[1],
                                                 rcoords[2] - rstop / self.contigDict[subject].scaledown, rcoords[3],
                                                  fill=colour, tags=('c' + subject, 'sblast', 'map', 'c' + subject + 'h', 'self' + str(hitnum)))
                    endx1 = rcoords[2] - rstart / self.contigDict[subject].scaledown
                    endx2 = rcoords[2] - rstop / self.contigDict[subject].scaledown
                starty = qcoords[1] + self.contigheight
                endy = rcoords[1] + self.contigheight
                self.canvas.create_line(startx1, starty, (startx1 + endx1) / 2, abs(startx1 - endx1) /4 + (starty + endy) / 2, endx1, endy,\
                                         smooth=True, width=1, tags=('c' + query + 's', 'c' + subject + 'e', 'arc'))
                self.canvas.create_line(startx2, starty, (startx2 + endx2) / 2, abs(startx2 - endx2) /4 + (starty + endy) / 2, endx2, endy,\
                                         smooth=True, width=1, tags=('c' + query + 's', 'c' + subject + 'e', 'arc'))
        self.canvas.tag_raise('arc')
        self.canvas.tag_raise('text')
        self.blast_options.destroy()

    def view_options(self):
        try:
            self.view_options.destroy()
        except:
            pass
        self.view_options = Toplevel()
        self.frame3 = Frame(self.view_options)
        self.view_options.geometry('+20+30')
        self.view_options.title('View contig graph')
        self.viewlabel = Label(self.frame3, text='Contigs to view:')
        self.viewlabel.grid(column=0, row=0)
        self.viewentry = OptionMenu(self.frame3, self.view, 'All', 'BLAST', 'List')#, 'Coverage', 'BLAST and Coverage', 'Blast or Coverage')
        self.viewentry.grid(column=1, row=0, columnspan=2, sticky=EW)
        self.viewlistlabel = Label(self.frame3, text='List of contigs:')
        self.viewlistlabel.grid(column=0, row=1)
        self.viewlistentry = Entry(self.frame3, textvariable=self.viewlist)
        self.viewlistentry.grid(column=1, row=1)
        self.viewlistbutton = Button(self.frame3, text='...', command=self.loadview)
        self.viewlistbutton.grid(column=2, row=1)
        self.viewreflabel = Label(self.frame3, text='View reference and BLAST hits:')
        self.viewreflabel.grid(column=0, row=2)
        self.viewrefentry = Checkbutton(self.frame3, variable=self.viewref)
        self.viewrefentry.grid(column=1, row=2, columnspan=2)
        self.shortenlabel = Label(self.frame3, text='Shorten long contigs:')
        self.shortenlabel.grid(column=0, row=3)
        self.shortenentry = OptionMenu(self.frame3, self.shorten, 'No', 'Log', 'Max', 'Min')
        self.shortenentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.scaledownlabel = Label(self.frame3, text='Scale: 1 pixel =')
        self.scaledownlabel.grid(column=0, row=4)
        self.scaledownenentry = Entry(self.frame3, textvariable=self.scaledown)
        self.scaledownenentry.grid(column=1, row=4)
        self.zebplabel = Label(self.frame3, text='bp')
        self.zebplabel.grid(column=2, row=4)
        self.okview = Button(self.frame3, text='Ok', command=self.ok_view)
        self.okview.grid(column=2, row=5, sticky=E)
        self.frame3.grid(column=0, row=0)

    def ok_view(self):
        self.canvas.delete(ALL)
        self.newscaledown = self.scaledown.get()
        self.visible = set()
        self.contigheight = 25
        self.currxscroll = self.originalxscroll
        self.curryscroll = self.originalyscroll
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        if self.shorten.get() == 'No':
            for i in self.contigDict:
                self.contigDict[i].xlength = self.contigDict[i].length / self.scaledown.get()
                self.contigDict[i].scaledown = self.scaledown.get()
        if self.shorten.get() == 'Log':
            for i in self.contigDict:
                self.contigDict[i].xlength = (2 ** math.log10(self.contigDict[i].length)) * (self.scaledown.get() / 4)
                self.contigDict[i].scaledown = self.contigDict[i].length * 1.0 / ((2 ** math.log10(self.contigDict[i].length)) * (self.scaledown.get() / 4))
        if self.shorten.get() == 'Max':
            for i in self.contigDict:
                self.contigDict[i].xlength = min([self.contigDict[i].length / self.scaledown.get(), 100])
                self.contigDict[i].scaledown = self.contigDict[i].length * 1.0 / min([self.contigDict[i].length / self.scaledown.get(), 100])
        if self.shorten.get() == 'Min':
            for i in self.contigDict:
                self.contigDict[i].xlength = max([self.contigDict[i].length / self.scaledown.get(), 100])
                self.contigDict[i].scaledown = self.contigDict[i].length * 1.0 / max([self.contigDict[i].length / self.scaledown.get(), 100])
        if self.view.get() == 'BLAST':
            if self.hitlist is None:
                tkMessageBox.showerror('Comparison not found.', 'Please perform comparison before choosing this method.')
                return
            elif self.hitlist == []:
                tkMessageBox.showerror('No hits found.', 'No blast hits were found, please perform ' +
                                                         'comparison again with different parameters.')
                self.view_options.destroy()
                return
            self.visible = self.getContigsWithHit()
            self.orderContigsBlast()
            self.drawContigs()
            self.drawEdges()
        elif self.view.get() == 'All':
            self.visible = set(self.contigDict)
            self.orderContigsAll()
            self.drawContigs()
            self.drawEdges()
        elif self.view.get() == 'List':
            self.orderContigsList()
            self.drawContigs()
            self.drawEdges()
        if self.viewref.get() and not self.hitlist is None:
            self.drawRefHits()
        self.canvas.move(ALL, -self.leftmost + 50, 0)
        self.currxscroll = self.rightmost - self.leftmost + 200
        self.curryscroll = 1000
        self.canvas.configure(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.view_options.destroy()

    def drawEdges(self):
        for i in self.edgelist:
            if i[0] in self.visible and i[2] in self.visible:
                if i[1]:
                    startx = self.contigDict[i[0]].xpos + self.contigDict[i[0]].xlength
                else:
                    startx = self.contigDict[i[0]].xpos
                starty = self.contigDict[i[0]].ypos + self.contigheight / 2
                if i[3]:
                    endx = self.contigDict[i[2]].xpos
                else:
                    endx = self.contigDict[i[2]].xpos + self.contigDict[i[2]].xlength
                endy = self.contigDict[i[2]].ypos + self.contigheight / 2
                self.canvas.create_line(startx, starty, (startx + endx) / 2, abs(startx - endx) /4 + (starty + endy) / 2, endx, endy,\
                                         smooth=True, width=3, tags=('c' + i[0] + 's', 'c' + i[2] + 'e', 'arc'))

    def getContigsWithHit(self):
        outset = set()
        for i in self.hitlist:
            outset.add(i[0])
        return outset
    
    def drawRefHits(self):
        self.hitlist.sort(key=lambda x: x[3], reverse=True)
        for i in self.hitlist:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = i
            if query in self.visible:
                if rstart < rstop:
                    colour = '#F85931'
                else:
                    colour = '#EDB92E'
                self.canvas.create_polygon(self.refpos[subject] + rstart / self.scaledown.get(), self.refline + self.contigheight,
                                           self.refpos[subject] + rstop / self.scaledown.get(), self.refline + self.contigheight,
                                            self.contigDict[query].xpos + qstop/self.contigDict[query].scaledown, self.contigline,
                                             self.contigDict[query].xpos + qstart/self.contigDict[query].scaledown, self.contigline,
                                             fill=colour, outline="black", tags=('c' + query + 'b', 'r' + subject + 'b', 'blast'))
        for i in self.reforder:
            if self.leftmost is None or self.refpos[i] < self.leftmost:
                self.leftmost = self.refpos[i]
            if self.rightmost is None or self.refpos[i] + len(self.refDict[i]) / self.scaledown.get() > self.rightmost:
                self.rightmost = self.refpos[i] + len(self.refDict[i]) / self.scaledown.get()
            self.canvas.create_rectangle(self.refpos[i], self.refline, \
                                         self.refpos[i] + len(self.refDict[i]) / self.scaledown.get(), self.refline + self.contigheight, \
                                         fill="#CE1836", tags=('r' + i, "ref", "map"))

    def drawContigs(self):
        for i in self.visible:
            if self.leftmost is None or self.contigDict[i].xpos < self.leftmost:
                self.leftmost = self.contigDict[i].xpos
            if self.rightmost is None or self.contigDict[i].xpos + self.contigDict[i].xlength > self.rightmost:
                self.rightmost = self.contigDict[i].xpos + self.contigDict[i].xlength
            self.canvas.create_rectangle(self.contigDict[i].xpos, self.contigDict[i].ypos,
                                         self.contigDict[i].xpos + self.contigDict[i].xlength, self.contigDict[i].ypos + self.contigheight,\
                                          fill="#009989", tags=('c' + i, 'contig', 'map'))
            if self.contigDict[i].orient[0]:
                dir = '+'
            else:
                dir = '-'
            thetext = i + ' ' + dir + ' ' + self.contigDict[i].strlen
            text = self.canvas.create_text(self.contigDict[i].xpos + 2, self.contigDict[i].ypos + self.contigheight/2,\
                                            fill='white', font=self.customFont, anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
            if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                self.canvas.delete(text)
                thetext = i + ' ' + dir
                text = self.canvas.create_text(self.contigDict[i].xpos + 2, self.contigDict[i].ypos + self.contigheight/2,\
                                                fill='white', font=self.customFont, anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                    self.canvas.delete(text)
                    thetext = i
                    text = self.canvas.create_text(self.contigDict[i].xpos + 2, self.contigDict[i].ypos + self.contigheight/2,\
                                                    fill='white', font=self.customFont, anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                    if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                        self.canvas.delete(text)

    def orderContigsList(self):
        listfile = open(self.viewlist.get())
        currx = 10
        for line in listfile:
            name = line.rstrip()
            self.visible.add(name)
            self.contigDict[name].xpos = currx
            self.contigDict[name].ypos = self.contigline
            curr += self.contigDict[name].xlength + 10
        listfile.close()

    def orderContigsAll(self):
        todo = set(self.contigDict)
        listoflists = []
        while len(todo) != 0:
            longest = None
            length = 0
            for i in todo:
                if self.contigDict[i].length > length:
                    length = self.contigDict[i].length
                    longest = i
            todo.remove(longest)
            currlist = [longest]
            currnode = longest
            currdir = True
            goforward = True
            while goforward:
                goforward = False
                if currdir:
                    for i in self.contigDict[currnode].to:
                        if i[0] in todo:
                            todo.remove(i[0])
                            currlist.append(i[0])
                            currnode = i[0]
                            currdir = i[1]
                            goforward = True
                            break
                else:
                    for i in self.contigDict[currnode].fr:
                        if i[0] in todo:
                            todo.remove(i[0])
                            currlist.append(i[0])
                            currnode = i[0]
                            currdir = i[1]
                            goforward = True
                            break
            currnode = longest
            currdir = False
            goforward = True
            while goforward:
                goforward = False
                if currdir:
                    for i in self.contigDict[currnode].to:
                        if i[0] in todo:
                            todo.remove(i[0])
                            currlist = [i[0]] + currlist
                            currnode = i[0]
                            currdir = i[1]
                            goforward = True
                            break
                else:
                    for i in self.contigDict[currnode].fr:
                        if i[0] in todo:
                            todo.remove(i[0])
                            currlist = [i[0]] + currlist
                            currnode = i[0]
                            currdir = i[1]
                            goforward = True
                            break
            listoflists.append(currlist)
        currpos = 0
        qpos = []
        for i in listoflists:
            for j in i:
                qpos.append((currpos, j))
                currpos += self.contigDict[j].xlength + 10
        self.startviewpos = self.currxscroll/2 - currpos/2
        qpos = map(lambda y: (y[0] + self.startviewpos, y[1]), qpos)
        for i in qpos:
            self.contigDict[i[1]].xpos = i[0]
            self.contigDict[i[1]].ypos = self.contigline

    def orderContigsBlast(self):
        refpos = []
        currx = 0
        for i in self.reforder:
            refpos.append(currx)
            currx += len(self.refDict[i]) / self.scaledown.get() + 10
        refpos = map(lambda x: x + max((self.currxscroll/2 - currx/2, 10)), refpos)
        self.refpos = {}
        for i in range(len(refpos)):
            self.refpos[self.reforder[i]] = refpos[i]
        besthit = {}
        for i in self.hitlist:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = i
            if query in besthit:
                if length > besthit[query][1]:
                    besthit[query] = [self.refpos[subject] + min((rstart, rstop)) / self.scaledown.get(), length]
            else:
                besthit[query] = [self.refpos[subject] + min((rstart, rstop)) / self.scaledown.get(), length]
        qpos = []
        for i in besthit:
            qpos.append((besthit[i][0], i))
        qpos.sort()
        currxpos = qpos[0][0]
        self.startviewpos = currxpos - 20
        for i in qpos:
            self.contigDict[i[1]].xpos = currxpos
            self.contigDict[i[1]].ypos = self.contigline
            currxpos += self.contigDict[i[1]].xlength + 10


    def findPaths(self):
        try:
            self.find_paths.destroy()
        except:
            pass
        self.find_paths = Toplevel()
        self.frame4 = Frame(self.find_paths)
        self.find_paths.geometry('+20+30')
        self.find_paths.title('Find paths between contigs')
        self.maxbplabel = Label(self.frame4, text='Max. bp')
        self.maxbplabel.grid(column=0, row=1)
        self.maxbpentry = Entry(self.frame4, textvariable=self.maxbp)
        self.maxbpentry.grid(column=1, row=1)
        self.maxnodelabel = Label(self.frame4, text='Max. node')
        self.maxnodelabel.grid(column=0, row=2)
        self.maxnodeentry = Entry(self.frame4, textvariable=self.maxnode)
        self.maxnodeentry.grid(column=1, row=2)
        self.maxpathlabel = Label(self.frame4, text='Max. paths')
        self.maxpathlabel.grid(column=0, row=3)
        self.maxpathentry = Entry(self.frame4, textvariable=self.maxpath)
        self.maxpathentry.grid(column=1, row=3)
        self.onlyshortlabel = Label(self.frame4, text='Only find shortest path')
        self.onlyshortlabel.grid(column=0, row=4)
        self.onlyshortentry = Checkbutton(self.frame4, variable=self.onlyshort)
        self.onlyshortentry.grid(column=1, row=4)
        self.okedges = Button(self.frame4, text='Ok', command=self.ok_paths)
        self.okedges.grid(row=5, column=0, sticky=E)
        self.frame4.grid(row=0, column=0)

    def ok_paths(self):
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        self.find_paths.destroy()
        if self.namelist.size() == 0:
            tkMessageBox.showerror('Error', 'Need to select contigs before searching for paths.')
            return
        self.queue = Queue.Queue()
        self.thethread = threading.Thread(target=self.findpathsthread)
        self.thethread.start()
        self.update_path()

    def update_path(self):
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text != 'Path search complete.':
                    root.after(1000, self.update_path)
                else:
                    self.path_add()
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            self.update_console('Path search failed, please check console output.')
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_path)
            return

    def path_add(self):
        for i in self.contigstoadd:
            self.add_contig(i)


    def findpathsthread(self):
        maxbp = self.maxbp.get()
        maxnode = self.maxnode.get()
        maxpath = self.maxpath.get()
        phantomlist = []
        phantomset = set()
        for i in phantomlist:
            phantomset.add(i[0])
        candContigs = set()
        for i in self.namelist.get(0, END):
            candContigs.add(i)
        outpaths = []
        for i in candContigs:
            todo = []
            for j in self.contigDict[i].to:
                if not j[1]:
                    if type(j[2]) == int:
                        overlap = - j[2]
                    else:
                        overlap = len(j[2])
                    todo.append((((i, True), (j[0], False)), self.contigDict[j[0]].length + overlap))
                else:
                    if type(j[2]) == int:
                        overlap = - j[2]
                    else:
                        overlap = len(j[2])
                    todo.append((((i, True), (j[0], True)), self.contigDict[j[0]].length + overlap))
            for j in self.contigDict[i].fr:
                if not j[1]:
                    if type(j[2]) == int:
                        overlap = - j[2]
                    else:
                        overlap = len(j[2])
                    todo.append((((i, False), (j[0], False)), self.contigDict[j[0]].length + overlap))
                else:
                    if type(j[2]) == int:
                        overlap = - j[2]
                    else:
                        overlap = len(j[2])
                    todo.append((((i, False), (j[0], True)), self.contigDict[j[0]].length + overlap))
            while todo != [] and len(todo) < maxpath:
                currPath = todo.pop(0)
                notphant = True
                if currPath[0][-1][0] in phantomset:
                    phantcount = 0
                    for j in currPath[0][:-1]:
                        if j[0] in phantomset:
                            phantcount += 1
                            if phantcount == 4:
                                notphant = False
                                break
                if currPath[0][-1][0] in candContigs:
                    outpaths.append(currPath)
                elif len(currPath[0]) <= maxnode and currPath[1] <= maxbp and notphant:
                    if currPath[0][-1][1]:
                        for j in self.contigDict[currPath[0][-1][0]].to:
                            if not j[1]:
                                if type(j[2]) == int:
                                    overlap = - j[2]
                                else:
                                    overlap = len(j[2])
                                todo.append((currPath[0] + ((j[0], False),), currPath[1] + self.contigDict[j[0]].length + overlap))
                            else:
                                if type(j[2]) == int:
                                    overlap = - j[2]
                                else:
                                    overlap = len(j[2])
                                todo.append((currPath[0] + ((j[0], True),), currPath[1] + self.contigDict[j[0]].length + overlap))
                    else:
                        for j in self.contigDict[currPath[0][-1][0]].fr:
                            if not j[1]:
                                if type(j[2]) == int:
                                    overlap = - j[2]
                                else:
                                    overlap = len(j[2])
                                todo.append((currPath[0] + ((j[0], False),), currPath[1] + self.contigDict[j[0]].length + overlap))
                            else:
                                if type(j[2]) == int:
                                    overlap = - j[2]
                                else:
                                    overlap = len(j[2])
                                todo.append((currPath[0] + ((j[0], True),), currPath[1] + self.contigDict[j[0]].length + overlap))
        if len(outpaths) == 0:
            self.queue.put('No paths found.')
            self.queue.put('Path search complete.')
        paths = []
        pathlengths = []
        for i in outpaths:
            paths.append(i[0])
            pathlengths.append(i[1] - self.contigDict[i[0][-1][0]].length)
        if self.onlyshort.get():
            newpaths = {}
            for i in range(len(paths)):
                pathstart = paths[i][0]
                pathend = paths[i][-1]
                pathlength = pathlengths[i]
                if (pathstart, pathend) in newpaths:
                    if pathlength < newpaths[(pathstart, pathend)][0]:
                        newpaths[(pathstart, pathend)] = (pathlength, paths[i])
                elif ((pathend[0], not pathend[1]), (pathstart[0], not pathstart[1])) in newpaths:
                    if pathlength < newpaths[((pathend[0], not pathend[1]), (pathstart[0], not pathstart[1]))][0]:
                        newpaths[((pathend[0], not pathend[1]), (pathstart[0], not pathstart[1]))] = (pathlength, paths[i])
                else:
                    newpaths[(pathstart, pathend)] = (pathlength, paths[i])
            outpaths = []
            for i in newpaths:
                outpaths.append((newpaths[i][1], newpaths[i][0]))
        self.contigstoadd = []
        for i in outpaths:
            paths, pathlength = i
            firstcoords = self.canvas.coords('c' + paths[0][0])
            secondcoords = self.canvas.coords('c' + paths[-1][0])
            ystart = min([firstcoords[3], secondcoords[3]]) + 35
            if secondcoords[0] < firstcoords[0]:
                paths = list(paths)
                paths.reverse()
                xstart = (secondcoords[2] + firstcoords[0]) /2 - (pathlength / self.newscaledown + 10 * (len(paths) -3)) /2
            else:
                xstart = (firstcoords[2] + secondcoords[0]) /2 - (pathlength / self.newscaledown + 10 * (len(paths) -3)) /2
            for j in paths[1:-1]:
                if not j[0] in self.visible and not j[0] in self.contigstoadd:
                    self.contigDict[j[0]].xpos = xstart
                    self.contigDict[j[0]].ypos = ystart
                    self.contigDict[j[0]].xlength = self.contigDict[j[0]].length / self.newscaledown
                    self.contigstoadd.append(j[0])
                xstart += self.contigDict[j[0]].length / self.newscaledown
        self.queue.put(str(len(outpaths)) + ' paths found.')
        self.queue.put('Path search complete.')

    def writeFasta(self):
        try:
            self.write_fasta.destroy()
        except:
            pass
        self.write_fasta = Toplevel()
        self.frame5 = Frame(self.write_fasta)
        self.write_fasta.geometry('+20+30')
        self.write_fasta.title('Write FASTA')
        self.outfilelabel = Label(self.frame5, text='Write to file')
        self.outfilelabel.grid(column=0, row=1)
        self.outfileentry = Entry(self.frame5, textvariable=self.outfile)
        self.outfileentry.grid(column=1, row=1)
        self.outfilebutton = Button(self.frame5, text='...', command=self.loadoutfile)
        self.outfilebutton.grid(column=2, row=1)
        self.bufferlabel = Label(self.frame5, text='Buffer sequence')
        self.bufferlabel.grid(column=0, row=2)
        self.bufferentry = Entry(self.frame5, textvariable=self.buffer)
        self.bufferentry.grid(column=1, row=2, columnspan=2, sticky=EW)
        self.okedges = Button(self.frame5, text='Ok', command=self.ok_fasta)
        self.okedges.grid(row=3, column=0, sticky=E)
        self.frame5.grid(padx=20, pady=20)


    def ok_fasta(self):
        if self.outfile.get() == '':
            tkMessageBox.showerror('Error', 'Please specify file to write to.')
            return
        self.write_fasta.destroy()
        if self.namelist.size() == 0:
            tkMessageBox.showerror('Error', 'Need to select contigs before writing fasta.')
            return
        out = open(self.outfile.get(), 'w')
        name = self.namelist.get(0)
        if self.dirlist.get(0) == '+':
            lastdir = True
            seq = self.contigDict[name].forseq
        else:
            lastdir = False
            seq = self.contigDict[name].revseq
        for i in range(1, self.namelist.size()):
            lastname = name
            name = self.namelist.get(i)
            overlap = self.buffer.get()
            dir = self.dirlist.get(i) == '+'
            if dir:
                if lastdir:
                    for j in self.contigDict[lastname].to:
                        if j[0] == name and j[1] == dir:
                            overlap = j[2]
                else:
                    for j in self.contigDict[lastname].fr:
                        if j[0] == name and j[1] == dir:
                            overlap = j[2]
                if type(overlap) == int:
                    seq += self.contigDict[name].forseq[overlap:]
                else:
                    seq += overlap + self.contigDict[name].forseq
                lastdir = True
            else:
                if lastdir:
                    for j in self.contigDict[lastname].to:
                        if j[0] == name and j[1] == dir:
                            overlap = j[2]
                else:
                    for j in self.contigDict[lastname].fr:
                        if j[0] == name and j[1] == dir:
                            overlap = j[2]
                if type(overlap) == int:
                    seq += self.contigDict[name].revseq[overlap:]
                else:
                    seq += overlap + self.contigDict[name].revseq
                lastdir = False
        out.write('>contiguity_scaff\n')
        for j in range(0, len(seq), 60):
            out.write(seq[j:j+60] + '\n')
        out.close()
        self.update_console('File written.')

    def writeMultiFasta(self):
        if self.namelist.size() == 0:
            tkMessageBox.showerror('Error', 'Need to select contigs before writing fasta.')
            return
        outfile = tkFileDialog.asksaveasfilename()
        if outfile == '':
            return
        out = open(outfile, 'w')
        for i in range(self.namelist.size()):
            name = self.namelist.get(i)
            out.write('>' + name + '\n')
            if self.dirlist.get(i) == '+':
                seq = self.contigDict[name].forseq
            else:
                seq = self.contigDict[name].revseq
            for j in range(0, len(seq), 60):
                out.write(seq[j:j+60] + '\n')
        out.close()
        self.update_console('File written.')

        
root = Tk()
root.title('Contiguity')
root.option_add('*Font', 'TkDefaultFont 12')
app = App(root)
root.mainloop()