#!/usr/bin/env python

# Contiguity - Tool for visualising assemblies
# Copyright (C) 2013-2015 Mitchell Sullivan
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Mitchell Sullivan
# mjsull@gmail.com
# School of Chemistry & Molecular Biosciences
# The University of Queensland
# Brisbane, QLD 4072.
# Australia


__title__ = 'Contiguity'
__version__ = '1.0.4'
__description__ = "Tool for visualising assemblies"
__author__ = 'Mitchell Sullivan'
__license__ = 'GPLv3'
__author_email__ = "mjsull@gmail.com"
__url__ = 'https://github.com/BeatsonLab-MicrobialGenomics/Contiguity'


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
import argparse
import platform

transtab = string.maketrans('atcgATCG', 'tagcTAGC')

# Instead of pushing a messaged to another thread (reported in console) command line version pushes the message here
# which will print the message to the console
class clqueue:
    def put(self, theval):
        sys.stdout.write(theval + '\n')

# Data from the GUI uses special data types, this mimics it's functionality for the command line
class dummyVar:
    def __init__(self, theval):
        self.theval = theval
    def get(self):
        return self.theval
    def set(self, theval):
        self.theval = theval
    def put(self, theval):
        self.theval = theval

# Contig class, on construction calculates GC content, Skew etc.
class contig:
    def __init__(self, name, shortname, sequence, revseq=None, coverage=None):
        self.name = name
        self.shortname = shortname
        self.forseq = sequence.upper()
        if revseq == None:
            tempseq = self.forseq[::-1]
            self.revseq = tempseq.translate(transtab)
        else:
            self.revseq = revseq
        self.length = len(sequence)
        self.xlength = None
        if self.length >= 1000000000:
            self.strlen = str(round(self.length * 1.0 / 1000000000, 2)) + 'Gb'
        elif self.length >= 1000000:
            self.strlen = str(round(self.length * 1.0 / 1000000, 2)) + 'Mb'
        elif self.length >= 1000:
            self.strlen = str(self.length / 1000) + 'Kb'
        else:
            self.strlen = str(self.length) + 'bp'
        self.visible = False
        self.to = []
        self.fr = []
        if coverage is None:
            self.coverage = 'N/A'
        else:
            self.coverage = round(coverage, 2)
        gcount = self.forseq.count('G')
        ccount = self.forseq.count('C')
        acount = self.forseq.count('A')
        tcount = self.forseq.count('T')
        self.gccontent = round((gcount + ccount) * 100.0 / self.length, 2)
        try:
            self.gcskew = round((gcount - ccount) * 1.0 / (gcount + ccount), 2)
        except ZeroDivisionError:
            self.gcskew = 0
        try:
            self.atskew = round((acount - tcount) * 1.0 / (acount + tcount), 2)
        except ZeroDivisionError:
            self.atskew = 0
        self.customval = 'N/A'
        self.customcol = '#9933CC'
        self.xpos = None
        self.ypos = None
        self.orient = [True]
        self.dupnumber = 1
        self.startanchor = 0
        self.stopanchor = 0


# checks to see if program is in the users path
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True
    return False


# the application
class App:
    def __init__(self, master):
        # if master is a Tk app create the application, otherwise run in command line mode
        if not master is None:
            self.menubar = Menu(master)
            self.filemenu = Menu(self.menubar, tearoff=0)
            self.filemenu.add_command(label="Create CAG file", command=self.create_edges)
            self.filemenu.add_command(label="Create comparison", command=self.create_comp)
            self.filemenu.add_separator()
            self.filemenu.add_command(label="Load assembly", command=self.load_assembly)
            self.filemenu.add_command(label="Save image", command=self.write_canvas_ps)
            self.filemenu.add_separator()
            self.filemenu.add_command(label="Change working directory", command=self.loadwork)
            self.filemenu.add_command(label="Cancel running process", command=self.cancelprocess)
            self.filemenu.add_separator()
            self.filemenu.add_command(label="Exit", command=self.quit)
            self.menubar.add_cascade(label="File", menu=self.filemenu)
            self.viewmenu = Menu(self.menubar, tearoff=0)
            self.viewmenu.add_command(label="View assembly", command=self.view_options)
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
            self.selectmenu = Menu(self.menubar, tearoff=0)
            self.selectmenu.add_command(label="Select all", command=self.select_all)
            self.selectmenu.add_command(label="Clear selected", command=self.clear_lists)
            self.selectmenu.add_separator()
            self.selectmenu.add_command(label="Find paths", command=self.findPaths)
            self.selectmenu.add_command(label="Write fasta", command=self.writeFasta)
            self.selectmenu.add_command(label="Write multifasta", command=self.writeMultiFasta)
            self.selectmenu.add_separator()
            self.selectmenu.add_command(label="Colour options", command=self.colour_options)
            self.selectmenu.add_command(label="Load custom", command=self.load_special)
            self.menubar.add_cascade(label="Tools", menu=self.selectmenu)
            self.currxscroll = 1000000
            self.curryscroll = 1000000
            self.fontsize = 12 # When zooming in/out font does scale, we create a custom font and then change the font size when we zoom
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
            self.utility.add(self.console, sticky=NSEW, height=30, width=180, minsize=100)
            self.contiglist = Frame(self.utility)
            self.utility.add(self.contiglist, sticky=NSEW, minsize=120)
            self.dirframe = Frame(self.utility)
            self.utility.add(self.dirframe, sticky=NSEW, minsize=70)
            self.lengthframe = Frame(self.utility)
            self.utility.add(self.lengthframe, sticky=NSEW, minsize=70)
            self.specialframe = Frame(self.utility)
            self.utility.add(self.specialframe, sticky=NSEW, minsize=70)
            self.contiglist.grid_rowconfigure(1, weight=1, minsize=30)
            self.contiglist.grid_columnconfigure(0, weight=1)
            self.dirframe.grid_rowconfigure(1, weight=1, minsize=30)
            self.dirframe.grid_columnconfigure(0, weight=1)
            self.lengthframe.grid_rowconfigure(1, weight=1, minsize=30)
            self.lengthframe.grid_columnconfigure(0, weight=1)
            self.specialframe.grid_rowconfigure(1, weight=1, minsize=30)
            self.specialframe.grid_columnconfigure(0, weight=1)
            self.clscroll = Scrollbar(self.lengthframe, orient=VERTICAL)
            self.namelabel = Label(self.contiglist, text='Contig name', bg='#FFFFFF',  relief=SUNKEN, pady=4, bd=1)
            self.namelabel.grid(row=0, column=0, sticky=EW, pady=2)
            self.namelist = Listbox(self.contiglist, yscrollcommand=self.clscroll.set, exportselection=0, width=30, height=5)
            self.namelist.bind('<Button-1>', self.setselectedcontig)
            self.namelist.bind('<Double-Button-1>', self.doublecontig)
            self.namelist.bind('<MouseWheel>', self.onmousewheel)
            self.namelist.bind('<Button-4>', self.onmousewheel)
            self.namelist.bind('<Button-5>', self.onmousewheel)
            self.namelist.bind('MouseWheel', self.onmousewheel)
            self.namelist.grid(row=1, column=0, sticky=NSEW)
            self.dirlabel = Label(self.dirframe, bg='#FFFFFF', text='Strand', relief=SUNKEN, pady=4, bd=1)
            self.dirlabel.grid(row=0, column=0, sticky=NSEW, pady=2)
            self.dirlist = Listbox(self.dirframe, yscrollcommand=self.clscroll.set, exportselection=0, width=7, height=5)
            self.dirlist.bind('<Button-1>', self.setselectedcontig)
            self.dirlist.bind('<Double-Button-1>', self.doublecontig)
            self.dirlist.bind('<MouseWheel>', self.onmousewheel)
            self.dirlist.bind('<Button-4>', self.onmousewheel)
            self.dirlist.bind('<Button-5>', self.onmousewheel)
            self.dirlist.bind('MouseWheel', self.onmousewheel)
            self.dirlist.grid(row=1, column=0, sticky=NSEW)
            self.lengthlabel = Label(self.lengthframe, text='Length', bg='#FFFFFF',  relief=SUNKEN, pady=4, bd=1)
            self.lengthlabel.grid(row=0, column=0, sticky=NSEW, pady=2)
            self.lengthlist = Listbox(self.lengthframe, yscrollcommand=self.clscroll.set, exportselection=0, width=7, height=5)
            self.lengthlist.bind('<Button-1>', self.setselectedcontig)
            self.lengthlist.bind('<Double-Button-1>', self.doublecontig)
            self.lengthlist.bind('<MouseWheel>', self.onmousewheel)
            self.lengthlist.bind('<Button-4>', self.onmousewheel)
            self.lengthlist.bind('<Button-5>', self.onmousewheel)
            self.lengthlist.bind('MouseWheel', self.onmousewheel)
            self.lengthlist.grid(row=1, column=0, sticky=NSEW)
            self.currspecial = StringVar(value='None')
            self.speciallabel = OptionMenu(self.specialframe, self.currspecial, 'Coverage', 'GC content', 'Custom', 'GC skew', 'AT skew', 'None', command=self.changespecial)
            self.speciallabel.config(bg='#FFFFFF',  relief=SUNKEN, pady=4, bd=1)
            self.speciallabel.grid(row=0, column=0, sticky=NSEW)
            self.speciallist = Listbox(self.specialframe, yscrollcommand=self.clscroll.set, exportselection=0, width=14, height=5)
            self.speciallist.bind('<Button-1>', self.setselectedcontig)
            self.speciallist.bind('<Double-Button-1>', self.doublecontig)
            self.speciallist.bind('<MouseWheel>', self.onmousewheel)
            self.speciallist.bind('<Button-4>', self.onmousewheel)
            self.speciallist.bind('<Button-5>', self.onmousewheel)
            self.speciallist.bind('MouseWheel', self.onmousewheel)
            self.speciallist.grid(row=1, column=0, sticky=NSEW)
            self.clscroll = Scrollbar(self.specialframe, orient=VERTICAL)
            self.clscroll.config(command=self.yview)
            self.clscroll.grid(row=1, column=1, sticky=NS)
            master.geometry('+10+20')
            root.minsize(600, 400)
            self.canvas.tag_bind('map', '<B1-Motion>', self.dragMove)
            self.canvas.tag_bind('map', '<Button-1>', self.recordMark)
            self.canvas.tag_bind('map', '<Double-Button-1>', self.addtolist)

            self.blastmenu = Menu(root, tearoff=0)
            self.querymenu = Menu(self.blastmenu, tearoff=0)
            self.querymenu.add_command(label="Move", command=self.move_query)
            self.querymenu.add_command(label="Goto", command=self.goto_query)
            self.querymenu.add_command(label="Write", command=self.write_query)
            self.subjectmenu = Menu(self.blastmenu, tearoff=0)
            self.subjectmenu.add_command(label="Move", command=self.move_subject)
            self.subjectmenu.add_command(label="Goto", command=self.goto_subject)
            self.subjectmenu.add_command(label="Write", command=self.write_subject)
            self.blastmenu.add_cascade(label="Query", menu=self.querymenu)
            self.blastmenu.add_cascade(label="Subject", menu=self.subjectmenu)
            self.blastmenu.add_command(label="Remove", command=self.delete_blast)
            self.blastmenu.add_command(label="Details", command=self.showhitblast)
            self.blastmenu.add_command(label="To front", command=self.blasttofront)
            self.blastmenu.add_command(label="To back", command=self.blasttoback)
            self.blasthlmenu = Menu(self.blastmenu, tearoff=0)
            self.blasthlmenu.add_command(label='black', command=self.blast_black)
            self.blasthlmenu.add_command(label='red', command=self.blast_red)
            self.blasthlmenu.add_command(label='blue', command=self.blast_blue)
            self.blasthlmenu.add_command(label='green', command=self.blast_green)
            self.blastmenu.add_cascade(label="Highlight", menu=self.blasthlmenu)
            self.rcmenu = Menu(root, tearoff=0)
            self.rcmenu.add_command(label="Reverse", command=self.reverse_contig)
            self.rcmenu.add_command(label="Remove", command=self.remove_contig)
            self.rcmenu.add_command(label="Duplicate", command=self.duplicate_contig)
            self.rcmenu.add_command(label="Details", command=self.contig_details)
            self.rcmenu.add_separator()
            self.tomenu = Menu(self.rcmenu, tearoff=0)
            self.tohlmenu = Menu(self.tomenu, tearoff=0)
            self.tohlmenu.add_command(label='black', command=self.to_black)
            self.tohlmenu.add_command(label='red', command=self.to_red)
            self.tohlmenu.add_command(label='blue', command=self.to_blue)
            self.tohlmenu.add_command(label='green', command=self.to_green)
            self.tomenu.add_cascade(label="highlight", menu=self.tohlmenu)
            self.tomenu.add_command(label="add", command=self.show_to)
            self.tomenu.add_command(label="move", command=self.move_to)
            self.tomenu.add_command(label="duplicate", command=self.dupe_to)
            self.rcmenu.add_cascade(label="to", menu=self.tomenu)
            self.frmenu = Menu(self.rcmenu, tearoff=0)
            self.frhlmenu = Menu(self.frmenu, tearoff=0)
            self.frhlmenu.add_command(label='black', command=self.fr_black)
            self.frhlmenu.add_command(label='red', command=self.fr_red)
            self.frhlmenu.add_command(label='blue', command=self.fr_blue)
            self.frhlmenu.add_command(label='green', command=self.fr_green)
            self.frmenu.add_cascade(label="highlight", menu=self.frhlmenu)
            self.frmenu.add_command(label="add", command=self.show_fr)
            self.frmenu.add_command(label="move", command=self.move_fr)
            self.frmenu.add_command(label="duplicate", command=self.dupe_fr)
            self.rcmenu.add_cascade(label="fr", menu=self.frmenu)
            if platform.system() in ['Windows', 'Linux']:
                self.canvas.tag_bind('blast', '<Button-3>', self.rcblast)
                self.canvas.tag_bind('selfhit', '<Button-3>', self.rcblast)
                self.canvas.tag_bind('map', '<Button-3>', self.rightClick)
                self.canvas.bind('<Button-2>', self.beginDrag)
                self.canvas.bind('<B2-Motion>', self.dragCanvas)
            else:
                self.canvas.tag_bind('blast', '<Button-2>', self.rcblast)
                self.canvas.tag_bind('selfhit', '<Button-2>', self.rcblast)
                self.canvas.tag_bind('map', '<Button-2>', self.rightClick)
                self.canvas.bind('<Button-3>', self.beginDrag)
                self.canvas.bind('<B3-Motion>', self.dragCanvas)
            root.bind('w', self.zoomin)
            root.bind('s', self.zoomout)
            root.bind('a', self.shrink)
            root.bind('d', self.stretch)
            root.bind('<Left>', self.scrollleft)
            root.bind('<Right>', self.scrollright)
            root.bind('<Up>', self.scrollup)
            root.bind('<Down>', self.scrolldown)
            root.bind('<MouseWheel>', self.zoomcanvas)
            self.canvas.bind('<Button-5>', self.zoomout)
            self.canvas.bind('<Button-4>', self.zoomin)
            self.canvas.bind('<MouseWheel>', self.zoomcanvas)
            self.canvas.bind('<Button-1>', self.removerc)
            self.csagfile = StringVar(value='')
            self.selfcomparefile = StringVar(value='')
            self.outfile = StringVar(value='')
            self.buffer = StringVar(value='nnnnnnnnnn')
            self.reffile = StringVar(value='')
            self.blastfile = StringVar(value='')
            self.maxdist = IntVar(value=300)
            self.nmercut = IntVar(value=4)
            self.nmerave = IntVar(value=10)
            self.nmersize = IntVar(value=31)
            self.ht_size = StringVar(value=args.ht_size)
            self.ht_number = IntVar(value=args.ht_number)
            self.num_threads = IntVar(value=args.num_threads)
            self.minoverlap = IntVar(value=30)
            self.contigfile = StringVar(value='')
            self.readfile = StringVar(value='')
            self.workingDir = StringVar(value='.contiguity_wd')
            if os.path.exists(self.workingDir.get()):
                if not os.path.isdir(self.workingDir.get()):
                    tkMessageBox.showerror('Error', 'File not folder selected.')
                    sys.exit()
            else:
                os.makedirs(self.workingDir.get())
            self.insertsize = IntVar(value=600)
            self.minrl = IntVar(value=75)
            self.minpairedge = IntVar(value=5)
            self.maxmm = IntVar(value=2)
            self.longoveralpident = IntVar(value=85)
            self.minlengthblast = IntVar(value=100)
            self.intra = IntVar(value=0)
            self.onlyedge = IntVar(value=0)
            self.maxedge = IntVar(value=10)
            self.minlengthratio = DoubleVar(value=0.0)
            self.minident = DoubleVar(value=95.0)
            self.minbitscore = DoubleVar(value=0)
            self.maxevalue = DoubleVar(value=0.05)
            self.getOverlap = IntVar(value=1)
            self.getdb = IntVar(value=1)
            self.cutauto = IntVar(value=1)
            self.getPaired = IntVar(value=1)
            self.orientation = StringVar(value='--> <--')
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
            self.gcmin = DoubleVar(value=40.0)
            self.gcmax = DoubleVar(value=60.0)
            self.skewmin = DoubleVar(value=-0.2)
            self.skewmax = DoubleVar(value=0.2)
            self.covmin = DoubleVar(value=0.0)
            self.covmax = DoubleVar(value=200.0)
            self.covlog = IntVar(value=0)
            self.usbase = IntVar(value=0)
            self.dsbase = IntVar(value=0)
            self.bofile = StringVar(value='')
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
            self.reforder = None
            self.abortqueue = Queue.Queue()
            self.abortqueue.put(False)
            self.abort = False
        else: # command line mode
            if args.no_overlap_edges:
                self.getOverlap = dummyVar(0)
            else:
                self.getOverlap = dummyVar(1)
            if args.no_db_edges:
                self.getdb = dummyVar(0)
            else:
                self.getdb = dummyVar(1)
            if args.no_paired_edges:
                self.getPaired = dummyVar(0)
            else:
                self.getPaired = dummyVar(1)
            self.nmercut = dummyVar(args.kmer_average)
            self.nmerave = dummyVar(args.kmer_cutoff)
            self.nmersize = dummyVar(args.kmer_size)
            self.ht_size = dummyVar(args.ht_size)
            self.ht_number = dummyVar(args.ht_number)
            self.num_threads = dummyVar(args.num_threads)
            self.maxdist = dummyVar(args.max_distance)
            if args.kmer_cutoff == -1 or args.kmer_average == -1:
                self.cutauto = dummyVar(1)
            else:
                self.cutauto = dummyVar(0)
            self.contigfile = dummyVar(args.contig_file)
            self.readfile = dummyVar(args.read_file)
            self.insertsize = dummyVar(args.max_insert_size)
            self.minrl = dummyVar(args.min_read_length)
            self.minpairedge = dummyVar(args.minimum_pairs_edge)
            self.maxmm = dummyVar(args.max_mismatch)
            self.longoveralpident = dummyVar(args.long_overlap_ident)
            self.workingDir = dummyVar(args.output_folder)
            if args.overlap is None:
                self.minoverlap = dummyVar(args.kmer_size-1)
            else:
                self.minoverlap = dummyVar(args.minoverlap)
            if os.path.exists(self.workingDir.get()):
                if not os.path.isdir(self.workingDir.get()):
                    sys.stderr.write('Output is not a folder.\n')
                    sys.exit()
            else:
                os.makedirs(self.workingDir.get())
            self.queue = clqueue()
            self.abortqueue = Queue.Queue()
            self.abortqueue.put(False)
            self.abort = False
            self.edge_thread() # Create the CAG

    # when the optonal value in the table is changed this method is run
    def changespecial(self, theopt):
        self.speciallist.delete(0, END)
        if theopt == 'None':
            for i in self.namelist.get(0, END):
                self.speciallist.insert(END, ' ')
            for i in self.canvas.find_withtag('contig'):
                self.canvas.itemconfig(i, fill='#9933CC')
        elif theopt == 'GC content':
            for i in self.namelist.get(0, END):
                self.speciallist.insert(END, self.contigDict[i].gccontent)
            if self.gcmin.get() >= self.gcmax.get() or self.gcmin.get() < 0 or self.gcmax.get() > 100:
                tkMessageBox.showerror('Invalid GC min or max value.', 'Please ensure min < max and values in range 0..100')
                return
            for i in self.canvas.find_withtag('contig'):
                contigname = self.canvas.gettags(i)[0][1:]
                if '_dup' in contigname:
                    contigname = contigname.split('_dup')[0]
                gccontent = self.contigDict[contigname].gccontent
                if gccontent <= self.gcmin.get():
                    colour = '#9933CC'
                elif gccontent >= self.gcmax.get():
                    colour = '#CC0000'
                else:
                    ratio = (gccontent - self.gcmin.get()) / (self.gcmax.get() - self.gcmin.get())
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
                self.canvas.itemconfig(i, fill=colour)
        elif theopt == 'GC skew':
            for i in self.namelist.get(0, END):
                self.speciallist.insert(END, self.contigDict[i].gcskew)
            if self.skewmin.get() >= self.skewmax.get() or self.skewmin.get() < -1 or self.skewmax.get() > 1:
                tkMessageBox.showerror('Invalid skew min or max value.', 'Please ensure min < max and values in range -1..1')
                return
            for i in self.canvas.find_withtag('contig'):
                contigname = self.canvas.gettags(i)[0][1:]
                if '_dup' in contigname:
                    contigname = contigname.split('_dup')[0]
                    orient = self.contigDict[contigname].orient[int(contigname.split('_dup')[1])]
                else:
                    orient = self.contigDict[contigname].orient[0]
                if orient:
                    gcskew = self.contigDict[contigname].gcskew
                else:
                    gcskew = - self.contigDict[contigname].gccontent
                if gcskew <= self.skewmin.get():
                    colour = '#9933CC'
                elif gcskew >= self.skewmax.get():
                    colour = '#CC0000'
                else:
                    ratio = (gcskew - self.skewmin.get()) / (self.skewmax.get() - self.skewmin.get())
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
                self.canvas.itemconfig(i, fill=colour)
        elif theopt == 'AT skew':
            for i in self.namelist.get(0, END):
                self.speciallist.insert(END, self.contigDict[i].atskew)
            if self.skewmin.get() >= self.skewmax.get() or self.skewmin.get() < -1 or self.skewmax.get() > 1:
                tkMessageBox.showerror('Invalid skew min or max value.', 'Please ensure min < max and values in range -1..1')
                return
            for i in self.canvas.find_withtag('contig'):
                contigname = self.canvas.gettags(i)[0][1:]
                if '_dup' in contigname:
                    contigname = contigname.split('_dup')[0]
                    orient = self.contigDict[contigname].orient[int(contigname.split('_dup')[1])]
                else:
                    orient = self.contigDict[contigname].orient[0]
                if orient:
                    atskew = self.contigDict[contigname].atskew
                else:
                    atskew = - self.contigDict[contigname].gccontent
                if atskew <= self.skewmin.get():
                    colour = '#9933CC'
                elif atskew >= self.skewmax.get():
                    colour = '#CC0000'
                else:
                    ratio = (atskew - self.skewmin.get()) / (self.skewmax.get() - self.skewmin.get())
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
                self.canvas.itemconfig(i, fill=colour)
        elif theopt == 'Coverage':
            for i in self.namelist.get(0, END):
                self.speciallist.insert(END, self.contigDict[i].coverage)
            if self.covmin.get() >= self.covmax.get() or self.covmin.get() < 0:
                tkMessageBox.showerror('Invalid skew min or max value.', 'Please ensure min < max and values in range -1..1')
                return
            if self.covmin.get() == 0 and self.covlog.get() == 1:
                self.covmin.set(1)
            for i in self.canvas.find_withtag('contig'):
                contigname = self.canvas.gettags(i)[0][1:]
                if '_dup' in contigname:
                    contigname = contigname.split('_dup')[0]
                coverage = self.contigDict[contigname].coverage
                if coverage <= self.covmin.get() or coverage == 'N/A':
                    colour = '#9933CC'
                elif coverage >= self.covmax.get():
                    colour = '#CC0000'
                else:
                    if self.covlog.get() == 0:
                        ratio = (coverage - self.covmin.get()) / (self.covmax.get() - self.covmin.get())
                    else:
                        ratio = (math.log10(coverage) - math.log10(self.covmin.get())) / (math.log10(self.covmax.get()) - math.log10(self.covmin.get()))
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
                self.canvas.itemconfig(i, fill=colour)
        elif theopt == 'Custom':
            for i in self.namelist.get(0, END):
                self.speciallist.insert(END, self.contigDict[i].customval)
            for i in self.canvas.find_withtag('contig'):
                self.canvas.itemconfig(i, fill='#9933CC')

    # send message to running threads to terminate
    def cancelprocess(self):
        self.abort = True
        self.abortqueue.put(True)

    # on quit terminate running processes
    def quit(self):
        self.cancelprocess()
        self.update_console('Terminating threads.')
        self.quitpoll()

    # check to see if threads are dead before quiting
    def quitpoll(self):
        self.dot_console()
        try:
            if self.thethread.is_alive():
                root.after(1000, self.quitpoll) # check to see if thread is still alive in 1 second
            else:
                root.quit()
        except: # or if thread is not found quit
            root.quit()

    # select all visible contigs
    def select_all(self):
        # remove all currently selected contigs
        self.namelist.delete(0, END)
        self.dirlist.delete(0, END)
        self.lengthlist.delete(0, END)
        self.speciallist.delete(0, END)
        self.selected = []
        contiglist = self.canvas.find_withtag('contig')
        newlist = []
        # highlight all contigs and then find order
        for i in contiglist:
            self.canvas.itemconfig(i, outline='#FFCC00', width=3)
            coords = self.canvas.coords(i)
            thetag = self.canvas.gettags(i)[0]
            newlist.append((coords, thetag))
        newlist.sort() # sort contigs from left to right
        # add all contigs to the table
        for i in newlist:
            thetag = i[1]
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
            if self.currspecial.get() == 'None':
                self.speciallist.insert(END, ' ')
            elif self.currspecial.get() == 'GC content':
                self.speciallist.insert(END, self.contigDict[contig].gccontent)
            elif self.currspecial.get() == 'GC skew':
                self.speciallist.insert(END, self.contigDict[contig].gcskew)
            elif self.currspecial.get() == 'AT skew':
                self.speciallist.insert(END, self.contigDict[contig].atskew)
            elif self.currspecial.get() == 'Coverage':
                self.speciallist.insert(END, self.contigDict[contig].coverage)
            elif self.currspecial.get() == 'Custom':
                self.speciallist.insert(END, self.contigDict[contig].custom)

    # unselect all contigs
    def clear_lists(self):
        self.namelist.delete(0, END)
        self.dirlist.delete(0, END)
        self.lengthlist.delete(0, END)
        self.speciallist.delete(0, END)
        self.selected = []
        contiglist = self.canvas.find_withtag('contig')
        for i in contiglist:
            self.canvas.itemconfig(i, outline='#000000', width=1) # unhighlight all contigs


    # load special values and colours for each contig
    def load_special(self):
        filename = tkFileDialog.askopenfilename()
        if filename == '' or filename == ():
            return
        cset = set('abcdef0123456789ABCDEF')
        specfile = open(filename)
        for line in specfile:
            splitline = line.split()
            if splitline[0] in conitgDict:
                self.contigDict[splitline[0]].customval = splitline[1]
            else:
                self.update_console(splitline[0] + ' not found.')
            if len(splitline) == 3:
                if len(splitline[2]) == 7 and splitline[2][0] == '#' and set(splitline[2][1:]) <= cset:
                    self.contigDict[splitline[0]].customcol = splitline[2]
                else:
                    self.update_console('Invalid colour')
            elif len(splitline) != 2:
                self.update_console('Invalid column number')

    # select contig
    def addtolist(self, event):
        thetag = self.canvas.gettags(CURRENT)[0]
        thecontig = self.canvas.find_withtag(thetag)[0]
        thecolour = self.canvas.itemcget(thecontig, "outline")
        # check the colour of the contigs outline to see if already selected
        if thecolour == "#000000":
            self.canvas.itemconfig(thecontig, outline='#FFCC00', width=3)
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
            if self.currspecial.get() == 'None':
                self.speciallist.insert(END, ' ')
            elif self.currspecial.get() == 'GC content':
                self.speciallist.insert(END, self.contigDict[contig].gccontent)
            elif self.currspecial.get() == 'GC skew':
                self.speciallist.insert(END, self.contigDict[contig].gcskew)
            elif self.currspecial.get() == 'AT skew':
                self.speciallist.insert(END, self.contigDict[contig].atskew)
            elif self.currspecial.get() == 'Coverage':
                self.speciallist.insert(END, self.contigDict[contig].coverage)
            elif self.currspecial.get() == 'Custom':
                self.speciallist.insert(END, self.contigDict[contig].custom)
            self.selected.append(thetag)
        else: # if contig selected find position in list and remove from table
            self.canvas.itemconfig(thecontig, outline='#000000', width=1)
            for i in range(len(self.selected)):
                if self.selected[i] == thetag:
                    self.namelist.delete(i)
                    self.dirlist.delete(i)
                    self.lengthlist.delete(i)
                    self.speciallist.delete(i)
                    self.selected.pop(i)
                    break

    # Post the blast menu if hit is right clicked
    def rcblast(self, event):
        self.rctag = self.canvas.gettags(CURRENT)[0]
        self.rcitem = self.canvas.find_withtag(CURRENT)[0]
        self.blastmenu.post(event.x_root, event.y_root)
        self.rcpos = (event.x_root, event.y_root)

    # move the query of the hit to be aligned with the reference
    def move_query(self):
        bcoords = self.canvas.coords(self.rcitem)
        if len(bcoords) == 8:
            qleftmost = min([bcoords[4], bcoords[6]])
            sleftmost = min([bcoords[0], bcoords[2]])
            ctag = self.rctag[:-1]
            ymod = 0
        if len(bcoords) == 20:
            qleftmost = min([bcoords[0], bcoords[16]])
            sleftmost = min([bcoords[6], bcoords[10]])
            ctag = self.rctag[:-2]
            if bcoords[1] == bcoords[19]:
                ymod = self.contigheight * 2
            else:
                ymod = 0
        movex = sleftmost - qleftmost
        ccords = self.canvas.coords(ctag)
        self.move_contig(ctag[1:], ccords[0] + movex, ccords[1] + ymod)

    # move to where the query is on the canvas
    def goto_query(self):
        if self.canvas.gettags(self.rcitem)[-2] == 'selfhit':
            ctag = self.rctag[1:-2]
        else:
            ctag = self.rctag[1:-1]
        self.goto(ctag)

    # Open a window with the options for writing hit sequence to file
    def write_blast_window(self):
        try:
            self.write_blast_top.destroy()
        except AttributeError:
            pass
        self.write_blast_top = Toplevel()
        self.write_blast_top.grab_set()
        self.write_blast_top.wm_attributes("-topmost", 1)
        x = root.winfo_rootx() + 15
        y = root.winfo_rooty() - 20
        geom = "+%d+%d" % (x,y)
        self.write_blast_top.geometry(geom)
        self.write_blast_top.title('Write BLAST')
        self.write_blast_frame = Frame(self.write_blast_top)
        self.bofilelabel = Label(self.write_blast_frame, text='Output file:', anchor=E)
        self.bofilelabel.grid(column=0, row=1, sticky=E)
        self.bofileentry = Entry(self.write_blast_frame, textvariable=self.bofile)
        self.bofileentry.grid(column=1, row=1)
        self.bofilebut = Button(self.write_blast_frame, text='...', command=self.loadbo)
        self.bofilebut.grid(column=2, row=1)
        self.usbaselabel = Label(self.write_blast_frame, text='Upstream bases:', anchor=E)
        self.usbaselabel.grid(column=0, row=2, sticky=E)
        self.usbaseentry = Entry(self.write_blast_frame, textvariable=self.usbase)
        self.usbaseentry.grid(column=1, row=2, columnspan=2, sticky=EW)
        self.dsbaselabel = Label(self.write_blast_frame, text='Downstream bases', anchor=E)
        self.dsbaselabel.grid(column=0, row=3, sticky=E)
        self.dsbaseentry = Entry(self.write_blast_frame, textvariable=self.dsbase)
        self.dsbaseentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.closeblastbut = Button(self.write_blast_frame, text='Ok', command=self.okwriteblast)
        self.closeblastbut.grid(column=2, row=4, sticky=E, pady=5)
        self.write_blast_frame.grid(padx=10, pady=10)

    # if the ok button is hit write sequence to file
    def okwriteblast(self):
        try:
            out = open(self.bofile.get(), 'w')
        except IOError:
            tkMessageBox.showerror('File not valid', 'Please choose a new file to write to.')
        try:
            self.write_blast_top.destroy()
        except AttributeError:
            pass
        start = self.bwstart - self.usbase.get()
        stop = self.bwstop + self.dsbase.get()
        if start < 0:
            start = 0
        if self.bwtype == 'query':
            seq = self.contigDict[self.bwcontig].forseq[start:stop]
        else:
            seq = self.refDict[self.bwcontig][start:stop]
        out.write('>contiguity_blast_hit\n')
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + '\n')

    # get the contig name, start, stop of the hit and open the write blast to fasta window
    def write_query(self):
        thetag = self.canvas.gettags(self.rcitem)[-1]
        if self.canvas.gettags(self.rcitem)[-2] == 'selfhit':
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = self.selfhit[int(thetag) -1]
        else:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = self.hitlist[int(thetag) -1]
        self.bwtype = 'query'
        self.bwstart = qstart
        self.bwstop = qstop
        self.bwcontig = query
        self.write_blast_window()

    # move the subject of the hit so it is in line with the query
    def move_subject(self):
        bcoords = self.canvas.coords(self.rcitem)
        thetag = self.canvas.gettags(self.rcitem)[1]
        if len(bcoords) == 8:
            qleftmost = min([bcoords[4], bcoords[6]])
            sleftmost = min([bcoords[0], bcoords[2]])
            ctag = thetag[:-1]
            ymod = 0
        if len(bcoords) == 20:
            qleftmost = min([bcoords[0], bcoords[16]])
            sleftmost = min([bcoords[6], bcoords[10]])
            ctag = thetag[:-2]
            if bcoords[1] == bcoords[19]:
                ymod = self.contigheight * 2
            else:
                ymod = 0
        movex = qleftmost - sleftmost
        ccords = self.canvas.coords(ctag)
        # check to see if self hit, or hit against reference
        if len(bcoords) == 20:
            self.move_contig(ctag[1:], ccords[0] + movex, ccords[1] + ymod)
        else:
            self.move_contig(ctag[1:], ccords[0] + movex, ccords[1] + ymod, True)

    # move canvas to the subject of the hit
    def goto_subject(self):
        thetag = self.canvas.gettags(self.rcitem)[1]
        if self.canvas.gettags(self.rcitem)[-2] == 'selfhit':
            ctag = thetag[1:-2]
            self.goto(ctag)
        else:
            ctag = thetag[1:-1]
            self.goto(ctag, True)

    # Find the subject name, start, stop and open write hit to fasta window
    def write_subject(self):
        thetag = self.canvas.gettags(self.rcitem)[-1]
        if self.canvas.gettags(self.rcitem)[-2] == 'selfhit':
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = self.selfhit[int(thetag) -1]
            self.bwtype = 'query' # tells write fasta to get sequence from contig sequences not reference sequences (bit confusing really)
        else:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = self.hitlist[int(thetag) -1]
            self.bwtype = 'subject'
        self.bwstart = rstart
        self.bwstop = rstop
        self.bwcontig = subject
        self.write_blast_window()

    # move hit to the foreground
    def blasttofront(self):
        self.canvas.lift(self.rcitem)

    # move hit to the background
    def blasttoback(self):
        self.canvas.lower(self.rcitem)

    # make blast outline black (normal)
    def blast_black(self):
        self.canvas.itemconfig(self.rcitem, outline='black', width=1)

    # highlight blast outline (blue)
    def blast_blue(self):
        self.canvas.itemconfig(self.rcitem, outline='blue', width=3)

    def blast_red(self):
        self.canvas.itemconfig(self.rcitem, outline='red', width=3)

    def blast_green(self):
        self.canvas.itemconfig(self.rcitem, outline='green', width=3)

    # show details about the blast hit
    def showhitblast(self):
        try:
            self.hitwindow.destroy()
        except:
            pass
        self.hitwindow = Toplevel()
        self.hitwindow.wm_attributes("-topmost", 1)
        self.hitframe = Frame(self.hitwindow)
        self.hitwindow.geometry('+20+30')
        self.hitwindow.title('BLAST hit')
        thetag = self.canvas.gettags(self.rcitem)[-1]
        if thetag == 'current':
            thetag = self.canvas.gettags(self.rcitem)[-2]
        if self.canvas.gettags(self.rcitem)[-2] == 'selfhit':
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = map(str, self.selfhit[int(thetag) -1])
        else:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = map(str, self.hitlist[int(thetag) -1])
        self.hl1 = Label(self.hitframe, text='Query:', anchor=E)
        self.hl1.grid(column=0, row=1)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=query), state='readonly')
        self.he1.grid(column=1, row=1)
        self.hl1 = Label(self.hitframe, text='Subject:', anchor=E)
        self.hl1.grid(column=0, row=2)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=subject), state='readonly')
        self.he1.grid(column=1, row=2)
        self.hl1 = Label(self.hitframe, text='Identity:', anchor=E)
        self.hl1.grid(column=0, row=3)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=ident), state='readonly')
        self.he1.grid(column=1, row=3)
        self.hl1 = Label(self.hitframe, text='Length:', anchor=E)
        self.hl1.grid(column=0, row=4)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=length), state='readonly')
        self.he1.grid(column=1, row=4)
        self.hl1 = Label(self.hitframe, text='Mismatches:', anchor=E)
        self.hl1.grid(column=0, row=5)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=mm), state='readonly')
        self.he1.grid(column=1, row=5)
        self.hl1 = Label(self.hitframe, text='Indels:', anchor=E)
        self.hl1.grid(column=0, row=6)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=indel), state='readonly')
        self.he1.grid(column=1, row=6)
        self.hl1 = Label(self.hitframe, text='Query start:', anchor=E)
        self.hl1.grid(column=0, row=7)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=qstart), state='readonly')
        self.he1.grid(column=1, row=7)
        self.hl1 = Label(self.hitframe, text='Query end:', anchor=E)
        self.hl1.grid(column=0, row=8)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=qstop), state='readonly')
        self.he1.grid(column=1, row=8)
        self.hl1 = Label(self.hitframe, text='Subject start:', anchor=E)
        self.hl1.grid(column=0, row=9)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=rstart), state='readonly')
        self.he1.grid(column=1, row=9)
        self.hl1 = Label(self.hitframe, text='Subject end:', anchor=E)
        self.hl1.grid(column=0, row=10)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=rstop), state='readonly')
        self.he1.grid(column=1, row=10)
        self.hl1 = Label(self.hitframe, text='E. value:', anchor=E)
        self.hl1.grid(column=0, row=11)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=eval), state='readonly')
        self.he1.grid(column=1, row=11)
        self.hl1 = Label(self.hitframe, text='Bitscore:', anchor=E)
        self.hl1.grid(column=0, row=12)
        self.he1 = Entry(self.hitframe, textvariable=StringVar(value=bitscore), state='readonly')
        self.he1.grid(column=1, row=12)
        self.hitframe.grid(padx=10, pady=10)

    # show details about the contig
    def contig_details(self):
        try:
            self.cd_top.destroy()
        except AttributeError:
            pass
        self.cd_top = Toplevel()
        self.cd_top.grab_set()
        self.cd_top.wm_attributes("-topmost", 1)
        x = root.winfo_rootx() + 15
        y = root.winfo_rooty() - 20
        geom = "+%d+%d" % (x,y)
        self.cd_top.geometry(geom)
        self.cd_top.title('Contig details')
        if '_dup' in self.rctag:
            contig = self.rctag.split('_dup')[0][1:]
        else:
            contig = self.rctag[1:]
        self.cd_frame = Frame(self.cd_top)
        self.cl1 = Label(self.cd_frame, text='Name:', anchor=E)
        self.cl1.grid(column=0, row=1)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=contig), state='readonly')
        self.ce1.grid(column=1, row=1)
        self.cl1 = Label(self.cd_frame, text='Short name:', anchor=E)
        self.cl1.grid(column=0, row=2)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=self.contigDict[contig].shortname), state='readonly')
        self.ce1.grid(column=1, row=2)
        self.cl1 = Label(self.cd_frame, text='GC content:', anchor=E)
        self.cl1.grid(column=0, row=3)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=self.contigDict[contig].gccontent), state='readonly')
        self.ce1.grid(column=1, row=3)
        self.cl1 = Label(self.cd_frame, text='GC Skew:', anchor=E)
        self.cl1.grid(column=0, row=4)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=self.contigDict[contig].gcskew), state='readonly')
        self.ce1.grid(column=1, row=4)
        self.cl1 = Label(self.cd_frame, text='AT Skew:', anchor=E)
        self.cl1.grid(column=0, row=5)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=self.contigDict[contig].atskew), state='readonly')
        self.ce1.grid(column=1, row=5)
        self.cl1 = Label(self.cd_frame, text='Coverage:', anchor=E)
        self.cl1.grid(column=0, row=6)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=self.contigDict[contig].coverage), state='readonly')
        self.ce1.grid(column=1, row=6)
        self.cl1 = Label(self.cd_frame, text='Length:', anchor=E)
        self.cl1.grid(column=0, row=7)
        self.ce1 = Entry(self.cd_frame, textvariable=StringVar(value=self.contigDict[contig].length), state='readonly')
        self.ce1.grid(column=1, row=7)
        self.cl1 = Label(self.cd_frame, text='To contigs:')
        self.cl1.grid(column=1, row=8)
        self.tolist = Listbox(self.cd_frame, height=max([len(self.contigDict[contig].to), len(self.contigDict[contig].fr), 10]))
        for i in self.contigDict[contig].to:
            if i[1]:
                self.tolist.insert(END, str(self.contigDict[i[0]].shortname))
            else:
                self.tolist.insert(END, '-' + str(self.contigDict[i[0]].shortname))
        self.frlist = Listbox(self.cd_frame, height=max([len(self.contigDict[contig].to), len(self.contigDict[contig].fr), 10]))
        for i in self.contigDict[contig].fr:
            if i[1]:
                self.frlist.insert(END, str(self.contigDict[i[0]].shortname))
            else:
                self.frlist.insert(END, '-' + str(self.contigDict[i[0]].shortname))
        self.tolist.grid(column=1, row=9)
        self.cl1 = Label(self.cd_frame, text='From Contigs:')
        self.cl1.grid(column=0, row=8)
        self.frlist.grid(column=0, row=9)
        self.cd_frame.grid(padx=10, pady=10)

    # delte blast hit
    def delete_blast(self):
        self.canvas.delete(self.rcitem)

    # write canvas image to selcript file
    def write_canvas_ps(self):
        try:
            saveas = tkFileDialog.asksaveasfilename(parent=root)
        except IOError:
            tkMessageBox.showerror('File not valid', 'Please choose another file.')
        self.canvas.postscript(file=saveas, colormode='color')

    # Couldn't get scroll bar working properly with windows, bit of a hackity hack unfortuanately
    def yview(self, *args):
        self.namelist.yview_scroll(args[0], 'units')
        self.dirlist.yview_scroll(args[0], 'units')
        self.lengthlist.yview_scroll(args[0], 'units')
        self.speciallist.yview_scroll(args[0], 'units')

    # scroll all lists in table on mouse wheel
    def onmousewheel(self, event):
        if event.num == 5 or event.delta == -120:
            self.namelist.yview_scroll(1, 'units')
            self.dirlist.yview_scroll(1, 'units')
            self.lengthlist.yview_scroll(1, 'units')
            self.speciallist.yview_scroll(1, 'units')
        if event.num == 4 or event.delta == 120:
            self.namelist.yview_scroll(-1, 'units')
            self.dirlist.yview_scroll(-1, 'units')
            self.lengthlist.yview_scroll(-1, 'units')
            self.speciallist.yview_scroll(-1, 'units')
        return "break" # scrolling in a window automatically scrolls that window by one, returning break stops this behaviour

    # what to do when contig in list is clicked (may have some functionality later)
    def setselectedcontig(self, event):
        pass

    # when contig in list is double clicked move to contig on canvas
    def doublecontig(self, event):
        x = self.namelist.nearest(event.y)
        self.goto(self.namelist.get(x))

    # on starting mouse drag record original position of contig
    def beginDrag(self, event):
        self.canvas.scan_mark(event.x, event.y)  # initial middle-mouse click

    # capture dragging
    def dragCanvas(self, event):
        self.canvas.scan_dragto(event.x, event.y, 1)

    # if contig is right clicked post menu
    def rightClick(self, event):
        self.rctag = self.canvas.gettags(CURRENT)[0]
        self.rcmenu.post(event.x_root, event.y_root)
        self.rcpos = (event.x_root, event.y_root)

    # if canvas is clicked remove posted menus
    def removerc(self, event):
        self.rcmenu.unpost()
        self.blastmenu.unpost()

    # remove contig and assosciated edges, blast hits, if contig is selected unselect
    def remove_contig(self):
        thetag = self.rctag
        starttag = thetag + 's'
        endtag = thetag + 'e'
        btag = thetag + 'b'
        sstag = thetag + 'ss'
        setag = thetag + 'se'
        thecolour = self.canvas.itemcget(thetag, "outline")
        if thetag[0] == 'c':
            self.visible.remove(thetag[1:])
        self.canvas.delete(thetag)
        self.canvas.delete(starttag)
        self.canvas.delete(endtag)
        self.canvas.delete(btag)
        self.canvas.delete(sstag)
        self.canvas.delete(setag)
        if thecolour == "#FFCC00":
            for i in range(len(self.selected)):
                if self.selected[i] == thetag:
                    self.namelist.delete(i)
                    self.dirlist.delete(i)
                    self.lengthlist.delete(i)
                    self.speciallist.delete(i)
                    self.selected.pop(i)
                    break

    # flip the contig on the canvas
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
        if self.currspecial.get() == 'GC skew':
            if self.contigDict[contig].orient[dupnum]:
                gcskew = self.contigDict[contig].gcskew
            else:
                gcskew = - self.contigDict[contig].gcskew
            if gcskew <= self.skewmin.get():
                colour = '#9933CC'
            elif gcskew >= self.skewmax.get():
                colour = '#CC0000'
            else:
                ratio = (gcskew - self.skewmin.get()) / (self.skewmax.get() - self.skewmin.get())
                colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
            itema = self.canvas.find_withtag(thetag)
            itemb = self.canvas.find_withtag('contig')
            for item in itema:
                if item in itemb:
                    self.canvas.itemconfig(item, fill=colour)
        if self.currspecial.get() == 'AT skew':
            if self.contigDict[contig].orient[dupnum]:
                atskew = self.contigDict[contig].atskew
            else:
                atskew = - self.contigDict[contig].atskew
            if atskew <= self.skewmin.get():
                colour = '#9933CC'
            elif atskew >= self.skewmax.get():
                colour = '#CC0000'
            else:
                ratio = (atskew - self.skewmin.get()) / (self.skewmax.get() - self.skewmin.get())
                colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
            itema = self.canvas.find_withtag(thetag)
            itemb = self.canvas.find_withtag('contig')
            for item in itema:
                if item in itemb:
                    self.canvas.itemconfig(item, fill=colour)
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
        htag = thetag + 'ss'
        selfhits = self.canvas.find_withtag(htag)
        for i in selfhits:
            x = self.canvas.coords(i)
            newx1 = start + end - x[0]
            newx2 = start + end - x[16]
            self.canvas.coords(i, newx1, x[1], newx1, x[1], (newx1 + x[6]) / 2, abs(newx1 - x[6]) /4 + (x[1] + x[7]) / 2, x[6], x[7],
                                x[6], x[7], x[10], x[11], x[10], x[11], (newx2 + x[10]) / 2, abs(newx2 - x[10]) /4 + (x[1] + x[7]) / 2,
                                newx2, x[17], newx2, x[17])
            thecolour = self.canvas.itemcget(i, "fill")
            if thecolour.startswith('#ff'):
                ratio = (131 - int(thecolour[5:], 16)) * 1.0 / 131
                thecolour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
            else:
                ratio = (168 - int(thecolour[1:3], 16)) * 1.0 / 168
                thecolour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
            self.canvas.itemconfig(i, fill=thecolour)
        htag = thetag + 'se'
        selfhits = self.canvas.find_withtag(htag)
        for i in selfhits:
            x = self.canvas.coords(i)
            newx1 = start + end - x[6]
            newx2 = start + end - x[10]
            self.canvas.coords(i, x[0], x[1], x[0], x[1], (x[0] + newx1) / 2, abs(x[0] - newx1) /4 + (x[1] + x[7]) / 2, newx1, x[7],
                                newx1, x[7], newx2, x[11], newx2, x[11], (x[16] + newx2) / 2,
                                abs(x[16] - (newx2)) /4 + (x[1] + x[7]) / 2, x[16], x[17], x[16], x[17])
            thecolour = self.canvas.itemcget(i, "fill")
            if thecolour.startswith('#ff'):
                ratio = (131 - int(thecolour[5:], 16)) * 1.0 / 131
                thecolour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
            else:
                ratio = (168 - int(thecolour[1:3], 16)) * 1.0 / 168
                thecolour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
            self.canvas.itemconfig(i, fill=thecolour)
        endtag = thetag + 'b'
        bhits = self.canvas.find_withtag(endtag)
        if endtag[0] == 'r':
            for i in bhits:
                x = self.canvas.coords(i)
                self.canvas.coords(i, end - (x[0] - start), x[1], end - (x[2] - start), x[3], x[4], x[5], x[6], x[7])
                thecolour = self.canvas.itemcget(i, "fill")
                if thecolour.startswith('#ff'):
                    ratio = (131 - int(thecolour[5:], 16)) * 1.0 / 131
                    thecolour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
                else:
                    ratio = (168 - int(thecolour[1:3], 16)) * 1.0 / 168
                    thecolour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
                self.canvas.itemconfig(i, fill=thecolour)
        else:
            for i in bhits:
                x = self.canvas.coords(i)
                self.canvas.coords(i, x[0], x[1], x[2], x[3], end - (x[4] - start), x[5], end - (x[6] - start), x[7])
                thecolour = self.canvas.itemcget(i, "fill")
                if thecolour.startswith('#ff'):
                    ratio = (131 - int(thecolour[5:], 16)) * 1.0 / 131
                    thecolour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
                else:
                    ratio = (168 - int(thecolour[1:3], 16)) * 1.0 / 168
                    thecolour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
                self.canvas.itemconfig(i, fill=thecolour)
        texttag = thetag + 't'
        try:
            textitem = self.canvas.find_withtag(texttag)[0]
            thetext = self.canvas.itemcget(textitem, "text")
            namelen = len(self.contigDict[contig].shortname)
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


    # duplicate the contig on the canvas
    def duplicate_contig(self, zecontig=None, zecoords=None):
        if zecontig == None:
            thetag = self.rctag
        else:
            thetag = 'c' + zecontig
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
        if not zecoords is None:
            modx = zecoords[0] - bb[0]
            mody = zecoords[1] - bb[1]
        else:
            modx = int(1.5*self.contigheight)
            mody = int(1.5*self.contigheight)
        newtag = 'c' + contigtag + '_dup' + str(self.contigDict[contigtag].dupnumber)
        self.visible.add(newtag[1:])
        self.contigDict[contigtag].dupnumber += 1
        self.contigDict[contigtag].orient.append(self.contigDict[contigtag].orient[dupnum])
        btag = thetag + 'b'
        bhits = self.canvas.find_withtag(btag)
        if btag[0] == 'r':
            pass
        else:
            for i in bhits:
                x = self.canvas.coords(i)
                colour = self.canvas.itemcget(i, "fill")
                self.canvas.create_polygon(x[0], x[1], x[2], x[3], x[4] + modx, x[5] + mody, x[6] + modx, x[7] + mody, \
                                           fill=colour, outline="black", tags=(newtag + 'b', self.canvas.gettags(i)[1], 'blast', self.canvas.gettags(i)[-1]))
        colour = self.canvas.itemcget(thetag, "fill")
        self.canvas.create_rectangle(bb[0] + modx, bb[1] + mody, bb[2] + modx, bb[3] + mody, fill=colour, tags=(newtag, 'contig', 'map'))
        if bb[3] + mody > self.curryscroll:
            self.curryscroll = bb[3] + mody + 20
            self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        starttag = thetag + 'ss'
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            colour = self.canvas.itemcget(i, "fill")
            thetags = self.canvas.itemcget(i, "tags")
            x = self.canvas.coords(i)
            self.canvas.create_polygon(x[0] + modx, x[1] + mody, x[0] + modx, x[1] + mody, (x[0] + modx + x[6]) / 2, abs(x[0] + modx - x[6]) /4 + (x[1] + mody + x[7]) / 2, x[6], x[7],
                                x[6], x[7], x[10], x[11], x[10], x[11], (x[16] + modx + x[10]) / 2, abs(x[16] + modx - x[10]) /4 + (x[1] + mody + x[7]) / 2,
                                x[16] + modx, x[17] + mody, x[16] + modx, x[17] + mody, smooth=1, fill=colour, outline='black', tags=(newtag + 'ss', thetags[0], thetags[2], thetags[3]))
        endtag = thetag + 'se'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            colour = self.canvas.itemcget(i, "fill")
            thetags = self.canvas.itemcget(i, "tags")
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0], x[1], x[0], x[1], (x[0] + x[6] + modx) / 2, abs(x[0] - (x[6] + modx)) /4 + (x[1] + x[7] + mody) / 2, x[6] + modx, x[7] + mody,
                                x[6] + modx, x[7] + mody, x[10] + modx, x[11] + mody, x[10] + modx, x[11] + mody, (x[16] + x[10] + modx) / 2,
                                abs(x[16] - (x[10] + modx)) /4 + (x[1] + x[7] + mody) / 2, x[16], x[17], x[16], x[17], smooth=1, fill=colour, outline='black',
                                tags=(thetags[0], newtag + 'se', thetags[2], thetags[3]))
        starttag = thetag + 's'
        arcs= self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.create_line(x[0] + modx, x[1] + mody, (x[0] + modx + x[4]) / 2, abs(x[0] + modx - x[4]) / 4 + (x[1] + mody + x[5]) / 2, x[4], x[5],\
                                     smooth=True, width=3, tags=(newtag + 's', self.canvas.gettags(i)[1] , 'arc'))
        endtag = thetag + 'e'
        arcs= self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.create_line(x[0], x[1], (x[0] + x[4] + modx) / 2, abs(x[0] - x[4] - modx) / 4 + (x[1] + x[5] + mody) / 2, x[4] + modx, x[5] + mody,\
                                    smooth=True, width=3, tags=(self.canvas.gettags(i)[0], newtag + 'e', 'arc'))
        texttag = thetag + 't'
        try:
            textitem = self.canvas.find_withtag(texttag)[0]
            thetext = self.canvas.itemcget(textitem, "text")
            x = self.canvas.coords(textitem)
            self.canvas.create_text(x[0] + modx, x[1] + mody, fill='white', font=self.customFont,
                                    anchor=W, text=thetext, tags=(newtag, 'map', 'text', newtag + 't'))
        except IndexError:
            pass

    # add contig to the canvas
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
            if self.currspecial.get() == 'GC content':
                gccontent = self.contigDict[i].gccontent
                if gccontent <= self.gcmin.get():
                    colour = '#9933CC'
                elif gccontent >= self.gcmax.get():
                    colour = '#CC0000'
                else:
                    ratio = (gccontent - self.gcmin.get()) / (self.gcmax.get() - self.gcmin.get())
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
            elif self.currspecial.get() == 'GC skew':
                gcskew = self.contigDict[i].gcskew
                if gcskew <= self.skewmin.get():
                    colour = '#9933CC'
                elif gcskew >= self.skewmax.get():
                    colour = '#CC0000'
                else:
                    ratio = (gcskew - self.skewmin.get()) / (self.skewmax.get() - self.skewmin.get())
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
            elif self.currspecial.get() == 'AT skew':
                atskew = self.contigDict[i].atskew
                if atskew <= self.skewmin.get():
                    colour = '#9933CC'
                elif atskew >= self.skewmax.get():
                    colour = '#CC0000'
                else:
                    ratio = (atskew - self.skewmin.get()) / (self.skewmax.get() - self.skewmin.get())
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
            elif self.currspecial.get() == 'Coverage':
                coverage = self.contigDict[i].coverage
                if coverage <= self.covmin.get() or coverage == 'N/A':
                    colour = '#9933CC'
                elif coverage >= self.covmax.get():
                    colour = '#CC0000'
                else:
                    if self.covlog.get() == 0:
                        ratio = (coverage - self.covmin.get()) / (self.covmax.get() - self.covmin.get())
                    else:
                        ratio = (math.log10(coverage) - math.log10(self.covmin.get())) / (math.log10(self.covmax.get()) - math.log10(self.covmin.get()))
                    colour = '#%02x%02x%02x' % (int(ratio * 51) + 153, 51 - int(ratio * 51), 204 - int(ratio * 204))
            elif self.currspecial.get() == 'Custom':
                colour = self.contigDict[i].customcol
            elif self.currspecial.get() == 'None':
                colour = '#9933CC'
            self.canvas.create_rectangle(x1, y1, x2, y2, fill=colour, tags=('c' + i, 'contig', 'map'))
            if y2 > self.curryscroll:
                self.curryscroll = y2 + 20
                self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
            self.contigDict[i].orient[0] = True
            dir = '+'
            thetext = self.contigDict[i].shortname + ' ' + dir + ' ' + self.contigDict[i].strlen
            text = self.canvas.create_text(x1 + 2, y1 + self.contigheight/2, fill='white', font=self.customFont,
                                           anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
            if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                self.canvas.delete(text)
                thetext = self.contigDict[i].shortname + ' ' + dir
                text = self.canvas.create_text(x1 + 2, y1 + self.contigheight/2, fill='white', font=self.customFont,
                                               anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                    self.canvas.delete(text)
                    thetext = self.contigDict[i].shortname
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

    # bring up messsage asking which contig to be added
    def add_contig_dialogue(self):
        contig = tkSimpleDialog.askstring('Add Contig', 'Contig name')
        if contig is None:
            return
        if not contig in self.contigDict:
            for i in self.contigDict:
                if self.contigDict[i].shortname == contig:
                    contig = i
                    break
            if contig != i:
                tkMessageBox.showerror('Contig not found.', 'Please try again.')
                return
        if contig in self.visible: # if contig is already on canvas go to it
            self.goto(contig)
        else:
            self.contigDict[contig].visible = True
            self.contigDict[contig].xpos = 50
            self.contigDict[contig].ypos = 275
            self.contigDict[contig].xlength = self.contigDict[contig].length / self.newscaledown
            self.add_contig(contig)
            self.canvas.xview_moveto(0)
            self.canvas.yview_moveto(0)

    # add not yet visible contig adjacent to the 3' end of the selected contig
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
                    starty += self.contigheight + 10
        else:
            for i in self.contigDict[contig].fr:
                if not i[0] in self.visible:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                    starty += self.contigheight + 10

    # move contig or reference to coordinates x, y
    def move_contig(self, contig, x, y, ref=False):
        if ref:
            thetag = 'r' + contig
        else:
            thetag = 'c' + contig
        thecoords = self.canvas.coords(thetag)
        diffx = x - thecoords[0]
        diffy = y - thecoords[1]
        self.oldx = x
        self.oldy = y
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
            self.canvas.coords(i, x[0], x[1], x[2] + diffx/2, abs(x[0] - (x[4] + diffx)) / 4 + (x[1] + x[5] + diffy) / 2, x[4] + diffx, x[5] + diffy)
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
        starttag = thetag + 'ss'
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0] + diffx, x[1] + diffy, x[0] + diffx, x[1] + diffy, (x[0] + diffx + x[6]) / 2, abs(x[0] + diffx - x[6]) /4 + (x[1] + diffy + x[7]) / 2, x[6], x[7],
                                x[6], x[7], x[10], x[11], x[10], x[11], (x[16] + diffx + x[10]) / 2, abs(x[16] + diffx - x[10]) /4 + (x[1] + diffy + x[7]) / 2,
                                x[16] + diffx, x[17] + diffy, x[16] + diffx, x[17] + diffy)
        endtag = thetag + 'se'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0], x[1], x[0], x[1], (x[0] + x[6] + diffx) / 2, abs(x[0] - (x[6] + diffx)) /4 + (x[1] + x[7] + diffy) / 2, x[6] + diffx, x[7] + diffy,
                                x[6] + diffx, x[7] + diffy, x[10] + diffx, x[11] + diffy, x[10] + diffx, x[11] + diffy, (x[16] + x[10] + diffx) / 2,
                                abs(x[16] - (x[10] + diffx)) /4 + (x[1] + x[7] + diffy) / 2, x[16], x[17], x[16], x[17])

    # move all contigs adjacent to selected contigs 3' end
    def move_to(self):
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
                if i[0] in self.visible:
                    self.move_contig(i[0], startx, starty)
                    starty += self.contigheight + 10
        else:
            for i in self.contigDict[contig].fr:
                if i[0] in self.visible:
                    self.move_contig(i[0], startx, starty)
                    starty += self.contigheight + 10

    # duplicated all contigs adjacent to the 3' end of the selected contig, if contig is not already visible add it
    def dupe_to(self):
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
                if i[0] in self.visible:
                    self.duplicate_contig(i[0], (startx, starty))
                else:
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                starty += self.contigheight + 10
        else:
            for i in self.contigDict[contig].fr:
                if i[0] in self.visible:
                    self.duplicate_contig(i[0], (startx, starty))
                else:
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                starty += self.contigheight + 10

    # unhighlight edges
    def to_black(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='black')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='black')

    # highlight all 3' edges
    def to_red(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='red')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='red')

    def to_blue(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='blue')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='blue')

    def to_green(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='green')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='green')

    def to_yellow(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='green')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[2], 4):
                self.canvas.itemconfig(i, fill='green')

    # duplicates of above commands except for contigs adjacent to the 5' end of the contig n.b. 3' 5' relative to orientation of contig on canvas
    def show_fr(self):
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
                    self.add_contig(i[0])
                    starty += self.contigheight + 10
        else:
            for i in self.contigDict[contig].fr:
                if not i[0] in self.visible:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.contigDict[i[0]].xpos = startx - 10 - self.contigDict[i[0]].xlength
                    self.add_contig(i[0])
                    starty += self.contigheight + 10

    def move_fr(self):
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
        if self.contigDict[contig].orient[dupnum]:
            for i in self.contigDict[contig].fr:
                if i[0] in self.visible:
                    self.move_contig(i[0], startx  - 10 - self.contigDict[i[0]].xlength, starty)
                    starty += self.contigheight + 10
        else:
            for i in self.contigDict[contig].to:
                if i[0] in self.visible:
                    self.move_contig(i[0], startx  - 10 - self.contigDict[i[0]].xlength, starty)
                    starty += self.contigheight + 10

    def dupe_fr(self):
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
        if self.contigDict[contig].orient[dupnum]:
            for i in self.contigDict[contig].fr:
                if i[0] in self.visible:
                    self.duplicate_contig(i[0], (startx - 10 - self.contigDict[i[0]].xlength, starty))
                else:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                starty += self.contigheight + 10
        else:
            for i in self.contigDict[contig].to:
                if i[0] in self.visible:
                    self.duplicate_contig(i[0], (startx - 10 - self.contigDict[i[0]].xlength, starty))
                else:
                    self.contigDict[i[0]].visible = True
                    self.contigDict[i[0]].xpos = startx
                    self.contigDict[i[0]].ypos = starty
                    self.contigDict[i[0]].xlength = self.contigDict[i[0]].length / self.newscaledown
                    self.add_contig(i[0])
                starty += self.contigheight + 10

    def fr_black(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='black')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='black')


    def fr_red(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='red')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='red')

    def fr_blue(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='blue')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='blue')

    def fr_green(self):
        thetag = self.rctag
        starttag = thetag + 's'
        bbox = self.canvas.coords(thetag)
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[0], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='green')
        endtag = thetag + 'e'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            if round(x[4], 4) == round(bbox[0], 4):
                self.canvas.itemconfig(i, fill='green')

    # Ask for contig name then find on canvas
    def find_contig(self):
        contig = tkSimpleDialog.askstring('Find Contig', 'Contig name')
        if contig != None:
            if not contig in self.contigDict:
                for i in self.contigDict:
                    if self.contigDict[i].shortname == contig:
                        contig = i
                        break
                if contig != i:
                    tkMessageBox.showerror('Contig not found.', 'Please try again.')
                    return
            self.goto(contig)

    # go to chosen contig on the canvas
    def goto(self, contig, ref=False):
        try:
            if ref:
                x = self.canvas.coords('r' + contig)
            else:
                x = self.canvas.coords('c' + contig)
            self.canvas.xview_moveto(max([0, x[0] - 200]) / self.currxscroll)
            self.canvas.yview_moveto(max([0, x[1] - 200]) / self.curryscroll)
        except:
            pass

    # zoom in (from menu)
    def zoominmenu(self):
        self.canvas.scale(ALL, 0, 0, 1.1, 1.1)
        self.currxscroll *= 1.1
        self.curryscroll *= 1.1
        self.contigheight *= 1.1
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.fontsize *= 1.1
        self.customFont.configure(size=int(round(self.fontsize)))

    # zoom out (from menu)
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

    # zoom in (mouse wheel)
    def zoomin(self, event):
        self.canvas.scale(ALL, 0, 0, 1.1, 1.1)
        self.currxscroll *= 1.1
        self.curryscroll *= 1.1
        self.contigheight *= 1.1
        canvx = self.canvas.canvasx(event.x)
        canvy = self.canvas.canvasy(event.y)
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.canvas.xview_moveto((canvx * 1.1 - event.x) / self.currxscroll)
        self.canvas.yview_moveto((canvy * 1.1 - event.y) / self.curryscroll)
        self.fontsize *= 1.1
        self.customFont.configure(size=int(round(self.fontsize)))

    # zoom out (mouse wheel)
    def zoomout(self, event):
        if self.fontsize < 1 or self.currxscroll *0.909090909 < self.canvas.winfo_width():
            pass
        else:
            self.canvas.scale(ALL, 0, 0, 0.909090909, 0.909090909)
            self.currxscroll *= 0.909090909
            self.curryscroll *= 0.909090909
            self.contigheight *= 0.909090909
            canvx = self.canvas.canvasx(event.x)
            canvy = self.canvas.canvasy(event.y)
            self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
            self.canvas.xview_moveto((canvx * 0.909090909 - event.x) / self.currxscroll)
            self.canvas.yview_moveto((canvy * 0.909090909 - event.y) / self.curryscroll)
            self.fontsize *= 0.909090909
            self.customFont.configure(size=int(round(self.fontsize)))

    # shorten contigs
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
                contigend = self.canvas.coords(contigtag)[2]
                if textend >= contigend - 2:
                    self.canvas.delete(z)
                    thetext = self.contigDict[i].shortname + ' ' + dir
                    text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                   anchor=W, text=thetext, tags=(contigtag, 'map', 'text', contigtag + 't'))
                    textend = self.canvas.bbox(text)[2]
                    if textend >= contigend - 2:
                        self.canvas.delete(text)
                        thetext = self.contigDict[i].shortname
                        text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                       anchor=W, text=thetext, tags=(contigtag, 'map', 'text', contigtag + 't'))
                        textend = self.canvas.bbox(text)[2]
                        if textend >= contigend:
                            self.canvas.delete(text)

    # lengthen contigs
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
                contigend = self.canvas.coords(z)[2]
                thetext = self.contigDict[i].shortname + ' ' + dir + ' ' + length
                text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                anchor=W, text=thetext, tags=(contigtag, 'map', 'text', contigtag + 't'))
                textend = self.canvas.bbox(text)[2]
                if textend >= contigend - 2:
                    self.canvas.delete(text)
                    thetext = self.contigDict[i].shortname + ' ' + dir
                    text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                   anchor=W, text=thetext, tags=(contigtag, 'map', 'text', contigtag + 't'))
                    textend = self.canvas.bbox(text)[2]
                    if textend >= contigend - 2:
                        self.canvas.delete(text)
                        thetext = self.contigDict[i].shortname
                        text = self.canvas.create_text(textcoords[0], textcoords[1], fill='white', font=self.customFont,
                                                       anchor=W, text=thetext, tags=('c' + i, 'map', 'text', contigtag + 't'))
                        textend = self.canvas.bbox(text)[2]
                        if textend >= contigend - 2:
                            self.canvas.delete(text)

    # mouse wheel zoom, mousewheel event, only for windows? need to test
    def zoomcanvas(self, event):
        if event.delta == 120:
            self.zoomin(event)
        elif event.delta == -120:
            self.zoomout(event)

    # move canvas with arrows
    def scrollleft(self, event):
        self.canvas.xview_scroll(-2, 'units')

    def scrollright(self, event):
        self.canvas.xview_scroll(2, 'units')

    def scrollup(self, event):
        self.canvas.yview_scroll(-2, 'units')

    def scrolldown(self, event):
        self.canvas.yview_scroll(2, 'units')

    # for drag moving
    def recordMark(self, event):
        self.oldx = self.canvas.canvasx(event.x)
        self.oldy = self.canvas.canvasy(event.y)

    # drag contigs places
    def dragMove(self, event):
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        if y < self.contigheight/2:
            y= self.contigheight/2
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
            self.canvas.coords(i, x[0], x[1], x[2] + diffx/2, abs(x[0] - (x[4] + diffx)) / 4 + (x[1] + x[5] + diffy) / 2, x[4] + diffx, x[5] + diffy)
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
        starttag = thetag + 'ss'
        arcs = self.canvas.find_withtag(starttag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0] + diffx, x[1] + diffy, x[0] + diffx, x[1] + diffy, (x[0] + diffx + x[6]) / 2, abs(x[0] + diffx - x[6]) /4 + (x[1] + diffy + x[7]) / 2, x[6], x[7],
                                x[6], x[7], x[10], x[11], x[10], x[11], (x[16] + diffx + x[10]) / 2, abs(x[16] + diffx - x[10]) /4 + (x[1] + diffy + x[7]) / 2,
                                x[16] + diffx, x[17] + diffy, x[16] + diffx, x[17] + diffy)
        endtag = thetag + 'se'
        arcs = self.canvas.find_withtag(endtag)
        for i in arcs:
            x = self.canvas.coords(i)
            self.canvas.coords(i, x[0], x[1], x[0], x[1], (x[0] + x[6] + diffx) / 2, abs(x[0] - (x[6] + diffx)) /4 + (x[1] + x[7] + diffy) / 2, x[6] + diffx, x[7] + diffy,
                                x[6] + diffx, x[7] + diffy, x[10] + diffx, x[11] + diffy, x[10] + diffx, x[11] + diffy, (x[16] + x[10] + diffx) / 2,
                                abs(x[16] - (x[10] + diffx)) /4 + (x[1] + x[7] + diffy) / 2, x[16], x[17], x[16], x[17])

    # push text here to update the console
    def update_console(self, text):
        if self.consoletext.get() == '':
            self.consoletext.set(text)
        else:
            self.consoletext.set(self.consoletext.get().strip('.') + '.\n' + text) # strip excess dots once process has stopped running

    # the . .. ... . .. ... effect in the console
    def dot_console(self):
        text = self.consoletext.get()
        if text[-3:] == '...':
            self.consoletext.set(text[:-2])
        else:
            self.consoletext.set(text + '.')

    # reset canvas to default settings
    def clear_all(self):
        self.canvas.delete(ALL)
        self.edgelist = []
        self.contigDict = {}
        self.visible = set()
        self.currxscroll = self.originalxscroll
        self.curryscroll = self.originalyscroll
        self.contigheight = 25
        self.fontsize = 12
        self.customFont.configure(size=12)
        self.clear_lists()
        self.hitlist = None
        self.reforder = None
        self.rightmost = None
        self.leftmost = None
        self.blastfile.set('')

    # load an assembly from a variety of files
    def load_assembly(self):
        filename = tkFileDialog.askopenfilename(title='select assembly graph file')
        if filename == '' or filename == ():
            return
        self.clear_all()
        self.csagfile.set(filename)
        whatsthis = open(self.csagfile.get()) # figure out filetype from first line of file
        what = whatsthis.readline()
        whatsthis.close()
        self.contigDict = {}
        self.edgelist = []
        if what[0] == '>':
            self.contigfile.set(filename)
            gotit = self.load_fasta()
            if gotit == 1:
                self.update_console('FASTA loaded.')
            elif gotit == 2:
                self.update_console('.dot loaded.')
        elif what.split()[0] == 'NODE':
            self.load_csag()
            self.update_console('CAG loaded.')
        elif len(what.split()) == 4 and what.split()[0].isdigit() and what.split()[1].isdigit()\
            and what.split()[2].isdigit() and what.split()[3].isdigit():
            self.load_lastgraph()
            self.update_console('LastGraph loaded.')
        elif what[:2] == 'AS':
            self.load_ace()
            self.update_console('ACE file loaded.')
        elif what.startswith('digraph adj {'):
            self.contigfile.set(filename)
            gotit = self.load_fasta()
            if gotit == 2:
                self.update_console('.dot loaded.')
        elif what.startswith('#FASTG'):
            self.load_fastg()
            self.update_console('FASTG loaded.')
        elif what.startswith('H\tVN:Z:'):
            self.load_gfa()
            self.update_console('GFA loaded.')
        else:
            tkMessageBox.showerror('Invalid format', 'Contiguity cannot recognise file type.')
        self.writeWorkCont()

    def load_gfa(self):
        seqdict = {}
        with open(self.csagfile.get()) as gfa:
            for line in gfa:
                if line.startswith('S\t'):
                    sname = line.split()[1]
                    seq = line.split()[2]
                    if seq == '*' and seqdict == {}:
                        fasta = tkFileDialog.askopenfilename(title='select FASTA file.')
                        with open(fasta) as f:
                            first = True
                            for faline in f:
                                if faline.startswith('>'):
                                    if first:
                                        first = False
                                    else:
                                        seqdict[name] = seq
                                    name = faline[1:].split()[0]
                                    seq = ''
                                else:
                                    seq += faline.rstrip()
                            seqdict[name] = seq
                    if seq == '*':
                        seq = seqdict[sname]
                    aninstance = contig(sname, sname, seq)
                    self.contigDict[sname] = aninstance
                if line.startswith('L\t') or line.startswith('C\t'):
                    splitline = line.split()
                    if splitline[2] == '+':
                        dira = True
                    else:
                        dira = False
                    if splitline[4] == '+':
                        dirb = True
                    else:
                        dirb = False
                    cigar = splitline[5]
                    intstring = ''
                    overlap = 0
                    for i in cigar:
                        if i.isdigit():
                            intstring += i
                        else:
                            if i in ['M', 'D', 'N', 'H', 'P']:
                                overlap += int(intstring)
                                intstring = ''
                    self.edgelist.append((splitline[1], dira, splitline[3], dirb, overlap))
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


    # load ace file
    def load_ace(self):
        ace = open(self.csagfile.get())
        getseq = False
        edgedict = {}
        name = None
        for line in ace:
            if line.startswith('CO '):
                if not name is None:
                    coverage = totalbase * 1.0 / len(seq)
                    aninstance = contig(name, shortname, seq, None, coverage)
                    self.contigDict[name] = aninstance
                name = line.split()[1]
                shortname = str(int(name[6:]))
                totalbase = 0
                getseq = True
                seq = ''
            elif line.startswith('BQ'):
                getseq = False
            elif line.startswith('AF '):
                readname = line.split()[1]
                if '.to' in readname and '.fr' in readname:
                    pass
                elif '.to' in readname:
                    readname = readname.split('.')[0]
                    if readname in edgedict:
                        if not (name, True, edgedict[readname][0], edgedict[readname][1], 0) in self.edgelist and not \
                                    (edgedict[readname][0], not edgedict[readname][1], name, False, 0) in self.edgelist:
                            self.edgelist.append((name, True, edgedict[readname][0], edgedict[readname][1], 0))
                    else:
                        edgedict[readname] = [name, False]
                elif '.fm' in readname:
                    readname = readname.split('.')[0]
                    if readname in edgedict:
                        if not (name, False, edgedict[readname][0], edgedict[readname][1], 0) in self.edgelist and not \
                                    (edgedict[readname][0], not edgedict[readname][1], name, True, 0) in self.edgelist:
                            self.edgelist.append((name, False, edgedict[readname][0], edgedict[readname][1], 0))
                    else:
                        edgedict[readname] = [name, True]
            elif getseq:
                seq += line.rstrip().replace('*', '').upper()
            elif line.startswith('QA '):
                totalbase += int(line.split()[2]) - int(line.split()[1])
        coverage = totalbase * 1.0 / len(seq)
        aninstance = contig(name, shortname, seq, None, coverage)
        self.contigDict[name] = aninstance
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

    # load LastGRaph file
    def load_lastgraph(self):
        lg = open(self.csagfile.get())
        getseq = 0
        for line in lg:
            if line.startswith('NODE'):
                name = line.split()[1]
                getseq = 1
            elif getseq == 1:
                forseq = line.rstrip()
                getseq = 2
            elif getseq == 2:
                revseq = line.rstrip()
                aninstance = contig(name, name, forseq, revseq)
                self.contigDict[name] = aninstance
                getseq = 0
            elif line.startswith('ARC'):
                nodea, nodeb = line.split()[1:3]
                if nodea[0] == '-':
                    nodea = nodea[1:]
                    dira = False
                else:
                    dira = True
                if nodeb[0] == '-':
                    nodeb = nodeb[1:]
                    dirb = False
                else:
                    dirb = True
                self.edgelist.append((nodea, dira, nodeb, dirb, 0))
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


    def load_fastg(self):
        fg = open(self.csagfile.get())
        name = None
        for line in fg:
            if line.startswith('#'):
                pass
            elif line.startswith('>'):
                if not name is None:
                    if '[' in seq:
                        newseq = ''
                        for i in seq.split('['):
                            newseq += i.split[']'][-1]
                        seq = newseq
                    aninstance = contig(entry, name, seq)
                    if name in self.contigDict:
                        tkMessageBox.showerror('FASTG error', 'Please make sure the section of the FASTA header preceeding the first whitespace is unique.')
                        return None
                    else:
                        self.contigDict[entry] = aninstance
                if ':' in line:
                    splitline = line.rstrip().split(':')
                    entry = splitline[0][1:]
                    for i in splitline[1:]:
                        if i[0] == '~':
                            if i[-1] == '\'':
                                self.edgelist.append((entry, False, i[1:-1], False, 0))
                            else:
                                self.edgelist.append((entry, False, i[1:], True, 0))
                        else:
                            if i[-1] == '\'':
                                self.edgelist.append((entry, True, i[:-1], False, 0))
                            else:
                                self.edgelist.append((entry, True, i, True, 0))
                else:
                    entry = line.rstrip()[1:-1]

                if entry.startwith('NODE_'):
                    name = entry.split('_')[1]
                else:
                    name = entry
                seq = ''
            else:
                seq += line.rstrip()
        if '[' in seq:
            newseq = ''
            for i in seq.split('['):
                newseq += i.split[']'][-1]
            seq = newseq
        aninstance = contig(entry, name, seq)
        if name in self.contigDict:
            tkMessageBox.showerror('FASTG error', 'Please make sure the section of the FASTA header preceeding the first whitespace is unique.')
            return None
        else:
            self.contigDict[entry] = aninstance
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




    # Load FASTA file
    def load_fasta(self):
        fastafile = open(self.contigfile.get())
        namelist = []
        maxlen = 0
        for line in fastafile:
            if line.startswith('>'):
                namelist.append(line.split()[0])
                if maxlen < len(namelist[-1]):
                    maxlen = len(namelist[-1])
        if namelist == []:
            tkMessageBox.showerror('.dot error', 'Please make sure the FASTA sequence of this file is concatenated to the end.')
            return None
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
                cuta = 1
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
        fastafile = open(self.contigfile.get())
        first = True
        getgraph = 0
        isgraph = False
        for line in fastafile:
            if line.startswith('digraph adj {'):
                getgraph = 1
                isgraph = True
            elif getgraph == 1:
                getgraph = 2
            elif getgraph == 2:
                getgraph =3
            elif getgraph == 3:
                if line.startswith('}'):
                    getgraph = 0
                else:
                    if len(line.split('"')) > 3:
                        nodea = line.split('"')[1]
                        nodeb = line.split('"')[3]
                        if nodea[-1] == '-':
                            dira = False
                        else:
                            dira = True
                        nodea = nodea[:-1]
                        if nodeb[-1] == '-':
                            dirb = False
                        else:
                            dirb = True
                        nodeb = nodeb[:-1]
                        self.edgelist.append((nodea, dira, nodeb, dirb, 0))
            elif line.startswith('>'):
                if first:
                    first = False
                else:
                    aninstance = contig(entry, name, seq)
                    if entry in self.contigDict:
                        tkMessageBox.showerror('FASTA error', 'Please make sure the section of the FASTA header preceeding the first whitespace is unique.')
                        return None
                    else:
                        self.contigDict[entry] = aninstance
                entry = line.split()[0][1:]
                if entry.startswith('NODE_'):
                    name = entry.split('_')[1]
                else:
                    name = line.split()[0][cuta:cutb]
                    if name == '':
                        name = entry
                seq = ''
            else:
                seq += line.rstrip()
        aninstance = contig(entry, name, seq)
        self.contigDict[entry] = aninstance
        fastafile.close()
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
        if isgraph:
            return 2
        else:
            return 1


    # Load contiguity's native filetype
    def load_csag(self):
        csag = open(self.csagfile.get())
        for line in csag:
            if line.split()[0] == 'NODE':
                if len(line.split()) == 4:
                    title, entry, name, seq = line.split()
                    aninstance = contig(entry, name, seq)
                else:
                    title, entry, name, coverage, seq = line.split()
                    if coverage == 'N/A':
                        coverage = None
                    else:
                        coverage = float(coverage)
                    aninstance = contig(entry, name, seq, None, coverage)
                self.contigDict[entry] = aninstance
            elif line.split()[0] == 'EDGE':
                title, n1, d1, n2, d2, overlap = line.split()
                if overlap == '.':
                    overlap = ''
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

    # file dialogues, ask for filenames then set approriate variable
    def loadref(self):
        filename = tkFileDialog.askopenfilename(parent=self.blast_options)
        if filename == '':
            return
        self.reffile.set(filename)

    def loadbo(self):
        filename = tkFileDialog.askopenfilename(parent=self.write_blast_top)
        if filename == '':
            return
        self.bofile.set(filename)

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
        filename = tkFileDialog.askopenfilename(parent=self.create_edges_top)
        if filename == '':
            return
        self.readfile.set(filename)

    def loadwork(self):
        filename = tkFileDialog.askdirectory(title='Please select a working directory', initialdir=self.workingDir.get())
        if filename == '':
            sys.exit()
        self.workingDir.set(filename)
        try:
            os.makedirs(filename)
        except:
            pass

    def loadoutfile(self):
        filename = tkFileDialog.asksaveasfilename(parent=self.write_fasta)
        if filename == '':
            return
        else:
            self.outfile.set(filename)

    # check to see if thread should be aborted
    def aborttime(self):
        try:
            abort = self.abortqueue.get(0)
            if abort:
                return True
            else:
                return False
        except Queue.Empty:
            return False

    # colour options window
    def colour_options(self):
        try:
            self.colour_options_top.destroy()
        except AttributeError:
            pass
        self.colour_options_top = Toplevel()
        self.colour_options_top.grab_set()
        self.colour_options_top.wm_attributes("-topmost", 1)
        x = root.winfo_rootx() + 15
        y = root.winfo_rooty() - 20
        self.colour_options_top.title('Colour options')
        geom = "+%d+%d" % (x,y)
        self.colour_options_top.geometry(geom)
        self.colour_options_frame = Frame(self.colour_options_top)
        self.gcminlabel = Label(self.colour_options_frame, text='GC content min:')
        self.gcminlabel.grid(row=0, column=0, sticky=E)
        self.gcminentry = Entry(self.colour_options_frame, textvariable=self.gcmin)
        self.gcminentry.grid(row=0, column=1)
        self.gcmaxlabel = Label(self.colour_options_frame, text='GC content max:')
        self.gcmaxlabel.grid(row=1, column=0, sticky=E)
        self.gcmaxentry = Entry(self.colour_options_frame, textvariable=self.gcmax)
        self.gcmaxentry.grid(row=1, column=1)
        self.skewminlabel = Label(self.colour_options_frame, text='Skew min:')
        self.skewminlabel.grid(row=2, column=0, sticky=E)
        self.skewminentry = Entry(self.colour_options_frame, textvariable=self.skewmin)
        self.skewminentry.grid(row=2, column=1)
        self.skewmaxlabel = Label(self.colour_options_frame, text='Skew max:')
        self.skewmaxlabel.grid(row=3, column=0, sticky=E)
        self.skewmaxentry = Entry(self.colour_options_frame, textvariable=self.skewmax)
        self.skewmaxentry.grid(row=3, column=1)
        self.covminlabel = Label(self.colour_options_frame, text='Coverage min:')
        self.covminlabel.grid(row=4, column=0, sticky=E)
        self.covminentry = Entry(self.colour_options_frame, textvariable=self.covmin)
        self.covminentry.grid(row=4, column=1)
        self.covmaxlabel = Label(self.colour_options_frame, text='Coverage max:')
        self.covmaxlabel.grid(row=5, column=0, sticky=E)
        self.covmaxentry = Entry(self.colour_options_frame, textvariable=self.covmax)
        self.covmaxentry.grid(row=5, column=1)
        self.covloglabel = Label(self.colour_options_frame, text='Log coverage:')
        self.covloglabel.grid(row=6, column=0, sticky=E)
        self.covlogentry = Checkbutton(self.colour_options_frame, variable=self.covlog)
        self.covlogentry.grid(row=6, column=1)
        self.okcolour = Button(self.colour_options_frame, text='Ok', command=self.ok_colour)
        self.okcolour.grid(row=7, column=2, sticky=E)
        self.colour_options_frame.grid(padx=10, pady=10)

    def ok_colour(self):
        self.colour_options_top.destroy()
        self.changespecial(self.currspecial.get())

    def create_edges(self):
        try:
            self.create_edges_top.destroy()
        except:
            pass
        self.create_edges_top = Toplevel()
        self.create_edges_top.grab_set()
        self.create_edges_top.wm_attributes("-topmost", 1)
        self.create_edges_top.geometry('+20+30')
        self.create_edges_top.title('Create CAG')
        self.create_edges_frame = Frame(self.create_edges_top)
        self.contigfilelabel = Label(self.create_edges_frame, text='Contig file:')
        self.contigfilelabel.grid(row=0, column=0, sticky=E)
        self.contigfileentry = Entry(self.create_edges_frame, textvariable=self.contigfile, justify=RIGHT, width=30)
        self.contigfileentry.grid(row=0, column=1)
        self.contigfileentrybutton = Button(self.create_edges_frame, text='...', command=self.loadcontig)
        self.contigfileentrybutton.grid(row=0, column=2)
        self.readfilelabel = Label(self.create_edges_frame, text='Read file:')
        self.readfilelabel.grid(row=1, column=0, sticky=E)
        self.readfileentry = Entry(self.create_edges_frame, textvariable=self.readfile, justify=RIGHT, width=30)
        self.readfileentry.grid(row=1, column=1)
        self.readfileentrybutton = Button(self.create_edges_frame, text='...', command=self.loadread)
        self.readfileentrybutton.grid(row=1, column=2)
        self.create_edges_frame.grid(padx=10, pady=10)
        self.overlapEdges = Frame(self.create_edges_frame, relief=SUNKEN, bd=2)
        self.overlapEdges.grid(row=2, column=0, padx=5, pady=3, sticky=EW, columnspan=3)
        self.getOverlapLabel = Label(self.overlapEdges, text='Get overlapping edges:', width=30, anchor=E)
        self.getOverlapLabel.grid(row=0, column=0)
        self.getOverlapCheck = Checkbutton(self.overlapEdges, variable=self.getOverlap)
        self.getOverlapCheck.grid(row=0, column=1)
        self.minoverlaplabel = Label(self.overlapEdges, text='Minimum overlap length:', width=30, anchor=E)
        self.minoverlaplabel.grid(row=1, column=0)
        self.minoverlapentry = Entry(self.overlapEdges, textvariable=self.minoverlap)
        self.minoverlapentry.grid(row=1, column=1)
        self.maxmmlabel = Label(self.overlapEdges, text='Maximum mismatches in overlap:', width=30, anchor=E)
        self.maxmmlabel.grid(row=2, column=0)
        self.maxmmentry = Entry(self.overlapEdges, textvariable=self.maxmm)
        self.maxmmentry.grid(row=2, column=1, sticky=W)
        self.dbEdges = Frame(self.create_edges_frame, relief=SUNKEN, bd=2)
        self.dbEdges.grid(row=3, column=0, padx=5, pady=3, sticky=EW, columnspan=3)
        self.getdblabel = Label(self.dbEdges, text='Get de bruijn edges:', width=30, anchor=E)
        self.getdblabel.grid(row=0, column=0)
        self.getdbentry = Checkbutton(self.dbEdges, variable=self.getdb)
        self.getdbentry.grid(row=0, column=1)
        self.nmersizelabel = Label(self.dbEdges, text='kmer size:', width=30, anchor=E)
        self.nmersizelabel.grid(row=1, column=0)
        self.nmersizeentry = Entry(self.dbEdges, textvariable=self.nmersize)
        self.nmersizeentry.grid(row=1, column=1)
        self.maxdistlabel = Label(self.dbEdges, text='Max. distance:', width=30, anchor=E)
        self.maxdistlabel.grid(row=2, column=0)
        self.maxdistentry = Entry(self.dbEdges, textvariable=self.maxdist)
        self.maxdistentry.grid(row=2, column=1)
        self.cutautolabel = Label(self.dbEdges, text='Auto detect cutoffs:', width=30, anchor=E)
        self.cutautolabel.grid(row=3, column=0)
        self.cutautoentry = Checkbutton(self.dbEdges, variable=self.cutauto)
        self.cutautoentry.grid(row=3, column=1)
        self.nmercutlabel = Label(self.dbEdges, text='kmer cutoff:', width=30, anchor=E)
        self.nmercutlabel.grid(row=4, column=0)
        self.nmercutentry = Entry(self.dbEdges, textvariable=self.nmercut)
        self.nmercutentry.grid(row=4, column=1)
        self.nmeravelabel = Label(self.dbEdges, text='kmer average:', width=30, anchor=E)
        self.nmeravelabel.grid(row=5, column=0)
        self.nmeraveentry = Entry(self.dbEdges, textvariable = self.nmerave)
        self.nmeraveentry.grid(row=5, column=1)
        self.pairedEdges = Frame(self.create_edges_frame, relief=SUNKEN, bd=2)
        self.pairedEdges.grid(row=4, column=0, padx=5, pady=3, sticky=NSEW, columnspan=3)
        self.getPairedlabel = Label(self.pairedEdges, text='Get edges using paired reads:', width=30, anchor=E)
        self.getPairedlabel.grid(row=0, column=0)
        self.getPairedentry = Checkbutton(self.pairedEdges, variable=self.getPaired)
        self.getPairedentry.grid(row=0, column=1)
        self.insertsizelabel = Label(self.pairedEdges, text='Max. insert size:', width=30, anchor=E)
        self.insertsizelabel.grid(row=1, column=0)
        self.insertsizeentry = Entry(self.pairedEdges, textvariable=self.insertsize)
        self.insertsizeentry.grid(row=1, column=1)
        #self.orientationlabel = Label(self.pairedEdges, text='Read orientation:')
        #self.orientationlabel.grid(row=2, column=0)
        #self.orientationentry = OptionMenu(self.pairedEdges, self.orientation, '--> <--', '<-- -->', '--> -->')
        #self.orientationentry.grid(row=2, column=1, sticky=W)
        self.minrllabel = Label(self.pairedEdges, text='Minimum read length:', width=30, anchor=E)
        self.minrllabel.grid(row=3, column=0,)
        self.minrlentry = Entry(self.pairedEdges, textvariable=self.minrl)
        self.minrlentry.grid(row=3, column=1)
        self.minpairedgelabel = Label(self.pairedEdges, text='Min. reads for edge:', width=30, anchor=E)
        self.minpairedgelabel.grid(row=4, column=0)
        self.minpairedgeentry = Entry(self.pairedEdges, textvariable=self.minpairedge)
        self.minpairedgeentry.grid(row=4, column=1)
        #self.longEdges = Frame(self.create_edges_frame, relief=SUNKEN, bd=2)
        #self.longEdges.grid(row=7, column=0, padx=5, pady=3, sticky=NSEW, columnspan=3)
        #self.getLonglabel = Label(self.longEdges, text='Get edges using long reads:')
        #self.getLonglabel.grid(row=0, column=0)
        #self.getLongentry = Checkbutton(self.longEdges, variable=self.getLong)
        #self.getLongentry.grid(row=0, column=1, sticky=W)
        #self.minlonglabel = Label(self.longEdges, text='Min. read number:')
        #self.minlonglabel.grid(row=1, column=0)
        #self.minlongentry = Entry(self.longEdges, textvariable=self.minlong)
        #self.minlongentry.grid(row=1, column=1, sticky=W)
        #self.minlongoverlabel = Label(self.longEdges, text='Min. overlap:')
        #self.minlongoverlabel.grid(row=2, column=0)
        #self.minlongoverentry = Entry(self.longEdges, textvariable=self.minlongover)
        #self.minlongoverentry.grid(row=2, column=1, sticky=W)
        #self.minlongidentlabel = Label(self.longEdges, text='Min. Identity:')
        #self.minlongidentlabel.grid(row=3, column=0)
        #self.minlongidententry = Entry(self.longEdges, textvariable=self.minlongident)
        #self.minlongidententry.grid(row=3, column=1, sticky=W)
        self.okedges = Button(self.create_edges_frame, text='Ok', command=self.ok_edges)
        self.okedges.grid(row=8, column=2, sticky=E)

    def ok_edges(self):
        if (not os.path.exists(self.readfile.get()) and (self.getdb.get() == 1 or self.getPaired.get() == 1 or self.getLong.get() == 1))\
          or not os.path.exists(self.contigfile.get()):
            tkMessageBox.showerror('File missing', 'Please include a contig and read file.', parent=self.create_edges_top)
            return
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        try:
            global khmer
            import khmer
        except ImportError:
            proceed_no_khmer = tkMessageBox.askyesno('Khmer not found.', 'Proceed without installing Khmer? (not recommended - see manual)?')
            if proceed_no_khmer:
                args.khmer = False
            else:
                tkMessageBox.showerror('Aborted', 'Please install Khmer')
                return
        self.create_edges_top.destroy()
        self.edgeconsole = []
        self.lasteclen = 0
        self.queue = Queue.Queue()
        self.thethread = threading.Thread(target=self.edge_thread)
        self.thethread.start()
        self.abort = False
        self.update_console('Finding edges now.')
        self.update_edges()

    def update_edges(self):
        self.dot_console()
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text != 'Creating CAG finished.':
                    root.after(1000, self.update_edges)
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            if self.abort:
                self.update_console('Thread stopped.')
            else:
                self.update_console('Edge creation failed, please check console output.')
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_edges)
            return

    def edge_thread(self):
        self.edgelist = []
        self.contigDict = {}
        self.load_fasta()
        self.writeWorkCont()
        if self.getOverlap.get() == 1:
            self.queue.put('Finding overlaps.')
            count = len(self.edgelist)
            self.get_overlap_edges()
            count = len(self.edgelist) - count
            self.queue.put('Overlaps found ' + str(count) + ' edges.')
        if self.aborttime():
            return
        if self.getdb.get() == 1:
            self.queue.put('Finding De Bruijn edges.')
            count = len(self.edgelist)
            if args.khmer:
                abort = self.get_nmer_freq_khmer()
                if self.aborttime() or self.abort:
                    return
                if self.cutauto.get() == 1:
                    self.getnmercutkhmer()
                    if self.aborttime():
                        return
                    self.queue.put('Nmer cutoff A set to ' + str(self.nmercut.get()))
                    self.queue.put('Nmer cutoff B set to ' + str(self.nmerave.get()))
                abort = self.get_db_edges()
                if abort:
                    return
            else:
                abort = self.get_nmer_freq()
                if self.aborttime() or self.abort:
                    return
                self.read_nmer_freq()
                if self.aborttime():
                    return
                if self.cutauto.get() == 1:
                    self.getnmercut()
                    if self.aborttime():
                        return
                    self.queue.put('Nmer cutoff A set to ' + str(self.nmercut.get()))
                    self.queue.put('Nmer cutoff B set to ' + str(self.nmerave.get()))
                abort = self.get_db_edges()
                del self.ht
                if abort:
                    return
            count = len(self.edgelist) - count
            self.queue.put('De Bruijn found ' + str(count) + ' edges.')
        if self.aborttime():
            return
        if self.getPaired.get() == 1:
            self.queue.put('Finding paired end edges.')
            count = len(self.edgelist)
            abort = self.get_paired_edge()
            if abort:
                return
            count = len(self.edgelist) - count
            self.queue.put('Paired end found ' + str(count) + ' edges.')
        if self.aborttime():
            return
        # if self.getLong.get() == 1:
        #     self.queue.put('Finding edges with long reads.')
        #     count = len(self.edgelist)
        #     self.get_long_edge()
        #     count = len(self.edgelist) - count
        #     self.queue.put('Long reads found '+ str(count) + 'edges.')
        if self.aborttime():
            return
        self.removeDuplicates()
        if self.aborttime():
            return
        self.untangleEdges()
        if self.aborttime():
            return
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
        self.queue.put('Creating CAG finished.')


    def get_overlap_edges(self):
        minoverlap, maxmm = self.minoverlap.get(), self.maxmm.get()
        stderr = open(self.workingDir.get() + '/bwaerr.txt', 'w')
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

    def dbpathkhmer(self, currpath, endnmer):
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
            thepath = todo.pop().upper()
            search = True
            while search:
                if len(thepath) > tempmaxdist + 20:
                    search = False
                if len(thepath) > nmersize + 20 and thepath[-nmersize - 20:] in endnmer:
                    paths.append(thepath)
                if search:
                    newn = thepath[-nmersize + 1:] + 'A'
                    acount = self.ht.get(newn)
                    newn = thepath[-nmersize + 1:] + 'T'
                    tcount = self.ht.get(newn)
                    newn = thepath[-nmersize + 1:] + 'C'
                    ccount = self.ht.get(newn)
                    newn = thepath[-nmersize + 1:] + 'G'
                    gcount = self.ht.get(newn)
                    if max([acount, tcount, ccount, gcount]) < nmercut:
                        search = False
                    else:
                        maxnucl = max([(acount, 'A'), (tcount, 'T'), (ccount, 'C'), (gcount, 'G')])[1]
                        bases = 'ATCG'.replace(maxnucl, '')
                        for i in bases:
                            newn = thepath[-nmersize + 1:] + i
                            if self.ht.get(newn) > nmerave:
                                todo.append(thepath + i)
                        thepath += maxnucl
        countdict = {}
        for i in paths:
            total = 0
            acount = 0
            for j in range(0, len(i) - nmersize):
                newn = i[j:j+nmersize]
                total += self.ht.get(newn)
                acount += 1
            if i[-nmersize - 20:] in countdict:
                if total * 1.0 / acount > countdict[i[-nmersize - 20:]][1]:
                    countdict[i[-nmersize - 20:]] = (i, total * 1.0 / acount)
            else:
                countdict[i[-nmersize - 20:]] = (i, total * 1.0 / acount)
        paths = []
        for i in countdict:
            paths.append(countdict[i][0])
        return paths



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
                    newn = thepath[-nmersize + 1:] + 'A'
                    if newn[nmersize/2] == 'A' or newn[nmersize/2] == 'C':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        acount = self.nmerdict[newn]
                    else:
                        acount = 0
                    newn = thepath[-nmersize + 1:] + 'T'
                    if newn[nmersize/2] == 'A' or newn[nmersize/2] == 'C':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        tcount = self.nmerdict[newn]
                    else:
                        tcount = 0
                    newn = thepath[-nmersize + 1:] + 'C'
                    if newn[nmersize/2] == 'A' or newn[nmersize/2] == 'C':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        ccount = self.nmerdict[newn]
                    else:
                        ccount = 0
                    newn = thepath[-nmersize + 1:] + 'G'
                    if newn[nmersize/2] == 'A' or newn[nmersize/2] == 'C':
                        newn = newn[::-1]
                        newn = newn.translate(transtab)
                    if newn in self.nmerdict:
                        gcount = self.nmerdict[newn]
                    else:
                        gcount = 0
                    if max([acount, tcount, ccount, gcount]) < nmercut:
                        search = False
                    else:
                        maxnucl = max([(acount, 'A'), (tcount, 'T'), (ccount, 'C'), (gcount, 'G')])[1]
                        bases = 'ATGC'.replace(maxnucl, '')
                        for i in bases:
                            newn = thepath[-nmersize + 1:] + i
                            if newn[nmersize/2] == 'A' or newn[nmersize/2] == 'C':
                                newn = newn[::-1]
                                newn = newn.translate(transtab)
                            if newn in self.nmerdict and self.nmerdict[newn] > nmerave:
                                todo.append(thepath + i)
                        thepath += maxnucl
        countdict = {}
        for i in paths:
            acount = 0
            total = 0
            for j in range(0, len(i) - nmersize):
                newn = i[j:j+nmersize]
                if newn[nmersize/2] == 'A' or newn[nmersize/2] == 'C':
                    newn = newn[::-1]
                    newn = newn.translate(transtab)
                acount += 1
                try:
                    total += self.nmerdict[newn]
                except:
                    pass
            if i[-nmersize - 20:] in countdict:
                if total * 1.0 / acount > countdict[i[-nmersize - 20:]][1]:
                    countdict[i[-nmersize - 20:]] = (i, total * 1.0 / acount)
            else:
                countdict[i[-nmersize - 20:]] = (i, total * 1.0 / acount)
        paths = []
        for i in countdict:
            paths.append(countdict[i][0])
        return paths

    def getPerms(self, nmer):
        outlist = set()
        for i in range(len(nmer)):
            outlist.add(nmer[:i] + 'A' + nmer[i+1:])
            outlist.add(nmer[:i] + 'T' + nmer[i+1:])
            outlist.add(nmer[:i] + 'C' + nmer[i+1:])
            outlist.add(nmer[:i] + 'G' + nmer[i+1:])
        return(list(outlist))

    def get_db_edges(self):
        nmersize = self.nmersize.get()
        ends = {}
        for i in self.contigDict:
            name = i
            seq = self.contigDict[i].forseq
            nmer = seq[:nmersize + 20]
            #nmerlist = self.getPerms(nmer)
            #for nmer in nmerlist:
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
            if self.aborttime():
                return True
            if args.khmer:
                forpaths = self.dbpathkhmer(self.contigDict[i].forseq[-nmersize:], ends)
                self.contigDict[i].coverage = self.ht.get_median_count(self.contigDict[i].forseq.upper())[0]
            else:
                forpaths = self.dbpath(self.contigDict[i].forseq[-nmersize:], ends)
                ncount = []
                for zz in range(0, len(self.contigDict[i].forseq) - nmersize):
                    nmer = self.contigDict[i].forseq[zz:zz+nmersize]
                    if nmer[nmersize/2] == 'A' or nmer[nmersize/2] == 'C':
                        nmer = nmer[::-1]
                        nmer = nmer.translate(transtab)
                    if nmer in self.nmerdict:
                        ncount.append(self.nmerdict[nmer])
                    else:
                        ncount.append(0)
                ncount.sort()
                self.contigDict[i].coverage = ncount[len(ncount)/2]
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

            if args.khmer:
                revpaths = self.dbpathkhmer(self.contigDict[i].revseq[-nmersize:], ends)
            else:
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
        if args.khmer:
            del self.ht
        else:
            del self.nmerdict
        return False

    def getflag(self, n, count=11):
        return "".join([str((int(n) >> y) & 1) for y in range(count-1, -1, -1)])

    def get_paired_edge(self):
        minoverlap = self.minoverlap.get()
        if self.readfile.get()[-3:] == '.gz':
            tempread = gzip.open(self.readfile.get())
        else:
            tempread = open(self.readfile.get())
        firstline = tempread.readline()
        if firstline[0] == '@':
            readType = 'fq'
        elif firstline[0] == '>':
            readType = 'fa'
        else:
            try:
                self.edgeconsole.append('Read file not valid type, skipping paired edge creation.')
            except NameError:
                sys.stderr.write('Read file not valid type, skipping paired end creation.')
            tempread.close()
            return
        tempread.close()
        insertsize, minpairedge, minrl = self.insertsize.get(), self.minpairedge.get(), self.minrl.get()
        thestderr = open(self.workingDir.get() + '/bwaerr.txt', 'w')
        subprocess.Popen('bowtie2-build -f ' + self.workingDir.get() + '/contigs.fa ' + self.workingDir.get() + '/fullcontigs', shell=True, stderr=thestderr, stdout=thestderr).wait()
        if readType == 'fa':
            subprocess.Popen('bowtie2 --local -x ' + self.workingDir.get() + '/fullcontigs -U ' + self.readfile.get() + ' -S ' + self.workingDir.get() + '/aln2.sam -f -a -p ' + str(self.num_threads.get()), shell=True, stderr=thestderr, stdout=thestderr).wait()
        elif readType == 'fq':
            subprocess.Popen('bowtie2 --local -x ' + self.workingDir.get() + '/fullcontigs -U ' + self.readfile.get() + ' -S ' + self.workingDir.get() + '/aln2.sam -q -a -p ' + str(self.num_threads.get()), shell=True, stderr=thestderr, stdout=thestderr).wait()
        thestderr.close()
        samfile = open(self.workingDir.get() + '/aln2.sam')
        zereadpair = False
        count = 0
        edgedict = {}
        for line in samfile:
            if not line.startswith('@'):
                if not self.abortqueue.empty():
                    if self.aborttime():
                        return True
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
        return False

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

    def get_nmer_freq_khmer(self):
        nmersize, reads, ht_size, ht_n = self.nmersize.get(), self.readfile.get(), self.ht_size.get(), self.ht_number.get()
        ht_size = float('2e9')
        ht_n = 4
        bigcount = True
        self.ht = khmer.new_counting_hash(nmersize, ht_size, ht_n) # HT_size, number ht, threads
        self.ht.set_use_bigcount(bigcount)
        rparser = khmer.ReadParser(reads)
        self.queue.put('consuming input ' + reads)
        self.ht.consume_fasta_with_reads_parser(rparser)
        fp_rate = khmer.calc_expected_collisions(self.ht)
        self.queue.put('fp rate estimated to be %1.3f' % fp_rate)
        if fp_rate > 0.20:
            self.queue.put("**\n** ERROR: the counting hash is too small for\n** this data set.  Increase hashsize/num ht.\n**")
            sys.exit(1)

    def get_nmer_freq(self):
        nmersize, reads  = self.nmersize.get(), self.readfile.get()
        nucl = set('ATGC')
        self.nmerdict = {}
        if reads[-3:] == '.gz':
            readfile = gzip.open(reads)
        else:
            readfile = open(reads)
        getseq = False
        getfaseq = False
        for line in readfile:
            if not self.abortqueue.empty():
                if self.aborttime():
                    return True
            if line.startswith('@'):
                getseq = True
            elif line.startswith('>'):
                if getfaseq:
                    for i in range(0, len(seq) - nmersize):
                        nmer = seq[i:i+nmersize]
                        if nmer[nmersize/2] == 'A' or nmer[nmersize/2] == 'C':
                            nmer = nmer[::-1]
                            nmer = nmer.translate(transtab)
                        if nmer in self.nmerdict:
                            self.nmerdict[nmer] += 1
                        else:
                            if set(nmer).issubset(nucl):
                                self.nmerdict[nmer] = 1
                getfaseq = True
                seq = ''
            elif getfaseq:
                seq += line.rstrip().upper()
            elif getseq:
                getseq = False
                seq = line.rstrip().upper()
                for i in range(0, len(seq) - nmersize + 1):
                    nmer = seq[i:i+nmersize]
                    if nmer[nmersize/2] == 'A' or nmer[nmersize/2] == 'C':
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
                if nmer[nmersize/2] == 'A' or nmer[nmersize/2] == 'C':
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
        return False # thread has not been aborted


    def read_nmer_freq(self):
        self.nmerdict = {}
        nmercut = self.nmercut.get()
        nmers = open(self.nmerfile)
        for line in nmers:
            nmer, freq = line.split()
            freq = int(freq)
            if freq > nmercut:
                self.nmerdict[nmer] = freq

    def getnmercutkhmer(self):
        freqDict = {}
        for i in self.contigDict:
            for j in range(0, len(self.contigDict[i].forseq.upper()) - self.nmersize.get()):
                nmer = self.contigDict[i].forseq[j:j+self.nmersize.get()].upper()
                count = self.ht.get(nmer)
                if count in freqDict:
                    freqDict[count] += 1
                else:
                    freqDict[count] = 1
        tf = None
        maxcount = 0
        total = 0
        for i in range(max(freqDict)):
            if i in freqDict:
                if freqDict[i] > maxcount:
                    maxcount = freqDict[i]
                    mf = i
                total += freqDict[i]
        runningtotal = 0
        for i in range(max(freqDict)):
            if i in freqDict:
                runningtotal += freqDict[i]
                if runningtotal >= total / 10000:
                    tf = i
                    break
        if self.nmercut.get() == -1:
            self.nmercut.set(tf)
        if self.nmerave.get() == -1:
            self.nmerave.set(mf/2)


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
        if self.nmercut.get() == -1 or self.nmerave.get() == -1:
            if self.nmercut.get() == -1:
                self.nmercut.set(tf)
            if self.nmerave.get() == -1:
                self.nmerave.set(int(mf)/2)
        else:
            self.nmerave.set(int(mf)/2)
            self.nmercut.set(tf)

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
        out = open(self.workingDir.get() + '/assembly.cag', 'w')
        for i in self.contigDict:
            out.write('NODE\t' + i + '\t' + self.contigDict[i].shortname + '\t' + str(self.contigDict[i].coverage) + '\t' + self.contigDict[i].forseq + '\n')
        for i in self.edgelist:
            if i[4] == '':
                overlap = '.'
            else:
                overlap = str(i[4])
            out.write('EDGE\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + str(i[3]) + '\t' + overlap + '\n')
        out.close()

    def create_comp(self):
        try:
            self.blast_options.destroy()
        except:
            pass
        self.blast_options = Toplevel()
        self.blast_options.grab_set()
        self.blast_options.wm_attributes("-topmost", 1)
        self.frame2 = Frame(self.blast_options)
        self.blast_options.geometry('+20+30')
        self.blast_options.title('Create comparison')
        self.reffilelabel = Label(self.frame2, text='Reference file:', anchor=E)
        self.reffilelabel.grid(column=0, row=1, sticky=E)
        self.reffileentry = Entry(self.frame2, textvariable=self.reffile)
        self.reffileentry.grid(column=1, row=1)
        self.reffilebutton = Button(self.frame2, text='...', command=self.loadref)
        self.reffilebutton.grid(column=2, row=1)
        self.blastfilelabel = Label(self.frame2, text='Comparison file:', anchor=E)
        self.blastfilelabel.grid(column=0, row=2, sticky=E)
        self.blastfileentry = Entry(self.frame2, textvariable=self.blastfile)
        self.blastfileentry.grid(column=1, row=2)
        self.blastfilebutton = Button(self.frame2, text='...', command=self.loadblast)
        self.blastfilebutton.grid(column=2, row=2)
        self.minlengthblastlabel = Label(self.frame2, text='Minimum length (bp):', anchor=E)
        self.minlengthblastlabel.grid(column=0, row=3, sticky=E)
        self.minlengthblastentry = Entry(self.frame2, textvariable=self.minlengthblast)
        self.minlengthblastentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.minlengthratiolabel = Label(self.frame2, text='Minimum length (ratio):', anchor=E)
        self.minlengthratiolabel.grid(column=0, row=4, sticky=E)
        self.minlengthratioentry = Entry(self.frame2, textvariable=self.minlengthratio)
        self.minlengthratioentry.grid(column=1, row=4, columnspan=2, sticky=EW)
        self.minidentlabel = Label(self.frame2, text='Minimum identity:', anchor=E)
        self.minidentlabel.grid(column=0, row=5, sticky=E)
        self.minidententry = Entry(self.frame2, textvariable=self.minident)
        self.minidententry.grid(column=1, row=5, columnspan=2, sticky=EW)
        self.minbitscorelabel = Label(self.frame2, text='Minimum bitscore:', anchor=E)
        self.minbitscorelabel.grid(column=0, row=6, sticky=E)
        self.minbitscoreentry = Entry(self.frame2, textvariable=self.minbitscore)
        self.minbitscoreentry.grid(column=1, row=6, columnspan=2, sticky=EW)
        self.maxevaluelabel = Label(self.frame2, text='Maximum evalue:', anchor=E)
        self.maxevaluelabel.grid(column=0, row=7, sticky=E)
        self.maxevalueentry = Entry(self.frame2, textvariable=self.maxevalue)
        self.maxevalueentry.grid(column=1, row=7, columnspan=2, sticky=EW)
        self.closeblastbut = Button(self.frame2, text='Ok', command=self.okblast)
        self.closeblastbut.grid(column=2, row=8, sticky=E, pady=5)
        self.frame2.grid(padx=10, pady=10)

    def okblast(self):
        if self.reffile.get() == '':
            tkMessageBox.showerror('File missing', 'Please choose a reference file.', parent=self.blast_options)
            return
        if os.path.exists(self.reffile.get()):
            temp = open(self.reffile.get())
            if temp.readline()[0] != '>':
                tkMessageBox.showerror('File not fasta', 'Please choose a valid reference file.', parent=self.blast_options)
                return
        else:
            tkMessageBox.showerror('File doesn\'t exist', 'Please choose a valid reference file.', parent=self.blast_options)
            return
        if self.contigDict == {}:
            tkMessageBox.showerror('Assembly not found', 'Please load a FASTA or CAG file using load assembly before creating comparison.', parent=self.blast_options)
            self.blast_options.destroy()
            return
        if self.blastfile.get() == '':
            self.blastit = tkMessageBox.askquestion('No blast files.', 'Create BLAST output?')
            if self.blastit == 'yes':
                if (not which('blastn.exe') or not which('makeblastdb.exe')) and (not which('blastn') or not which('makeblastdb')):
                    tkMessageBox.showerror('BLAST not found', 'Please install NCBI-BLAST to your path, or perform the comparison yourself.')
                    return
            else:
                tkMessageBox.showerror('No Comparison', 'Please choose valid comparison file.')
                return
        else:
            if not os.path.exists(self.blastfile.get()):
                tkMessageBox.showerror('No Comparison', 'Please choose valid comparison file.')
                return
            else:
                self.blastit = 'no'
        self.blast_options.destroy()
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        self.queue = Queue.Queue()
        self.thethread = threading.Thread(target=self.create_comp_thread)
        self.thethread.start()
        self.abort = False
        self.update_console('Creating comparison')
        self.update_create_comp()

    def update_create_comp(self):
        self.dot_console()
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text != 'Comparison created.':
                    root.after(1000, self.update_create_comp)
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            if self.abort:
                self.update_console('Thread stopped.')
            else:
                self.update_console('Creating comparison failed, please check console output.')
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_create_comp)
            return

    def create_comp_thread(self):
        ref = open(self.reffile.get())
        first = True
        self.refDict = {}
        self.reforder = []
        refout = open(self.workingDir.get() + '/reference.fa', 'w')
        for line in ref:
            if line.startswith('>'):
                if first:
                    first = False
                else:
                    self.refDict[name] = seq
                    self.reforder.append(name)
                    refout.write('>' + name + '\n' + seq + '\n')
                name = line.split()[0][1:]
                seq = ''
            else:
                seq += line.rstrip()
        self.refDict[name] = seq
        refout.write('>' + name + '\n' + seq + '\n')
        refout.close()
        self.reforder.append(name)
        if self.blastit == 'yes':
            self.queue.put('Running BLAST.')
            subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/tempdb -in ' + self.workingDir.get() + '/reference.fa', stdout=subprocess.PIPE, shell=True).wait()
            subprocess.Popen('blastn -task blastn -db ' + self.workingDir.get() + '/tempdb -outfmt 6 -query ' + self.workingDir.get() + '/contigs.fa -out ' + self.workingDir.get() + '/query_tempdb.out', shell=True).wait()
            self.blastfile.set(self.workingDir.get() + '/query_tempdb.out')
            self.queue.put('BLAST file created.')
        if self.blastfile.get() != '':
            blastfile = open(self.blastfile.get())
            self.hitlist = []
            for line in blastfile:
                query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
                qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
                ident, eval, bitscore = map(float, [ident, eval, bitscore])
                if ident >= self.minident.get() and length >= self.minlengthblast.get() and eval <= self.maxevalue.get() and bitscore >= self.minbitscore.get()\
                    and (length >= self.minlengthratio.get() * self.contigDict[query].length or length >= self.minlengthratio.get() * len(self.refDict[subject])):
                    self.hitlist.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
        self.queue.put('Comparison created.')

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
        self.reffilelabel = Label(self.frame2, text='Comparison file:', anchor=E)
        self.reffilelabel.grid(column=0, row=1, sticky=E)
        self.reffileentry = Entry(self.frame2, textvariable=self.selfcomparefile)
        self.reffileentry.grid(column=1, row=1)
        self.reffilebutton = Button(self.frame2, text='...', command=self.loadselfcomp)
        self.reffilebutton.grid(column=2, row=1)
        self.intralabel = Label(self.frame2, text='Show intra contig hits:', anchor=E)
        self.intralabel.grid(column=0, row=2, sticky=E)
        self.intraentry = Checkbutton(self.frame2, variable=self.intra)
        self.intraentry.grid(column=1, row=2, columnspan=2, sticky=EW)
        self.minlengthblastlabel = Label(self.frame2, text='Minimum length (bp):', anchor=E)
        self.minlengthblastlabel.grid(column=0, row=3, sticky=E)
        self.minlengthblastentry = Entry(self.frame2, textvariable=self.minlengthblast)
        self.minlengthblastentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.minlengthratiolabel = Label(self.frame2, text='Minimum length (ratio):', anchor=E)
        self.minlengthratiolabel.grid(column=0, row=4, sticky=E)
        self.minlengthratioentry = Entry(self.frame2, textvariable=self.minlengthratio)
        self.minlengthratioentry.grid(column=1, row=4, columnspan=2, sticky=EW)
        self.minidentlabel = Label(self.frame2, text='Minimum identity:', anchor=E)
        self.minidentlabel.grid(column=0, row=5, sticky=E)
        self.minidententry = Entry(self.frame2, textvariable=self.minident)
        self.minidententry.grid(column=1, row=5, columnspan=2, sticky=EW)
        self.minbitscorelabel = Label(self.frame2, text='Minimum bitscore:', anchor=E)
        self.minbitscorelabel.grid(column=0, row=6, sticky=E)
        self.minbitscoreentry = Entry(self.frame2, textvariable=self.minbitscore)
        self.minbitscoreentry.grid(column=1, row=6, columnspan=2, sticky=EW)
        self.maxevaluelabel = Label(self.frame2, text='Maximum evalue:', anchor=E)
        self.maxevaluelabel.grid(column=0, row=7, sticky=E)
        self.maxevalueentry = Entry(self.frame2, textvariable=self.maxevalue)
        self.maxevalueentry.grid(column=1, row=7, columnspan=2, sticky=EW)
        self.intralabel = Label(self.frame2, text='Only show edge hits:', anchor=E)
        self.intralabel.grid(column=0, row=8, sticky=E)
        self.intraentry = Checkbutton(self.frame2, variable=self.onlyedge)
        self.intraentry.grid(column=1, row=8, columnspan=2, sticky=EW)
        self.maxevaluelabel = Label(self.frame2, text='Max. distance from edge:', anchor=E)
        self.maxevaluelabel.grid(column=0, row=9, sticky=E)
        self.maxevalueentry = Entry(self.frame2, textvariable=self.maxedge)
        self.maxevalueentry.grid(column=1, row=9, columnspan=2, sticky=EW)
        self.closeblastbut = Button(self.frame2, text='Ok', command=self.ok_self_compare)
        self.closeblastbut.grid(column=2, row=10, sticky=E, pady=5)
        self.frame2.grid(padx=10, pady=10)

    def ok_self_compare(self):
        try:
            if self.thethread.is_alive():
                tkMessageBox.showerror('Already running process',
                                       'Please wait until current tasks have finished before running another process.')
                return
        except:
            pass
        self.blast_options.destroy()
        if self.selfcomparefile.get() == '':
            if (not which('blastn.exe') or not which('makeblastdb.exe')) and (not which('blastn') or not which('makeblastdb')):
                tkMessageBox.showerror('BLAST not found', 'Please install NCBI-BLAST or include own comparison file.')
                return
        if self.visible == set():
            tkMessageBox.showerror('No contigs visible.', 'Please add contigs to compare to main canvas.')
        self.canvas.delete('selfhit')
        self.queue = Queue.Queue()
        self.thethread = threading.Thread(target=self.self_hits_thread)
        self.thethread.start()
        self.abort = False
        self.update_console('Finding self hits.')
        self.update_self_hits()

    def update_self_hits(self):
        self.dot_console()
        while self.queue.qsize():
            try:
                text = self.queue.get(0)
                self.update_console(text)
                if text[-11:] != 'hits found.':
                    root.after(1000, self.update_self_hits)
                else:
                    self.draw_self_hits()
                return
            except Queue.Empty:
                pass
        if not self.thethread.is_alive():
            if self.abort:
                self.update_console('Thread stopped.')
            else:
                self.update_console('Self comparison failed, please check console output.')
            return
        elif not self.queue.qsize():
            root.after(1000, self.update_self_hits)
            return

    def self_hits_thread(self):
        if self.selfcomparefile.get() == '':
            stderr = open(self.workingDir.get() + '/bwaerr.txt', 'w')
            subprocess.Popen('makeblastdb -dbtype nucl -out ' + self.workingDir.get() + '/contigdb -in ' +
                             self.workingDir.get() + '/contigs.fa', shell=True, stdout=stderr).wait()
            subprocess.Popen('blastn -db ' + self.workingDir.get() + '/contigdb -outfmt 6 -num_threads 8 -query ' +
                             self.workingDir.get() + '/contigs.fa -out ' + self.workingDir.get() + '/contigscontigs.out', shell=True).wait()
            stderr.close()
            blastfile = open(self.workingDir.get() + '/contigscontigs.out')
        else:
            blastfile = open(self.selfcomparefile.get())
        self.queue.put('Comparison file created.')
        self.selfhit = []
        self.queue.put('Getting hits.')
        for line in blastfile:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length, mm = map(int, [qstart, qstop, rstart, rstop, length, mm])
            ident, eval, bitscore = map(float, [ident, eval, bitscore])
            if ident >= self.minident.get() and length >= self.minlengthblast.get() and eval <= self.maxevalue.get() and bitscore >= self.minbitscore.get()\
              and (length >= self.minlengthratio.get() * self.contigDict[query].length or length >= self.minlengthratio.get() * self.contigDict[subject].length):
                if query != subject:
                    if query < subject:
                        if self.onlyedge.get() == 1:
                            if (qstart <= self.maxedge.get() or qstop >= self.contigDict[query].length - self.maxedge.get()) and \
                              (rstart <= self.maxedge.get() or rstop <= self.maxedge.get() or
                              rstart >= self.contigDict[subject].length - self.maxedge.get() or
                              rstop >= self.contigDict[subject].length - self.maxedge.get()):
                                self.selfhit.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
                                if qstart <= self.maxedge.get() and rstart < rstop and rstop >= self.contigDict[subject].length -self.maxedge.get():
                                    self.contigDict[query].fr.append((subject, False, qstop + self.contigDict[subject].length - rstop))
                                    self.contigDict[subject].to.append((query, True, qstop + self.contigDict[subject].length - rstop))
                                if qstart <= self.maxedge.get() and rstart > rstop and rstop <= self.maxedge.get():
                                    self.contigDict[query].fr.append((subject, True, qstop + rstop -1))
                                    self.contigDict[subject].fr.append((query, True, qstop + rstop -1))
                                if qstop >= self.contigDict[query].length - self.maxedge.get() and rstart < rstop and rstart <= self.maxedge.get():
                                    self.contigDict[query].to.append((subject, True, rstop + self.contigDict[query].length - qstop))
                                    self.contigDict[subject].fr.append((query, False, rstop + self.contigDict[query].length - qstop))
                                if qstop >= self.contigDict[query].length - self.maxedge.get() and rstart > rstop and \
                                  rstart >= self.contigDict[subject].length -self.maxedge.get():
                                    self.contigDict[query].to.append((subject, False, self.contigDict[subject].length - rstop + self.contigDict[query].length - qstop))
                                    self.contigDict[subject].to.append((query, False, self.contigDict[subject].length - rstop + self.contigDict[query].length - qstop))
                        else:
                            self.selfhit.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
                elif self.intra.get() == 1:
                    if qstart < min([rstart, rstop]):
                        if self.onlyedge.get() == 1:
                            if (qstart <= self.maxedge.get() or qstop >= self.contigDict[query].length - self.maxedge.get()) and \
                              (rstart <= self.maxedge.get() or rstop <= self.maxedge.get() or
                              rstart >= self.contigDict[subject].length - self.maxedge.get() or
                              rstop >= self.contigDict[subject].length - self.maxedge.get()):
                                self.selfhit.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))
                        else:
                            self.selfhit.append((query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore))

        blastfile.close()
        self.queue.put(str(len(self.selfhit)) + ' hits found.')

    def draw_self_hits(self):
        # need to add in hits for duplicates
        self.canvas.delete('selfhit')
        hitnum = 0
        self.selfhit.sort(key=lambda x: x[3], reverse=True)
        for i in self.selfhit:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = i
            hitnum += 1
            if query in self.visible and subject in self.visible:
                ratio = (ident - self.minident.get()) / (100 - self.minident.get())
                if rstart < rstop and self.contigDict[query].orient[0] == self.contigDict[subject].orient[0]:
                    colour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
                elif rstart < rstop:
                    colour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
                elif rstop < rstart and self.contigDict[query].orient[0] == self.contigDict[subject].orient[0]:
                    colour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
                else:
                    colour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
                qcoords = self.canvas.coords('c' + query)
                rcoords = self.canvas.coords('c' + subject)
                if self.contigDict[query].orient[0]:
                    startx1 = qcoords[0] + qstart * 1.0 / self.contigDict[query].length * abs(qcoords[2] - qcoords[0])
                    startx2 = qcoords[0] + qstop * 1.0 / self.contigDict[query].length * abs(qcoords[2] - qcoords[0])
                else:
                    startx1 = qcoords[2] - qstart * 1.0 / self.contigDict[query].length * abs(qcoords[2] - qcoords[0])
                    startx2 = qcoords[2] - qstop * 1.0 / self.contigDict[query].length * abs(qcoords[2] - qcoords[0])
                if self.contigDict[subject].orient[0]:
                    endx1 = rcoords[0] + rstart * 1.0 / self.contigDict[subject].length * abs(rcoords[2] - rcoords[0])
                    endx2 = rcoords[0] + rstop * 1.0 / self.contigDict[subject].length * abs(rcoords[2] - rcoords[0])
                else:
                    endx1 = rcoords[2] - rstart * 1.0 / self.contigDict[subject].length * abs(rcoords[2] - rcoords[0])
                    endx2 = rcoords[2] - rstop * 1.0 / self.contigDict[subject].length * abs(rcoords[2] - rcoords[0])
                starty = qcoords[1] + self.contigheight / 2
                endy = rcoords[1] + self.contigheight / 2
                self.canvas.create_polygon(startx1, starty, startx1, starty, (startx1 + endx1) / 2, abs(startx1 - endx1) /4 + (starty + endy) / 2, endx1, endy,
                                           endx1, endy, endx2, endy, endx2, endy, (startx2 + endx2) / 2, abs(startx2 - endx2) /4 + (starty + endy) / 2,
                                           startx2, starty, startx2, starty, smooth=1, fill=colour, outline='black', tags=('c' + query + 'ss', 'c' + subject + 'se', 'selfhit', str(hitnum)))
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
        self.viewlabel = Label(self.frame3, text='Contigs to view:', anchor=E)
        self.viewlabel.grid(column=0, row=0, sticky=E)
        self.viewentry = OptionMenu(self.frame3, self.view, 'All', 'BLAST', 'List', 'Filter')#, 'Coverage', 'BLAST and Coverage', 'Blast or Coverage')
        self.viewentry.grid(column=1, row=0, columnspan=2, sticky=EW)
        self.viewlistlabel = Label(self.frame3, text='List of contigs:', anchor=E)
        self.viewlistlabel.grid(column=0, row=1, sticky=E)
        self.viewlistentry = Entry(self.frame3, textvariable=self.viewlist)
        self.viewlistentry.grid(column=1, row=1)
        self.viewlistbutton = Button(self.frame3, text='...', command=self.loadview)
        self.viewlistbutton.grid(column=2, row=1)
        self.viewreflabel = Label(self.frame3, text='View reference and BLAST hits:', anchor=E)
        self.viewreflabel.grid(column=0, row=2, sticky=E)
        self.viewrefentry = Checkbutton(self.frame3, variable=self.viewref)
        self.viewrefentry.grid(column=1, row=2, columnspan=2)
        self.shortenlabel = Label(self.frame3, text='Shorten long contigs:', anchor=E)
        self.shortenlabel.grid(column=0, row=3, sticky=E)
        self.shortenentry = OptionMenu(self.frame3, self.shorten, 'No', 'Log', 'Max', 'Min')
        self.shortenentry.grid(column=1, row=3, columnspan=2, sticky=EW)
        self.scaledownlabel = Label(self.frame3, text='Scale: 1 pixel =', anchor=E)
        self.scaledownlabel.grid(column=0, row=4, sticky=E)
        self.scaledownenentry = Entry(self.frame3, textvariable=self.scaledown)
        self.scaledownenentry.grid(column=1, row=4)
        self.zebplabel = Label(self.frame3, text='bp')
        self.zebplabel.grid(column=2, row=4)
        self.okview = Button(self.frame3, text='Ok', command=self.ok_view)
        self.okview.grid(column=2, row=5, sticky=E)
        self.frame3.grid(padx=10, pady=10)

    def ok_view(self):
        self.canvas.delete(ALL)
        self.newscaledown = self.scaledown.get()
        self.visible = set()
        self.clear_lists()
        self.contigheight = 25
        self.currxscroll = self.originalxscroll
        self.curryscroll = self.originalyscroll
        self.fontsize = 12
        self.customFont.configure(size=12)
        self.canvas.config(scrollregion=(0, 0, self.currxscroll, self.curryscroll))
        self.rightmost = None
        self.leftmost = None
        if self.shorten.get() == 'No':
            for i in self.contigDict:
                self.contigDict[i].xlength = self.contigDict[i].length / self.scaledown.get()
                self.contigDict[i].scaledown = self.scaledown.get()
        if self.shorten.get() == 'Log':
            for i in self.contigDict:
                self.contigDict[i].xlength = max([math.log10(self.contigDict[i].length / self.scaledown.get()), 1]) * self.scaledown.get()
                self.contigDict[i].scaledown = self.contigDict[i].length * 1.0 / (max([math.log10(self.contigDict[i].length / self.scaledown.get()), 1]) * self.scaledown.get())
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
            if self.viewref.get() and not self.hitlist is None and not self.reforder is None:
                self.orderContigsAllBlast()
            else:
                self.orderContigsVisible()
            self.drawContigs()
            self.drawEdges()
        elif self.view.get() == 'List':
            self.orderContigsList()
            self.drawContigs()
            self.drawEdges()
        elif self.view.get() == 'Filter':
            if self.hitlist is None:
                tkMessageBox.showerror('Comparison not found.', 'Please perform comparison before choosing this method.')
                return
            self.filterBlast()
            self.orderContigsVisible()
            self.drawContigs()
            self.drawEdges()
        if self.viewref.get() and not self.hitlist is None and not self.reforder is None and not self.view.get() == 'Filter':
            self.drawRefHits()
        self.canvas.move(ALL, -self.leftmost + 50, 0)
        self.currxscroll = self.rightmost - self.leftmost + 200
        self.curryscroll = self.currxscroll
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

    def filterBlast(self):
        self.visible = set(self.contigDict)
        for i in self.hitlist:
            if i[0] in self.visible:
                self.visible.remove(i[0])

    def drawRefHits(self):
        self.hitlist.sort(key=lambda x: x[3], reverse=True)
        hitcount = 0
        for i in self.hitlist:
            hitcount += 1
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = i
            if query in self.visible:
                ratio = (ident - self.minident.get()) / (100 - self.minident.get())
                if rstart < rstop:
                    colour = '#%02x%02x%02x' % (168 - int(ratio * 168), 223 - int(ratio * 70), 244 - int(ratio * 40))
                else:
                    colour = '#%02x%02x%02x' % (255, 216 - int(ratio * 80), 131 - int(ratio * 131))
                self.canvas.create_polygon(self.refpos[subject] + rstart / self.scaledown.get(), self.refline + self.contigheight,
                                           self.refpos[subject] + rstop / self.scaledown.get(), self.refline + self.contigheight,
                                            self.contigDict[query].xpos + qstop/self.contigDict[query].scaledown, self.contigline,
                                             self.contigDict[query].xpos + qstart/self.contigDict[query].scaledown, self.contigline,
                                             fill=colour, outline="black", tags=('c' + query + 'b', 'r' + subject + 'b', 'blast', str(hitcount)))
        for i in self.reforder:
            if self.leftmost is None or self.refpos[i] < self.leftmost:
                self.leftmost = self.refpos[i]
            if self.rightmost is None or self.refpos[i] + len(self.refDict[i]) / self.scaledown.get() > self.rightmost:
                self.rightmost = self.refpos[i] + len(self.refDict[i]) / self.scaledown.get()
            self.canvas.create_rectangle(self.refpos[i], self.refline, \
                                         self.refpos[i] + len(self.refDict[i]) / self.scaledown.get(), self.refline + self.contigheight, \
                                         fill="#669900", tags=('r' + i, "ref", "map"))
            if len(self.refDict[i]) > self.scaledown.get() * 100 * 2:
                markers = self.scaledown.get() * 500
                for j in range(markers, len(self.refDict[i]), markers):
                    self.canvas.create_line(self.refpos[i] + j / self.scaledown.get(), self.refline + 2, self.refpos[i] + j / self.scaledown.get(), self.refline + self.contigheight - 2,
                                            fill='white', tags=('r' + i, 'ref', 'map'))
                    thetext = str(j)
                    self.canvas.create_text(self.refpos[i] + j / self.scaledown.get() + 2, self.refline + self.contigheight/2,
                     fill='white', font=self.customFont, anchor=W, text=thetext, tags=('r' + i, 'ref', 'map'))
            thetext = i
            text = self.canvas.create_text(self.refpos[i] + 2, self.refline + self.contigheight/2,\
                                            fill='white', font=self.customFont, anchor=W, text=thetext, tags=('r' + i, 'ref', 'map'))
            if self.canvas.bbox(text)[2] >= self.refpos[i] + len(self.refDict[i]) / self.scaledown.get() - 2:
                self.canvas.delete(text)

    def drawContigs(self):
        for i in self.visible:
            if self.leftmost is None or self.contigDict[i].xpos < self.leftmost:
                self.leftmost = self.contigDict[i].xpos
            if self.rightmost is None or self.contigDict[i].xpos + self.contigDict[i].xlength > self.rightmost:
                self.rightmost = self.contigDict[i].xpos + self.contigDict[i].xlength
            self.canvas.create_rectangle(self.contigDict[i].xpos, self.contigDict[i].ypos,
                                         self.contigDict[i].xpos + self.contigDict[i].xlength, self.contigDict[i].ypos + self.contigheight,\
                                          fill="#9933CC", tags=('c' + i, 'contig', 'map'))
            if self.contigDict[i].orient[0]:
                dir = '+'
            else:
                dir = '-'
            thetext = self.contigDict[i].shortname + ' ' + dir + ' ' + self.contigDict[i].strlen
            text = self.canvas.create_text(self.contigDict[i].xpos + 2, self.contigDict[i].ypos + self.contigheight/2,\
                                            fill='white', font=self.customFont, anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
            if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                self.canvas.delete(text)
                thetext = self.contigDict[i].shortname + ' ' + dir
                text = self.canvas.create_text(self.contigDict[i].xpos + 2, self.contigDict[i].ypos + self.contigheight/2,\
                                                fill='white', font=self.customFont, anchor=W, text=thetext, tags=('c' + i, 'map', 'text', 'c' + i + 't'))
                if self.canvas.bbox(text)[2] >= self.contigDict[i].xpos + self.contigDict[i].xlength -2:
                    self.canvas.delete(text)
                    thetext = self.contigDict[i].shortname
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
            currx += self.contigDict[name].xlength + 10
        listfile.close()

    def orderContigsVisible(self):
        todo = set(self.visible)
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
                if length == besthit[query][1]:
                    if self.refpos[subject] + min((rstart, rstop)) / self.scaledown.get() < besthit[query][0]:
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

    def orderContigsAllBlast(self):
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
        todo = set(self.contigDict)
        for i in qpos:
            self.contigDict[i[1]].xpos = currxpos
            self.contigDict[i[1]].ypos = self.contigline
            currxpos += self.contigDict[i[1]].xlength + 10
            todo.remove(i[1])
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
        qpos = map(lambda y: (y[0] + currxpos, y[1]), qpos)
        for i in qpos:
            self.contigDict[i[1]].xpos = i[0]
            self.contigDict[i[1]].ypos = self.contigline

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
        self.okedges.grid(row=5, column=2, sticky=E)
        self.frame4.grid(padx=10, pady=10)

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
        self.outfilelabel = Label(self.frame5, text='Write to file:', anchor=E)
        self.outfilelabel.grid(column=0, row=1, sticky=E)
        self.outfileentry = Entry(self.frame5, textvariable=self.outfile)
        self.outfileentry.grid(column=1, row=1)
        self.outfilebutton = Button(self.frame5, text='...', command=self.loadoutfile)
        self.outfilebutton.grid(column=2, row=1)
        self.bufferlabel = Label(self.frame5, text='Buffer sequence:', anchor=E)
        self.bufferlabel.grid(column=0, row=2, sticky=E)
        self.bufferentry = Entry(self.frame5, textvariable=self.buffer)
        self.bufferentry.grid(column=1, row=2, columnspan=2, sticky=EW)
        self.okedges = Button(self.frame5, text='Ok', command=self.ok_fasta)
        self.okedges.grid(row=3, column=2, sticky=E)
        self.frame5.grid(padx=10, pady=10)

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


parser = argparse.ArgumentParser(prog='Contiguity', formatter_class=argparse.RawDescriptionHelpFormatter, description='''
Contiguity.py: A pairwise comparison and contig adjacency graph exploration tool.

Version: 1.0.4
License: GPLv3

USAGE: Contiguity.py -cl -c <contig_file.fa> -fq <read_file.fq> -o <output_folder>

contig file: FASTA file of contigs or scaffolds
read file: Interleaved fastq file - read1_left, read1_right, read2_left etc... orientated as such --> <--
output folder: folder to put output files in, can and will overwrite files in this folder, will create folder if folder doesn't exist

REQUIREMENTS: With default settings Contigutiy requires at least (Potentially more) 6gb of free memory (RAM). This is
because Contiguity uses Khmer to create a De Bruijn graph.

If you are running into memory issues (e.g. the process is being killed by the Operating system (killed 9))
you may want to reduce the hash table number or hash table size (using the flags -ht_n and -ht_s), or free up more memory.
For large datasets more memory may need to be used.
please read http://khmer.readthedocs.org/en/v1.1/choosing-table-sizes.html for more information about hash tables

The only other option to keep in mind is -rl if the read length is not 101bp, you may want to increase this value
for longer read lengths, or decrease for shorter read lengths (75% of maximum read length seems to work well).

''', epilog="Thanks for using Contiguity")
parser.add_argument('-co', '--contig_file', action='store', help='fasta file of assembled contigs or scaffolds')
parser.add_argument('-rf', '--read_file', action='store', help='read file')
parser.add_argument('-o', '--output_folder', action='store', help='output folder')
parser.add_argument('-k', '--kmer_size', action='store', type=int, default=31, help='k-mer size for finding adjacent contigs [31]')
parser.add_argument('-max_d', '--max_distance', action='store', type=int, default=200, help='maximum distance apart in the de bruijn graph for contigs to count as adjacent [300]')
parser.add_argument('-kmer_a', '--kmer_average', action='store', type=int, default=-1, help='All k-mers above half this value will be traversed [auto]')
parser.add_argument('-kmer_c', '--kmer_cutoff', action='store', type=int, default=-1, help='cutoff for k-mer values [auto]')
parser.add_argument('-ov', '--overlap', action='store', type=int, default=None, help='minimum overlap to create edge [kmer_size-1]')
parser.add_argument('-rl', '--min_read_length', action='store', type=int, default=75, help='Minimum read length [75]')
parser.add_argument('-max_mm', '--max_mismatch', action='store', type=int, default=2, help='maximum number of mismatches to count overlap [2]')
parser.add_argument('-lo', '--long_overlap_ident', action='store', type=int, default=85, help='minimum percent identity to create an edge where there is a long overlap [85]')
parser.add_argument('-mp', '--minimum_pairs_edge', action='store', type=int, default=2, help='Minimum pairs to create edge [2]')
parser.add_argument('-is', '--max_insert_size', action='store', type=int, default=600, help='Upper bound on insert size [600]')
parser.add_argument('-cl', '--command_line', action='store_true', default=False, help='Run contiguity in command line mode')
parser.add_argument('-no', '--no_overlap_edges', action='store_true', default=False, help='Don\'t get overlap edges')
parser.add_argument('-nd', '--no_db_edges', action='store_true', default=False, help='Don\'t get De Bruijn edges')
parser.add_argument('-np', '--no_paired_edges', action='store_true', default=False, help='Don\'t get paired-end edges')
parser.add_argument('-km', '--khmer', action='store_false', default=True, help='Don\'t use khmer for De Bruijn graph contruction (not recommended)')
parser.add_argument('-nt', '--num_threads', action='store', type=int, default=1, help='Number of threads to use for mapping reads with bowtie [1]')
parser.add_argument('-ht_s', '--ht_size', action='store', default='1e9', help='Hash table size.')
parser.add_argument('-ht_n', '--ht_number', action='store', type=int, default=4, help='Hash table number.')



args = parser.parse_args()

if args.command_line:
    if platform.system() == 'Windows':
        args.khmer = False
    if args.khmer:
        import khmer
    if args.contig_file is None or args.read_file is None or args.output_folder is None:
        sys.stdout.write("Command line CAG building requires a contig file [-co] a read file [-r] and an output folder [-o].\n")
        sys.exit()
    if args.khmer and args.kmer_size > 32:
        args.kmer_size = 32
        sys.stdout.write("Khmer's maximum kmer size is 32, setting kmer size to 32")
    if not args.overlap is None and args.overlap > args.kmer_size:
        sys.stdout.write("WARNING: It is strongly recommended you use a overlap size 1 less than your kmer size.")



    theapp = App(None)
else:
    root = Tk()
    root.title('Contiguity')
    root.option_add('*Font', 'Heveltica 10')
    root.option_add("*Background", "#E0E0FF")
    root.option_add("*Foreground", "#2B3856")
    root.option_add("*Listbox.Background", '#FFFFFF')
    root.option_add("*Scrollbar.Background", "#C0C0FF")
    root.option_add("*Entry.Background", "#FFFFFF")
    app = App(root)
    root.mainloop()