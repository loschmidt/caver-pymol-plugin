#!/usr/bin/env python

# CAVER Copyright Notice
# ============================
#

from __future__ import division
from __future__ import generators

import re
import os,math

try:
    import tkinter as tk
    from tkinter import *
except:
    import Tkinter as tk
    from Tkinter import *

try:
    from tkinter import filedialog
except:
    import tkFileDialog as filedialog

import Pmw
import distutils.spawn # used for find_executable
from pymol import cmd,selector
import sys
from pymol.cmd import _feedback,fb_module,fb_mask,is_list,_cmd
from pymol.cgo import *
from chempy.models import Indexed
from chempy import Bond, Atom
import threading
#import subprocess

import shutil
from pymol import stored
import time
#
# Global config variables
#

#win/linux ========================
# 1 for windows, 0 for linux

def import_file(full_path_to_module):
    import os
    module_dir, module_file = os.path.split(full_path_to_module)
    module_name, module_ext = os.path.splitext(module_file)

    if module_name in sys.modules:
        del(sys.modules[module_name])

    save_cwd = os.getcwd()
    os.chdir(module_dir)
    module_obj = __import__(module_name)
    module_obj.__file__ = full_path_to_module
    globals()[module_name] = module_obj
    os.chdir(save_cwd)

VERS_M = "3"
VERS_0 = "0"
VERS_1 = "3"
#JOPTS = "-Xmx@m" # @ is going to be replaced by user-specified value
#JHEAP = 1100

VERSION = "%s.%s.%s" % (VERS_M,VERS_0,VERS_1)
VERSION_ = "%s_%s_%s" % (VERS_M,VERS_0,VERS_1)

CAVER3_LOCATION = os.path.dirname(__file__)

OUTPUT_LOCATION = os.path.abspath(".")

LABEL_TEXT = "Caver directory:"

#
# Cheap hack for testing purposes
#
try:
    import pymol
    REAL_PYMOL = True
except ImportError:
    REAL_PYMOL = False
    class pymol:
        class cmd:
            def load(self,name,sel=''):
                pass
            def get_names(self):
                return ['mol1','mol2','map1','map2']
            def get_type(self,thing):
                if thing.startswith('mol'):
                    return 'object:molecule'
                else:
                    return 'object:map'
                f.close()
        cmd = cmd()
    pymol = pymol()


# pridani do menu
def __init__(self):
    lbb = "Caver %s" % (VERSION,)
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Launch Caver '  + VERSION,
                             label=lbb,
                             command = lambda s=self: AnBeKoM(s))




defaults = {
    "startingpoint": ('0','0','0'),
    "compute_command": 'Compute tunnels',
    "warn_command": 'Show warnings',
    "exit_command": 'Exit',
    "default_shell_depth": '2',
    "default_shell_radius": '3.0',
    "default_tunnels_probe": '0.7',
    "default_java_heap": '6000',
    "default_clustering_threshold": '1.5',
    "surroundings" : 'sele',
    "startingacids":('117','283','54'),
    "default_block": '10.0'
    }

url = "http://www.caver.cz/index.php?sid=123"

class MyThread (threading.Thread):
    def run (self):
        os.system("start " + url)


class DataStruct:
    def __init__(self):
        self.keys = []
        self.values = []
    def remove(self, key):
        idx = self.indexOf(key)
        if idx != -1: self.keys[idx] = "REMOVED"
        #self.keys.pop(idx)
        #self.values.pop(idx)
    def indexOf(self, key):
        #for idx in range (0, len(self.keys)):
        #  if self.keys[idx] == key:
        #    return idx
        #return -1
        if key in self.keys:
            return self.keys.index(key)
        else:
            return -1
    def add(self, key, value, isComment):
        idx = self.indexOf(key)
        if idx == -1 or isComment == 1:
            self.keys.append(key)
            self.values.append(value)
        else:
            self.values[idx] = self.values[idx] + " " + value
    def replace(self, key, value, isComment):
        idx = self.indexOf(key)
        if idx == -1 or isComment:
            #print("replacing " + key + " idx " + str(idx))
            self.keys.append(key)
            self.values.append(value)
            #print("size " + str(len(self.keys)))
        else:
            #print("Exists " + key + " idx " + str(idx))
            self.values[idx] = value
    def get(self, key):
        idx = self.indexOf(key)
        return self.values[idx]
    def getKeys(self):
        return self.keys
    def getValues(self):
        return self.values
    def clear(self):
        self.keys = []
        self.values = []

class PyJava:

    def status(self, r):
        if 0 == r:
            print("OK")
        else:
            print("FAIL")

    def __init__(self, maxXmx, caverfolder, caverjar, outdirInputs, cfgnew, out_dir):
        self.insufficient_memory = False
        self.jar = caverjar
        print("")
        print("*** Testing if Java is installed ***")
        r = self.java_present()
        self.status(r)

        self.java_missing = bool(r)
        if r:
            return

        print("")
        print("*** Optimizing memory allocation for Java ***")
        self.optimize_memory(maxXmx)
        self.cmd = [
            "java",
            "-Xmx%dm" % self.xmx,
            "-cp", os.path.join(caverfolder, "lib"),
            "-jar", caverjar,
            "-home", caverfolder,
            "-pdb", outdirInputs,
            "-conf", cfgnew,
            "-out", out_dir,
        ]
        print("*** Caver will be called using command ***")
        print(self.cmd)
        print("******************************************")

    def java_present(self):
        cmd = ["java", "-version"]
        r = self.execute(cmd, False)
        return r

    def run_caver(self):
        self.execute(self.cmd, False)

    def optimize_memory(self, s_max_xmx):
        max_xmx = int(s_max_xmx)
        values = [500, 800, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1400, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 14000, 16000, 20000, 32000, 48000, 64000]
        values.append(max_xmx)
        values.sort()
        #sorted(values)
        self.xmx = values[0]
        for xmx in values:
            if int(xmx) <= max_xmx:
                cmd = ["java", "-Xmx%dm" % xmx, "-jar", self.jar, "do_nothing"]
                code = self.execute(cmd, True)
                if 0 == code:
                    self.xmx = xmx
                    print("Xmx: " + str(self.xmx))
        print("*** Memory for Java: " + str(self.xmx) + " MB ***")
        print

    def execute_old(self, cmd):
        p = os.popen(cmd)
        for line in p.readlines():
            print(line.rstrip())
        p.close()
        return 1

    def execute(self, args, silent):
        if True:
            import subprocess
            try:
                p = subprocess.check_output(args, stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
                if not silent:
                    print(p)
            except subprocess.CalledProcessError as e:
                if not silent:
                    print(e)
                    print(e.cmd)
                    print(e.output)
                self.analyze(e.output.decode('UTF-8'))
                return e.returncode
            except OSError as e:
                error_dialog = Pmw.MessageDialog(title='Error',
                    message_text="Can't execute " + str(args) + "\n\n" + str(e))
                return -1
            except Exception as e:
                error_dialog = Pmw.MessageDialog(title='Error',
                    message_text="Unknown error: " + str(e))
                return -2
            return 0

    def analyze(self, output):
        if 'OutOfMemory' in output:
            self.insufficient_memory = True

class AnBeKoM:

    def pop_error(self, msg):
        error_dialog = Pmw.MessageDialog(self.parent, title = 'Error',message_text = msg)



    def __init__(self,app):
        parent = app.root
        self.parent = parent
        # workaround for list binding
        self.configJustLoaded = 0
        #by default select all
        self.xButton = "empty"

        self.dataStructure = DataStruct()

        self.optimizeNearValue = StringVar()
        self.optimizeNearValue.set("4.0")
        self.optimizeRadius = StringVar()
        self.optimizeRadius.set("1.8")
        self.AAKEY = "20_AA"
        self.inputsSubdir = "inputs"
        #ignore structures which match the follwing regexps
        self.ignoreStructures = [r"^origins$",r"_origins$", r"_v_origins$", r"_t\d\d\d_\d$"]

        # Create the dialog.
        self.dialog = Pmw.Dialog(parent,
                                 buttons = (defaults["compute_command"], defaults["exit_command"]),
                                 #defaultbutton = 'Run CAVER',
                                 title = 'Caver ' + VERSION,
                                 command = self.execute)
        self.dialog.withdraw()
        lbb = "Caver %s" % (VERSION,)

        w = tk.Label(self.dialog.interior(),
                                text = lbb ,
                                background = 'orange',
                                foreground = 'white',
                                #padx = 100,
                                )
        #w.config(font=labelfont)
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
        ww = tk.Button(self.dialog.interior(), text = 'Help and how to cite', command = self.launchHelp)
        ww.pack()


        #self.stdam_list = [ 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX', 'CYX', 'GLX', 'HI0', 'HID', 'HIE', 'HIM', 'HIP', 'MSE', 'ACE', 'ASH', 'CYM', 'GLH', 'LYN', 'NME']
        self.stdam_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        self.transparent_list = ['HOH', 'H20', 'WAT', 'DOD', 'D20']
        #transparent_list = ['HOH', 'H20', 'WAT', 'TIP', 'SOL', 'DOD', 'D20', 'EOH', 'MOH', 'PER', 'PO4', 'SO4', 'SUL', 'IOD', 'CL']

        def quickFileValidation(s):
            if s == '': return Pmw.PARTIAL
            elif os.path.isfile(s): return Pmw.OK
            elif os.path.exists(s): return Pmw.PARTIAL
            else: return Pmw.PARTIAL

        self.caver3locationAbsolute = CAVER3_LOCATION
        # hide the location field, not sure whether this is a good step, should not be visible at least read-only?
        if (0):
            self.caver3location = Pmw.EntryField(self.dialog.interior(),
                                         labelpos='w',
                                         value = CAVER3_LOCATION,
                                         label_text = LABEL_TEXT)


            self.caver3location.pack(fill='x',padx=4,pady=1) # vertical
#win/linux
        self.binlocation = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = OUTPUT_LOCATION,
                                     label_text = 'Output directories:')


        self.binlocation.pack(fill='x',padx=4,pady=1) # vertical
        self.configgroup = Pmw.Group(self.dialog.interior(), tag_text='Configuration save/load')
        self.conflocationDefault = os.path.join(self.caver3locationAbsolute,"config.txt")
        self.DEFCONF = "(default config used)"
        self.conflocation = tk.Label(self.configgroup.interior(),text = self.DEFCONF)
        #self.conflocation = Pmw.EntryField(self.configgroup.interior(),
        #                             labelpos='w',
        #                             value = self.DEFCONF,
        #                             label_text = 'Config location:')
        self.conflocation.pack(side=LEFT, padx=4,pady=1)

        self.configfilesave = tk.Button(self.configgroup.interior(), text = 'Save settings', command = self.configout)
        self.configfilesave.pack(side=RIGHT,padx=4,pady=1)
        self.configfileload = tk.Button(self.configgroup.interior(), text = 'Load settings', command = self.configin)
        self.configfileload.pack(side=RIGHT,padx=4,pady=1)

        self.configgroup.pack(expand="yes", fill="x")

        self.javaHeap = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = defaults["default_java_heap"],
                                     label_text = 'Maximum Java heap size (MB):')
        self.javaHeap.pack(fill='x',padx=4,pady=1) # vertical
        self.tunnelsProbe = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = defaults["default_tunnels_probe"],
                                     label_text = 'Minimum probe radius:')


        self.tunnelsProbe.pack(fill='x',padx=4,pady=1) # vertical
        self.shellDepth = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = defaults["default_shell_depth"],
                                     label_text = 'Shell depth:')


        self.shellDepth.pack(fill='x',padx=4,pady=1) # vertical
        self.shellRadius = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = defaults["default_shell_radius"],
                                     label_text = 'Shell radius:')


        self.shellRadius.pack(fill='x',padx=4,pady=1) # vertical
        self.clusteringThreshold = Pmw.EntryField(self.dialog.interior(),
                                     labelpos='w',
                                     value = defaults["default_clustering_threshold"],
                                     label_text = 'Clustering threshold:')


        self.clusteringThreshold.pack(fill='x',padx=4,pady=1) # vertical
        self.approxLbl = Label(self.dialog.interior(), text="Number of approximating balls:")
        self.approxLbl.pack()
        self.approxVar = StringVar()
        self.approxVar.set("4") #default value
        self.approxSph = OptionMenu(self.dialog.interior(), self.approxVar, "4", "6", "8", "12", "20")
        #self.approxVar.set(DEFAULTVALUE_OPTION)
        self.approxSph.pack()

        labframe0 = tk.Frame(self.dialog.interior())
        labframe0.pack(fill='x',padx=4,pady=2)

        self.varremovewater = IntVar()

        self.removewaterbutton = Checkbutton(labframe0, text="Ignore waters", variable=self.varremovewater)
        self.varremovewater.set(1)

        self.inModelGroup = Pmw.Group(self.dialog.interior(), tag_text='Input model:')
        self.listbox1 = tk.Listbox(self.inModelGroup.interior(), width=25, height=6,exportselection=0)
        self.listbox1.bind('<<ListboxSelect>>',self.inputAnalyseWrap)
        yscroll1 = tk.Scrollbar(self.inModelGroup.interior(),command=self.listbox1.yview, orient=tk.VERTICAL)
        self.listbox1.pack(side=LEFT)
        yscroll1.pack(side=LEFT, fill='y')
        self.listbox1.configure(yscrollcommand=yscroll1.set)
        self.reloadListButton = tk.Button(self.inModelGroup.interior(), text = 'Reload', command = self.updateList)
        self.reloadListButton.pack(side=LEFT)
        self.inModelGroup.pack()




        self.filterGroup = Pmw.Group(self.dialog.interior(), tag_text='Input atoms:')
        self.filterGroup.pack()
        self.checklist = []
        self.buttonlist = []

        self.updateList()
        #fill with data
        #self.listbox1.insert(0,"all")
        #self.listbox1.selection_set(0, 0) # Default sel
        #tindex = 1
        #for item in cmd.get_object_list():
        #  self.listbox1.insert(tindex,str(item))
        #  tindex = tindex + 1

        self.s = dict()
        self.s[self.AAKEY] = IntVar()
        #print("reinitialise&inputAnalyse")
        #self.reinitialise()
        #initialise should be done after config load




        groupstart = Pmw.Group(self.dialog.interior(),tag_text='Starting point')

        self.surroundingsvar = tk.IntVar()

        radioframe = tk.Frame(groupstart.interior())
        group1 = Pmw.Group(radioframe,
                tag_text='Convert surroundings to x,y,x coordinates of starting point')

        group1.pack(side='top',expand = 'yes',fill='x')
        self.selectionlist = Pmw.EntryField(group1.interior(),
                                  labelpos='w',
                                  label_text='Specify selection: ',
                                  value=defaults['surroundings'],
                                  entry_width=50
                                  )
        self.selectionlist.pack(fill='x',expand='yes',padx=4,pady=1) # vertical

        self.convertButton = tk.Button(group1.interior(), text = 'Convert to x,y,z', command = self.convert)
        self.convertButton.pack(fill='x',expand='yes',padx=4,pady=1)

        group2 = Pmw.Group(radioframe,
                tag_text='x, y, z coordinates of starting point')
        group2.pack(fill = 'x', expand = 1, side='top')
        radioframe.pack(side='left',expand='yes',fill='x')
#-------------
        groupstart.pack(padx=4,pady=1,expand='yes',fill='x')

        self.xlocvar=DoubleVar()
        self.xlocvar.set(float(defaults["startingpoint"][0]))
        self.ylocvar=DoubleVar()
        self.ylocvar.set(float(defaults["startingpoint"][1]))
        self.zlocvar=DoubleVar()
        self.zlocvar.set(float(defaults["startingpoint"][2]))

        self.xlocfr = tk.Frame(group2.interior())
        labX = Label(self.xlocfr,text="x")
        self.xlocation = Entry(self.xlocfr,textvariable=self.xlocvar,width=10)
        self.scrX=Scrollbar(self.xlocfr,orient="horizontal",command=self.changeValueX)

        self.ylocfr = tk.Frame(group2.interior())
        labY = Label(self.ylocfr,text="y")
        self.ylocation = Entry(self.ylocfr,textvariable=self.ylocvar,width=10)
        self.scrY=Scrollbar(self.ylocfr,orient="horizontal",command=self.changeValueY)

        self.zlocfr = tk.Frame(group2.interior())
        labZ = Label(self.zlocfr,text="z")
        self.zlocation = Entry(self.zlocfr,textvariable=self.zlocvar,width=10)
        self.scrZ=Scrollbar(self.zlocfr,orient="horizontal",command=self.changeValueZ)

        labX.pack(side=LEFT)
        self.xlocation.pack(side=LEFT)
        self.scrX.pack(side=LEFT)
        self.xlocfr.pack(side=LEFT,fill='x',padx=4,pady=1) # vertical
        labY.pack(side=LEFT)
        self.ylocation.pack(side=LEFT)
        self.scrY.pack(side=LEFT)
        self.ylocfr.pack(side=LEFT,fill='x',padx=4,pady=1) # vertical
        labZ.pack(side=LEFT)
        self.zlocation.pack(side=LEFT)
        self.scrZ.pack(side=LEFT)
        self.zlocfr.pack(side=LEFT,fill='x',padx=4,pady=1) # vertical

        self.OpGroup = Pmw.Group(radioframe,tag_text = "Starting point optimization")
        self.OpGroup.pack(fill='x')
        self.optimizeLabel = tk.Label(self.OpGroup.interior(),text = 'Maximum distance (A): ')
        self.optimizeLabel.pack(side=LEFT)

        self.optimizeNear = tk.Entry(self.OpGroup.interior(),textvariable=self.optimizeNearValue,justify='right', width=10)
        self.optimizeNear.pack(side=LEFT,padx=4,pady=1)
        self.optimizeLabel2 = tk.Label(self.OpGroup.interior(),text="Desired radius (A):")
        self.optimizeLabel2.pack(side=LEFT, padx=0, pady=1)
        self.optimizeNear = tk.Entry(self.OpGroup.interior(),textvariable=self.optimizeRadius,justify='right', width=10)
        self.optimizeNear.pack(side=LEFT,padx=4,pady=1)
        #self.optimizeButton = tk.Button(self.OpGroup.interior(), text = 'Optimize', command = self.optimize)
        #self.optimizeButton.pack(side=LEFT,padx=5,pady=1)
        #self.UoptimizeButton = tk.Button(self.OpGroup.interior(), text = 'Undo', command = self.uoptimize)
        #self.UoptimizeButton.pack(side=LEFT,padx=1,pady=1)
        self.egroup = Pmw.Group(self.dialog.interior(),tag_text = "Computation result")
        self.egroup.pack(fill='x')
        self.aftercomp = tk.Label(self.egroup.interior(),text="test",justify='right')
        self.aftercomp.pack(side=LEFT,padx=4,pady=1)
        self.afterbutt = tk.Button(self.egroup.interior(), text='Details', command=self.details, width = 5)
        self.afterbutt.pack(side=RIGHT,padx=4,pady=1)
        self.afterbutt.config(state=DISABLED)
    #hide group for now
        self.egroup.pack_forget()

        cf = self.getConfLoc()
        self.configLoad(cf)

        self.inputAnalyse()
        self.showAppModal()
    def getConfLoc(self):
        cf = self.conflocation.cget("text")
        if cf == self.DEFCONF:
            return self.conflocationDefault
        else:
            return cf
    def showCrisscross(self):
        startpoint=(float(self.xlocvar.get()),float(self.ylocvar.get()),float(self.zlocvar.get()))
        cmd.delete("crisscross")
        self.crisscross(startpoint[0],startpoint[1],startpoint[2],0.5,"crisscross")

#win/linux
    def changeValueX(self, *args):
            a = args[0] if len(args) == 1 else args[1]
            val=float(self.xlocvar.get())+float(a)*0.2
            self.xlocvar.set(val)
            self.showCrisscross()
    def changeValueY(self, *args):
            a = args[0] if len(args) == 1 else args[1]
            val=float(self.ylocvar.get())+float(a)*0.2
            self.ylocvar.set(val)
            self.showCrisscross()
    def changeValueZ(self, *args):
            a = args[0] if len(args) == 1 else args[1]
            val=float(self.zlocvar.get())+float(a)*0.2
            self.zlocvar.set(val)
            self.showCrisscross()



    def showAppModal(self):
        self.dialog.show()

    def structureIgnored(self, name):
        for key in self.ignoreStructures:
            if re.search(key, name):
                return 1
        return 0
    def updateList(self):
        #print("updateList")
        self.listbox1.delete(0, tk.END)
        #fill with data
        self.listbox1.selection_set(0, 0) # Default sel
        tindex = 0
        for item in cmd.get_object_list():
            stri = str(item)
            if not self.structureIgnored(stri):
                self.listbox1.insert(tindex,str(item))
                tindex = tindex + 1
        #select first by default
        self.listbox1.select_set(0)
        self.inputAnalyse()

    def launchHelp(self):
        import webbrowser
        webbrowser.open(url)

    def details(self):
        fc = self.loadFileContent("%s/warnings.txt" % (self.out_dir))
        error_dialog = Pmw.MessageDialog(self.parent,title = 'Information', message_text = fc,)

    def loadFileContent(self, file):
        handler = open(file)
        lines = handler.readlines()
        wresult = ""
        for line in lines:
            wresult += line
        return wresult

    def suitable(dir):
        if not os.path.exists(dir):
            return True
        else:
            if os.listdir(dir) == []:
                return True
            else:
                return False

    def initialize_out_dir(self):
        dir = self.binlocation.getvalue()
        if not os.path.exists(dir):
            os.mkdir(dir)
        dir = dir.replace("\\","/")
        if (dir.endswith("/")):
            dir = dir[:-1]
        out_home = dir + "/caver_output/"
        if not os.path.exists(out_home):
            os.mkdir(out_home)

        max = 0
        ls = os.listdir(out_home)
        for f in ls:
            fn = os.path.basename(f)
            if fn.isdigit():
                i = int(fn)
                if max < i:
                    max = i

        new_dir = out_home + str(max + 1)
        self.CreateDirectory(new_dir)
        self.out_dir = new_dir
        print("Output will be stored in " + self.out_dir)

    def coordinatesNotSet(self):
        b = float(self.xlocvar.get()) == 0 and float(self.ylocvar.get()) == 0 and float(self.zlocvar.get()) == 0
        return b

    def printErrorMessages(self, dir):
        f = dir + '/messages.txt'
        m = ""
        if os.path.exists(f):
            if 0 < os.path.getsize(f):
                handler = open(f)
                lines = handler.readlines()
                for line in lines:
                    m = m + line + "\n"
                handler.close()
                self.pop_error(m)


    def execute(self, result):
        #elif result == defaults["warn_command"]:
            #self.wtext = tk.Text(root, height=26, width=50)
            #scroll = Scrollbar(root, command=text.yview)
            #text.configure(yscrollcommand=scroll.set)
            #handler = open(self.caver3locationAbsolute + "/out/warnings.txt")
            #lines = handler.readlines()
            #wresult = ""
            #for line in lines:
            #  wresult += line
            #error_dialog = Pmw.MessageDialog(self.parent,title = 'Information', message_text = wresult,)
        if result == defaults["compute_command"]:

            if self.coordinatesNotSet():
                self.pop_error("Please specify starting point - e.g. by selecting atoms or residues and clicking at the button 'Convert to x, y, z'.")
                return


            self.showCrisscross()

            #input
            sel1index = self.listbox1.curselection()[0]
            sel1text = self.listbox1.get(sel1index)


            self.whichModelSelect = sel1text

            #print('selected ' + self.whichModelSelect)
            sel=cmd.get_model(self.whichModelSelect)

            self.initialize_out_dir()

            # create subdirectory for inputs
            outdirInputs = self.out_dir + "/" + self.inputsSubdir
            self.CreateDirectory(outdirInputs)

            self.stdamString = "+".join(self.stdam_list)
            # jen to zaskrtnute
            generatedString = ""
            for key in self.s:
                if self.s[key].get() == 1:
                # pak pouzit do vyberu:
                    if key == self.AAKEY:
                        generatedString = generatedString + "+" + self.stdamString
                    else:
                        generatedString = generatedString + "+" + key

            generatedString = generatedString[1:]
            #print("Checked: " + generatedString)

            mmodel = cmd.get_model(self.whichModelSelect)
            #print(self.whichModelSelect + " asize: " + str(len(mmodel.atom)))
            #newmodel = Indexed()
            #for matom in mmodel.atom:
                #if generatedString.find(matom.resn) > -1:
                    #print(matom.resn)
                    #newmodel.atom.append(matom)


            #cmd.load_model(newmodel,"tmpCaverModel")
            #cmd.label("example","name")

            input = "%s/%s.pdb" % (outdirInputs, self.whichModelSelect)
            cmd.set('retain_order',1)
            cmd.sort()
            cmd.save(input, self.whichModelSelect) # to by ulozilo cely model whichModelSelect.
            #cmd.save(input, "tmpCaverModel")

            #cmd.delete("tmpCaverModel")

            cesta = os.getcwd()


            # set ignore waters to false -- the model is already filtered by input model and aminos
            self.varremovewater.set(0)

            caverfolder = "%s" % (self.caver3locationAbsolute)
            caverjar = caverfolder + "/" + "caver.jar"

            cfg = self.getConfLoc()

            # create new config
            cfgTimestamp = time.strftime("%Y-%m-%d-%H-%M")
            cfgnew = outdirInputs + "/config_" + cfgTimestamp + ".txt"
            self.configSave(cfgnew, cfg)

            # set correct java options
            #javaOpts = JOPTS.replace("@", self.javaHeap.getvalue())

            pj = PyJava(self.javaHeap.getvalue(), caverfolder, caverjar, outdirInputs, cfgnew, self.out_dir)
            if pj.java_missing:
                return

            pj.run_caver()

            if pj.insufficient_memory:
                self.pop_error("Available memory (" + str(pj.xmx) + " MB) is not sufficient to analyze this structure. Try to allocate more memory. 64-bit operating system and Java are needed to get over 1200 MB. Using smaller 'Number of approximating balls' can also help, but at the cost of decreased accuracy of computation.")

            self.printErrorMessages(self.out_dir)
            prevDir = os.getcwd()
            print(prevDir)

            runview = "run " + self.out_dir + "/pymol/view_plugin.py"
            print(runview)
            cmd.do(runview)
            # adjust gui to display warnings & group
            self.egroup.pack(fill="x")

            err = "%s/warnings.txt" % (self.out_dir)
            if os.path.exists(err) and os.stat(err)[6] == 0:
                self.aftercomp.config(text="Computation finished succesfully")
                self.afterbutt.config(state=DISABLED)
            else:
                self.aftercomp.config(text="Warnings detected during computation")
                self.afterbutt.config(state=ACTIVE)

            #pass
            #self.deleteTemporaryFiles()
        else:
            #
            # Doing it this way takes care of clicking on the x in the top of the
            # window, which as result set to None.
            #
            if __name__ == '__main__':
                    #
                    # dies with traceback, but who cares
                    #
                self.parent.destroy()
            else:
                #self.dialog.deactivate(result)
                #global CAVER_BINARY_LOCATION
                #CAVER_BINARY_LOCATION = self.out_dir
                self.dialog.withdraw()

    def CreateDirectory(self,dir):
        if os.path.isdir(dir):
            return
        parent, base = os.path.split(dir)
        self.CreateDirectory(parent)
        os.mkdir(dir)
    # simple fix precision of a string or a number (conversion)
    def fixPrecision(self, numberStr):
        return math.floor(float(numberStr) * 1000) / 1000
    def convert(self):

        sel=cmd.get_model('(all)')
        cnt=0
        for a in sel.atom:
            cnt+=1
        if cnt == 0:
            error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',message_text = 'ERROR: No molecule loaded.',)
        #try:
        if 1:
            startpoint=[]
            s = self.selectionlist.getvalue()
            #startpoint=self.computecenter(s)
            startpoint = self.compute_center(s)
            if None == startpoint:
                return
            self.xlocvar.set(self.fixPrecision(startpoint[0]))
            self.ylocvar.set(self.fixPrecision(startpoint[1]))
            self.zlocvar.set(self.fixPrecision(startpoint[2]))
            self.crisscross(startpoint[0],startpoint[1],startpoint[2],0.5,"crisscross")
            self.showCrisscross()
        #except:
        #    error_dialog = Pmw.MessageDialog(self.parent,title = 'Error',message_text = 'ERROR: Invalid selection name',)

    def containsValue(self, array, value):
        for v in array:
            if (v == value):
                return 1
        return 0

    def configin(self):
        indi = os.path.dirname(self.getConfLoc())
        filepath = filedialog.askopenfilename(title="Open config file", initialdir=indi, filetypes=[("config txt file","*.txt"), ("all files","*.*")])
        if not filepath: return
        self.conflocation.config(text=filepath)
        self.configLoad(self.getConfLoc())
        self.configJustLoaded = 1
    def configout(self):
        indi = os.path.dirname(self.getConfLoc())
        filepath = filedialog.asksaveasfilename(title="Save config file", initialdir=indi,filetypes=[("config txt file","*.txt"), ("all files","*.*")], defaultextension='.txt')
        if not filepath: return
        self.configSave(filepath, self.getConfLoc())
    #perform actual config parse here
    def configLoad(self, file):
        self.dataStructure.clear()
        print('cleared datastruct')
        self.clearGUI()
        #Pmw.MessageDialog(self.parent,title = 'Information',message_text = file)
        # do nothing if file not exists
        if not os.path.isfile(file):
            return

        handler = open(file)
        lines = handler.readlines()
        for line in lines:
            liner = line.strip()
            # remove everything after last occurence of # char
            if '#' in liner and not liner.startswith('#'):
                liner = liner[0:liner.rfind("#")-1]
            if liner.startswith('#'):
                self.dataStructure.add("#", liner, 1)
            elif liner == "":
                self.dataStructure.add("#EMPTY#", liner, 1)
            else:
                parsed = liner.split(' ')
                key = parsed[0]
                if len(parsed) <= 1:
                    print('skipping ' + key)
                    val = ""
                else:
                    val = " ".join(parsed[1:len(parsed)])
                    self.dataStructure.add(key, val, 0)

        # check for specific problematic definitions in config


        if self.dataStructure.indexOf('starting_point_coordinates') != -1 and self.dataStructure.indexOf('starting_point_atom') != -1:
            Pmw.MessageDialog(self.parent,title = 'Information',message_text = 'Simultaneous usage of starting_point_coordinates parameter with starting_point_atom parameters is not supported by plugin. Please, use only one of these parameters. Now ignoring atom.')
            self.dataStructure.remove('starting_point_atom')
        if self.dataStructure.indexOf('starting_point_coordinates') != -1 and self.dataStructure.indexOf('starting_point_residue') != -1:
            Pmw.MessageDialog(self.parent,title = 'Information',message_text = 'Simultaneous usage of starting_point_coordinates parameter with starting_point_residue parameters is not supported by plugin. Please, use only one of these parameters. Now ignoring residue.')
            self.dataStructure.remove('starting_point_residue')
        if self.dataStructure.indexOf('starting_point_coordinates') == -1:
            #perform harakiri with selecting model and pre-loading coordinates with the command similar to the one below
            #cmd.select('starting_point','id 573+658 & structure | resi 120+24 & structure')
            selector = []
            rids = ""
            aids = ""
            if self.dataStructure.indexOf('starting_point_residue') != -1:
               rids = "+".join(self.dataStructure.get('starting_point_residue').split(" "))
            if self.dataStructure.indexOf('starting_point_atom') != -1:
               aids = "+".join(self.dataStructure.get('starting_point_atom').split(" "))
            #print(aids)
            #print(rids)
            sel1index = self.listbox1.curselection()
            if (sel1index):
                sel1text = self.listbox1.get(sel1index[0])
                model = self.listbox1.get(sel1index)
                if (aids):
                    selector.append("id " + aids + " & " + model)
                if (rids):
                    selector.append("resi " + rids + " & " + model)

                if len(selector) > 0:
                    selectorStr = " | ".join(selector)
                    print("updating starting point with: " + selectorStr)
                    cmd.select('starting_point', selectorStr)
                    # set field value
                    self.selectionlist.setvalue('starting_point')
                    # call conversion to xyz
                    self.convert()

        # test include/exclude
        if self.hasIncludeExclude():
           Pmw.MessageDialog(self.parent,title = 'Information',message_text = 'include_ and exclude_ parameters are not supported by plugin. Please, use the plugin to specify residues to be analyzed.')

        #print("reading done...")
        #now, all read in the structure. Multi-line params merged into one-liners
        #Traverse the structure and update gui controls
        self.structureLoad()
    def clearGUI(self):
        self.javaHeap.setvalue(defaults["default_java_heap"])
        self.tunnelsProbe.setvalue("")
        self.shellDepth.setvalue("")
        self.shellRadius.setvalue("")
        self.clusteringThreshold.setvalue("")
        self.approxVar.set("")
        self.optimizeNearValue.set("")
        self.optimizeRadius.set("")
        self.xlocvar.set(0)
        self.ylocvar.set(0)
        self.zlocvar.set(0)
    def hasIncludeExclude(self):
        notAllowed = [ "include_residue_names", "include_residue_ids", "include_atom_numbers", "exclude_residue_names", "exclude_residue_ids", "exclude_atom_numbers"]
        sk = self.dataStructure.getKeys()
        for i in range (0, len(sk)):
            key = sk[i]
            if key != "include *" and key in notAllowed:
                return 1
        return 0

    #consider all properties in the gui and store them into config file supplied
    # load file "readfile" and store params into new config file "file"
    def configSave(self, file, readfile):
        #load config to structure and then replace with gui params
        self.dataStructure.clear()
        handler = open(readfile)
        lines = handler.readlines()
        for line in lines:
            liner = line.strip()
            # remove everything after last occurence of # char
            if '#' in liner and not liner.startswith('#'):
                liner = liner[0:liner.rfind("#")-1]
            if liner.startswith('#'):
                self.dataStructure.add("#", liner, 1)
            elif liner == "":
                self.dataStructure.add("#EMPTY#", liner, 1)
            else:
                parsed = liner.split(' ')
                key = parsed[0]
                val = " ".join(parsed[1:len(parsed)])
                self.dataStructure.add(key, val, 0)
        handler.close()
        #print("reading done...")
        #now, all read in the structure. Multi-line params merged into one-liners
        #Traverse the structure and update gui controls

        self.structureUpdateFromGui()
        # now, gui is sync with data structure
        f = open(file, 'w')
        keys = self.dataStructure.getKeys()
        values = self.dataStructure.getValues()
        for idx in range (0, len(keys)):
            key = keys[idx]
            value = values[idx]
            #print("saving value/key " + key + " " + value + " " + str(len(keys)))
            if value == "":
                # do nothing for empty values (do not save to config, caver uses default values
                noop = 1 #noop
            elif key == "#":
                f.write(value)
                f.write("\n")
            elif key == "REMOVED":
                # do nothing
                noop = 1 #noop
            elif key == "#EMPTY#":
                f.write("\n")
            else:
                f.write(key)
                f.write(" ")
                f.write(value)
                f.write("\n")
        f.close()
    def structureLoad(self):
        keys = self.dataStructure.getKeys()
        values = self.dataStructure.getValues()
        for i in range(0, len(keys)):
            key = keys[i]
            val = values[i]
            #print(key + "->" + val)
            if key == "probe_radius":
                self.tunnelsProbe.setvalue(str(val))
            elif key == "java_heap":
                self.javaHeap.setvalue(str(val))
            elif key == "shell_depth":
                self.shellDepth.setvalue(str(val))
            elif key == "shell_radius":
                self.shellRadius.setvalue(str(val))
            elif key == "clustering_threshold":
                self.clusteringThreshold.setvalue(str(val))
            elif key == "number_of_approximating_balls":
                self.approxVar.set(str(val))
            elif key == "max_distance":
                self.optimizeNearValue.set(str(val))
            elif key == "desired_radius":
                self.optimizeRadius.set(str(val))
            elif key == "include_residue_names":
                self.s = dict()
                self.s.clear()
                ress = (str(val)).split(" ")
                aa_added = 0
                for idx in range(0, len(ress)):
                    fromconf = ress[idx]
                    if fromconf in self.stdam_list:
                        if aa_added == 0:
                            self.s[self.AAKEY] = IntVar()
                            self.s[self.AAKEY].set(1)
                            aa_added = 1
                    elif fromconf.strip() != "":
                        self.s[fromconf] = IntVar()
                        self.s[fromconf].set(1)
                self.reinitialiseFromConfig()
            elif key == "starting_point_coordinates":
                starr = (str(val)).split(" ")
                self.xlocvar.set(float(self.fixPrecision(starr[0])))
                self.ylocvar.set(float(self.fixPrecision(starr[1])))
                self.zlocvar.set(float(self.fixPrecision(starr[2])))
    def structureUpdateFromGui(self):
        self.dataStructure.replace("probe_radius", self.tunnelsProbe.getvalue(), 0)
        self.dataStructure.replace("java_heap", self.javaHeap.getvalue(), 0)
        self.dataStructure.replace("shell_depth",self.shellDepth.getvalue(), 0)
        self.dataStructure.replace("shell_radius",self.shellRadius.getvalue(), 0)
        self.dataStructure.replace("clustering_threshold",self.clusteringThreshold.getvalue(), 0)
        self.dataStructure.replace("number_of_approximating_balls",self.approxVar.get(), 0)
        #check-boxed residues
        result = ""
        for item in self.s.keys():
            #print("ITEM: " + item)
            # do not print(all amino acids, just the item 20_AA)
            #if item == self.AAKEY and self.s[item].get() == 1:
            #    result = result + " " + string.join(self.stdam_list, " ")
            #elif self.s[item].get() == 1:
            if self.s[item].get() == 1:
                result = result + " " + item
        self.dataStructure.replace("include_residue_names", result, 0)

        #active site:
        #remove other starting point definitions except those with atoms
        self.dataStructure.remove("starting_point_residue")
        self.dataStructure.remove("starting_point_atom")

        asit = str(self.xlocvar.get()) + " " + str(self.ylocvar.get()) + " " + str(self.zlocvar.get())
        self.dataStructure.replace("starting_point_coordinates",asit, 0)
        self.dataStructure.replace("max_distance",self.optimizeNearValue.get(), 0)
        self.dataStructure.replace("desired_radius",self.optimizeRadius.get(), 0)
        #print("len" + str(len(self.dataStructure.getKeys()))  + str(len(self.dataStructure.getValues())))
    def stdamMessage(self):
        Pmw.MessageDialog(self.parent,title = 'Information',message_text = self.AAKEY + ': Standard amino acids: \n ' + ", ".join(self.stdam_list))

    def inputAnalyseWrap(self, args):
            #print(self.listbox1.curselection()[0] # aby to fungovalo, musi byt bindnute na <<ListboxSelect>>)
            #print("calling from wrap")
        if self.configJustLoaded == 1:
            self.configJustLoaded = 0
        else:
            self.inputAnalyse()

    def inputAnalyse(self):
        sel1list = self.listbox1.curselection()
        if sel1list:
            sel1index = sel1list[0]
            sel1text = self.listbox1.get(sel1index)
            self.whichModelSelect = sel1text
            sel=cmd.get_model(self.whichModelSelect)
        #pripravit kontrolni strukturu pro nalezene
        self.s = dict()
        self.s.clear()
        if sel1list:
            #cntr = 0
            for a in sel.atom:
                if not a.resn in self.s:
                    if (self.containsValue(self.stdam_list, a.resn)):
                        self.s[self.AAKEY] = IntVar()
                        self.s[self.AAKEY].set(1)
                    else:
                        self.s[a.resn] = IntVar()
                        if (self.containsValue(self.transparent_list, a.resn)):
                            # uncheck waters and such
                            self.s[a.resn].set(0)
                        else:
                            self.s[a.resn].set(1)
        self.reinitialise()

    def reinitialiseFromConfig(self):
        ksorted = sorted(self.s.keys())
        #print("calling initialise from config" + str(len(ksorted)))
        for xs in self.checklist:
            xs.grid_remove()
        self.checklist = []

        for xs in self.buttonlist:
            xs.grid_remove()
        self.buttonlist = []
        cntr = 0

        if self.AAKEY in ksorted:
            tmpButton = tk.Checkbutton(self.filterGroup.interior(), text=self.AAKEY, variable=self.s[self.AAKEY])
            tmpButton.var = self.s[self.AAKEY]
            tmpButton.grid(sticky=W, row = int(cntr/5), column = (cntr % 5))
            self.checklist.append(tmpButton)
            cntr = cntr + 1
            tmpButton = tk.Button(self.filterGroup.interior(), text='?', command=self.stdamMessage, width = 5)
            tmpButton.grid(sticky=W, row = 0, column=1) # 0,1 = stdam, 0,2 = help
            tmpButton.var = self.s[self.AAKEY]
            tmpButton.grid(sticky=W, row = int(cntr/5), column = (cntr % 5))
            self.checklist.append(tmpButton)
            cntr = cntr + 4
        for key in ksorted:
            if key != self.AAKEY:
                    #print("adding button" + key)
                tmpButton = tk.Checkbutton(self.filterGroup.interior(), text=key, variable=self.s[key])
                tmpButton.var = self.s[key]
                tmpButton.grid(sticky=W, row = int(cntr/5), column = (cntr % 5))
                self.checklist.append(tmpButton)
                cntr = cntr + 1
    def reinitialise(self):
        #if 1: return
        #TODO: pouzivat zde uz setrideny
        ksorted = sorted(self.s.keys())

        #print("calling initialise")
        for xs in self.checklist:
            xs.grid_remove()
        self.checklist = []

        for xs in self.buttonlist:
            xs.grid_remove()
        self.buttonlist = []

        cntr = 0
        # tady uz setrideny, se STDAM a STDRNA na zacatku
        for key in ksorted:
            #
            #if cntr == 1:
            #    cntr = cntr + 4
            tmpButton = tk.Checkbutton(self.filterGroup.interior(), text=key, variable=self.s[key])
            tmpButton.var = self.s[key]
            tmpButton.grid(sticky=W, row = int(cntr/5), column = (cntr % 5))
            self.checklist.append(tmpButton)

                # zaridit balloon help -- neni potreba kdyz je tlacitko
                #if key == "AA":
                #  balloon = Pmw.Balloon(self.parent)
                #  balloon.bind(tmpButton, 'STanDard AMino acids: \n ' + string.join(self.stdam_list, ", "), 'STanDard AMino acids')

            # kdyz je tam pridano STDAM, vlozit tam  tedy i napovedu
            if key == self.AAKEY:
                self.xButton = tk.Button(self.filterGroup.interior(), text='?', command=self.stdamMessage, width = 5)
                self.xButton.grid(sticky=W, row = 0, column=1) # 0,1 = stdam, 0,2 = help
                self.buttonlist.append(self.xButton)
                cntr = cntr + 4


            cntr = cntr + 1


        ##print("size: " + str(len(self.checklist)))
        #aButton = tk.Button(self.filterGroup.interior(), text='Reload', command=self.inputAnalyse)



        ## zarovnat
        #if not cntr % 3 == 0:
        #  cntr = cntr + 3 - (cntr % 3)
        ##analyse,save,load
        #aButton.grid(row = int(cntr/3), column = (cntr % 3))
        #self.buttonlist.append(aButton)

    def getAtoms(self, selection="(all)"):
        return cmd.identify(selection, 0)

    def getResids(self, selection="(all)"):
        stored.list=[]
        cmd.iterate(selection,"stored.list.append((resi,chain))")
        return set(stored.list)


    def getObjectName(self, selection="(all)"):
        pairs = cmd.identify(selection, 1)
        name = None
        names = set([])
        for p in pairs:
            names.add(p[0])
        if 0 == len(names):
            self.pop_error("Selection is empty.")
        elif 1 == len(names):
            name = names.pop()
        else:
            s = "Starting point selection need to be limited to one object. Currently, it includes these objects: "
            for n in names:
                s += n + ' '
            self.pop_error(s)
        return name

    def compute_center(self,selection="(all)"):
        if not selection in cmd.get_names("selections") and not selection in cmd.get_names("objects"):
            self.pop_error("Selection '" + selection + "' does not exist, using all atoms.")
            selection = "all"
        object = self.getObjectName(selection)
        if None == object:
            return None
        Ts = []
        residues = self.getResids(selection) # SET1
        atoms = self.getAtoms(selection)     # SET2
        for r in residues:
            r_sel = 'resi ' + str(r[0]) + ' and chain ' + r[1] + ' and object ' + object
            residue_atoms = self.getAtoms(r_sel)
            all = []
            for a in residue_atoms:
                if a in atoms:
                    all = all + [a]
            if len(all) == len(residue_atoms):
                Ts = Ts + [self.computecenterRA(r_sel)]
            else:
                for a in all:
                    Ts = Ts + [self.computecenterRA('id ' + str(a)  + ' and object ' + object)]

        print('Centers: %s' % ', '.join(map(str, Ts)))
        #print('Centers: ' + Ts)
        sumx = 0
        sumy = 0
        sumz = 0
        if len(Ts) == 0:
            return (0, 0, 0)
        l = len(Ts)
        for center in Ts:
            sumx += center[0]
            sumy += center[1]
            sumz += center[2]
        print('Starting point: ' + str(sumx) + " " + str(sumy) + " " + str(sumz) + " " + str(l))
        return (sumx/l, sumy/l, sumz/l)
        # Ts = []
        # detect all residues in selection => SET1
        # detect all atoms in selection => SET2
        # foreach residue, detect all its atoms
        #    if all atoms are in SET2, attach computeCenter(S1) to Ts
        #    if at least one atom is not in SET1, traverse these atoms again and !only! those in SET2 add to Ts
        # computeCenter(Ts)
        # return

    # compute center for given selection
    def computecenterRA(self,selection="(all)"):
        stored.xyz = []
        cmd.iterate_state(1,selection,"stored.xyz.append([x,y,z])")
        centx=0
        centy=0
        centz=0
        cnt=0
        for a in stored.xyz:
            centx+=a[0]
            centy+=a[1]
            centz+=a[2]
            cnt+=1
        centx/=cnt
        centy/=cnt
        centz/=cnt
        return (centx,centy,centz)

    def computecenter(self,selection="(all)"):
        gcentx=0
        gcenty=0
        gcentz=0
        gcnt=0
        for selstr in selection.split():
            sel=cmd.get_model(selstr)

            centx=0
            centy=0
            centz=0
            cnt=len(sel.atom)
            if (cnt == 0):
                print('warning: selection used to compute starting point is empty')
                return (0, 0, 0)
            for a in sel.atom:
                centx+=a.coord[0]
                centy+=a.coord[1]
                centz+=a.coord[2]
            centx/=cnt
            centy/=cnt
            centz/=cnt
        #       fmttext="%lf\t%lf\t%lf\n" % (centx,centy,centz)
#               print(centx,centy,centz)
            gcentx+=centx
            gcenty+=centy
            gcentz+=centz
            gcnt+=1

        gcentx/=gcnt
        gcenty/=gcnt
        gcentz/=gcnt
        return (gcentx,gcenty,gcentz)

    def crisscross(self,x,y,z,d,name="crisscross"):

        obj = [
        LINEWIDTH, 3,

        BEGIN, LINE_STRIP,
        VERTEX, float(x-d), float(y), float(z),
        VERTEX, float(x+d), float(y), float(z),
        END,

        BEGIN, LINE_STRIP,
        VERTEX, float(x), float(y-d), float(z),
        VERTEX, float(x), float(y+d), float(z),
        END,

        BEGIN, LINE_STRIP,
        VERTEX, float(x), float(y), float(z-d),
        VERTEX, float(x), float(y), float(z+d),
        END

        ]
        view = cmd.get_view()
        cmd.load_cgo(obj,name)
        cmd.set_view(view)


# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = tk.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = AnBeKoM(app)
    exitButton = tk.Button(app.root, text = 'Exit', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()
