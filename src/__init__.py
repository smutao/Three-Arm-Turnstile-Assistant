'''
Three-Arm Turnstile Assistant 
by Yunwen Tao, Ph.D.

email: ywtao.smu[at]gmail.com 

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import os

import pymol
from pymol import cmd
from pymol.wizard import Wizard
from chempy import cpv
from pymol.cgo import *


# below is for the geometry calculation 
from math import pi ,sin, cos, sqrt

def R(theta, u):
    return [[cos(theta) + u[0]**2 * (1-cos(theta)), 
             u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta), 
             u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
            [u[0] * u[1] * (1-cos(theta)) + u[2] * sin(theta),
             cos(theta) + u[1]**2 * (1-cos(theta)),
             u[1] * u[2] * (1 - cos(theta)) - u[0] * sin(theta)],
            [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
             u[1] * u[2] * (1-cos(theta)) + u[0] * sin(theta),
             cos(theta) + u[2]**2 * (1-cos(theta))]]


def NormVec3pt(p1,p2,p3):
    x1=p1[0]
    x2=p2[0]
    x3=p3[0]
    y1=p1[1]
    y2=p2[1]
    y3=p3[1]
    z1=p1[2]
    z2=p2[2]
    z3=p3[2]
    a=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
    b=(z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
    c=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

    norm=sqrt(a*a+b*b+c*c)
    a=a/norm
    b=b/norm
    c=c/norm

    return [a,b,c]

def Rotate2(anchor,pointToRotate, point1, point2, point3, theta):

    u= []
    u=NormVec3pt(point1, point2, point3)

    r = R(theta, u)
    rotated = []

    for i in range(3):
        rotated.append((sum([r[j][i] * (pointToRotate[j]-anchor[j]) for j in range(3)])))
    for i in range(3):
        rotated[i] = rotated[i] + anchor[i]

    return rotated




# below is for the Wizard

object_prefix = "_pw"
object_subgroup_prefix = "_s"


class TurnstileWizard(Wizard):

    def __init__(self):
        Wizard.__init__(self)

        # some attributes to do with picking
        self.pick_count = 0
        self.subgroup_count = 0 # subgroup per arm 
        self.subgroup_sum = []
        self.object_count = 0
        self.object_prefix = object_prefix # put "_" can hide the explicit item from the list 
        self.object_subgroup_prefix = object_subgroup_prefix


        # the plane facet size (the 'radius' of the section of plane we show)
        self.facetSize = 5

        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect()

    def reset(self):
        cmd.delete(self.object_prefix + "*") # removed the selection of atoms "_pw*"
        #cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        self.pick_count = 0
        self.subgroup_count = 0 
        self.subgroup_sum = []

        cmd.refresh_wizard()

    def delete_all(self):
        cmd.delete("plane*")

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.reset()
        self.delete_all()

    def get_prompt(self): # print some info on the main screen
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = [ 'Please click on the anchor atom...']
        elif self.pick_count == 1:
            self.prompt = [ 'Please click on the first arm atom(s)...' ]
        elif self.pick_count == 2:
            self.prompt = [ 'Please click on the second arm atom(s)...' ]
        elif self.pick_count == 3:
            self.prompt = [ 'Please click on the third arm atom(s)...' ]    
        elif self.pick_count == 4:    
            self.prompt = [ 'Please click "Picking Finished" button.' ]   
        return self.prompt

    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
        #try:
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)
        #except pymol.CmdException, pmce:
        #    print pmce

    def pickNextAtom(self, atom_name):
        # transfer the click selection to a named selection
        cmd.select(atom_name, "(pk1)")
        print(atom_name) # pw0, pw1, ...

        # delete the click selection
        cmd.unpick()

        # using the magic of indicate, highlight stuff
        indicate_selection = "_indicate" + self.object_prefix # starting with "_" can hide from window 

        #print(indicate_selection)
        cmd.select(indicate_selection, atom_name)
        cmd.enable(indicate_selection)


        self.subgroup_count += 1

        #self.pick_count += 1
        #self.error = None

        # necessary to force update of the prompt
        #cmd.refresh_wizard()

    def finish_1arm(self):

        self.pick_count += 1
        self.error = None

        self.subgroup_sum.append(self.subgroup_count)
        #print(self.subgroup_sum)
        self.subgroup_count = 0

        # necessary to force update of the prompt
        cmd.refresh_wizard()        



    def do_pick(self, picked_bond):

        # this shouldn't actually happen if going through the "do_select"
        if picked_bond:
            self.error = "Error: please select bonds, not atoms"
            print(self.error)
            return

        atom_name = self.object_prefix + str(self.pick_count) + self.object_subgroup_prefix + str(self.subgroup_count) # combine string pw + [num]
        
        #cmd.show(spheres, atom_name)

        if self.pick_count < 4:
            self.pickNextAtom(atom_name)
            if self.pick_count == 0: # allow only one achor atom
               self.finish_1arm() 

        else:
            self.pickNextAtom(atom_name)

            #cmd.disable( atom_name)


            #point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
            #point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
            #point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))


            #plane = planeFromPoints(point1, point2, point3, self.facetSize)
            #planeName = "plane-%02d" % self.object_count
            self.object_count += 1
            #makePrimitive(plane, planeName)
            #cmd.show("cgo", "plane*")

            self.pick_count = 0
             
            print("print me ") 

            #cmd.alter_state(1,"(%s%s)" % (self.object_prefix, "0"),"x=x+1.0") # to test alter_state

            self.reset() # reset will clean up the selections....

    def get_panel(self): # show menu 
        return [
            [ 1, 'Three-Arm Turnstile Wizard',''],
            [ 2, 'Reset','cmd.get_wizard().reset()'],
            #[ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
            [ 2, 'Arm Atoms Selection Done', 'cmd.get_wizard().finish_1arm()'],
            [ 2, 'Done','cmd.set_wizard()'],
        ]









# below is for GUI 

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Turnstile Assistant', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    # entry point to PyMOL's API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt


    wiz = TurnstileWizard()
    initCoor = []





    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'demowidget.ui')
    form = loadUi(uifile, dialog)



    form.slider_angle.setDisabled(True)
    form.angle_text.setDisabled(True)
    form.set_angle.setDisabled(True)



    # callback for the "Ray" button
    def run():
        # get form data
        height = form.input_height.value()
        width = form.input_width.value()
        dpi = form.input_dpi.value()
        filename = form.input_filename.text()
        units = form.input_units.currentText()

        # calculate dots per centimeter or inch
        if units == 'cm':
            dots_per_unit = dpi * 2.54
        else:
            dots_per_unit = dpi

        # convert image size to pixels
        width *= dots_per_unit
        height *= dots_per_unit

        # render the image
        if filename:
            cmd.png(filename, width, height, dpi=dpi, ray=1, quiet=0)
        else:
            cmd.ray(width, height, quiet=0)
            print('No filename selected, only rendering on display')

    # callback for the "Browse" button
    def browse_filename():
        filename = getSaveFileNameWithExt(
            dialog, 'Save As...', filter='PNG File (*.png)')
        if filename:
            form.input_filename.setText(filename)



    ####################
    def start_wiz():

        wiz.reset()
        initCoor[:]=[]  


        #print('me clicked...')
        #wiz = TurnstileWizard()
        cmd.set_wizard(wiz)  

        #print('me clicked x2...')

        form.slider_angle.setDisabled(True)
        form.slider_angle.setValue(0.0)
        form.angle_text.setDisabled(True)
        form.set_angle.setDisabled(True)
        form.status_text.setText("Please specify the three-arm turnstile atoms")



    def picking_finish():
        
        #print(wiz.subgroup_sum)
        if len(wiz.subgroup_sum) != 4:
           form.status_text.setText("Please finish your selection first...")
        else:
           form.status_text.setText("You've selected "+str(wiz.subgroup_sum[1])+", "+str(wiz.subgroup_sum[2])+" and "+str(wiz.subgroup_sum[3])+" atoms for three arms.")    

           form.slider_angle.setEnabled(True)
           form.angle_text.setEnabled(True)
           form.set_angle.setEnabled(True)

           form.angle_text.setText(str(form.slider_angle.value()))

           
           # record the coordinates as initial state
           for i in range(4):
               initCoor.append([])

           str00 = object_prefix+"0"+object_subgroup_prefix+"0"
           initCoor[0].append(cmd.get_model(str00,1).get_coord_list()[0])
           #initCoor[0].append(cmd.get_atom_coords(object_prefix+"0"+object_subgroup_prefix+"0"))    
           for i in range(wiz.subgroup_sum[1]):
               str1i = object_prefix+"1"+object_subgroup_prefix+str(i)
               initCoor[1].append( cmd.get_model(str1i,1).get_coord_list()[0] )
               #initCoor[1].append(cmd.get_atom_coords(object_prefix+"1"+object_subgroup_prefix+str(i)))
           for i in range(wiz.subgroup_sum[2]):
               str2i = object_prefix+"2"+object_subgroup_prefix+str(i)
               initCoor[2].append( cmd.get_model(str2i,1).get_coord_list()[0] )           
               #initCoor[2].append(cmd.get_atom_coords(object_prefix+"2"+object_subgroup_prefix+str(i)))
           for i in range(wiz.subgroup_sum[3]):
               str3i = object_prefix+"3"+object_subgroup_prefix+str(i)
               initCoor[3].append( cmd.get_model(str3i,1).get_coord_list()[0] )              
               #initCoor[3].append(cmd.get_atom_coords(object_prefix+"3"+object_subgroup_prefix+str(i)))
           #print(initCoor) 




    def slider_move():

         form.angle_text.setText(str(form.slider_angle.value()))
    
         #cmd.alter_state(1,"_pw0_s0","(x,y,z)=(0.7755+"+str(form.slider_angle.value()/1800)+",0.0210+"+str(form.slider_angle.value()/1800)+",0.0+"+str(form.slider_angle.value()/1800)+")" )


         anchor = initCoor[0][0]
         point1 = initCoor[1][0]
         point2 = initCoor[2][0]
         point3 = initCoor[3][0]
         
         for j in range(3):
          for i in range(len(initCoor[j+1])):
            pp = initCoor[j+1][i]
            pp_new = Rotate2(anchor,pp, point1, point2, point3, form.slider_angle.value()*3.14159265/180.0)
            #print(pp_new)
            #t = "(x,y,z)=("+str(pp_new[0])+","+str(pp_new[1])+","+str(pp_new[2])+")"
            #print(t)
            cmd.alter_state(1,object_prefix+str(j+1)+object_subgroup_prefix+str(i),"(x,y,z)=("+str(pp_new[0])+","+str(pp_new[1])+","+str(pp_new[2])+")")
         cmd.rebuild() 
         

    def specify_angle():

        angle = int(form.angle_text.text())
        if angle < -180.0:
           angle = -180
           form.angle_text.setText("-180")
        if angle > 180.0:
           angle = 180
           form.angle_text.setText("180")  

        form.slider_angle.setValue(int(angle))    


    def revert_changes():
        form.angle_text.setText("0") 
        specify_angle()



    # hook up button callbacks
    #form.button_ray.clicked.connect(run)
    #form.button_browse.clicked.connect(browse_filename)
    form.button_close.clicked.connect(dialog.close)

    form.start.clicked.connect(start_wiz)
    form.pick_finish.clicked.connect(picking_finish)
    form.slider_angle.valueChanged.connect(slider_move)
    form.set_angle.clicked.connect(specify_angle)
    #form.angle_text.returnPressed.connect(specify_angle) # lead to undesired thing
    form.revert.clicked.connect(revert_changes)



    return dialog
