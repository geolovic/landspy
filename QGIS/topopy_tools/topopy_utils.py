# -*- coding: utf-8 -*-

"""
***************************************************************************
    topopy_utils.py
    ---------------------
    Date                 : December 2020
    Copyright            : (C) 2020 by J. Vicente Perez
    Email                : vperez at go dot ugr dot es

***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'J. Vicente Perez'
__date__ = 'December 2020'
__copyright__ = '(C) 2020, J. Vicente Perez'

import os
import inspect
import subprocess

TOPOPY_FOLDER = os.path.split(inspect.getfile(inspect.currentframe()))[0]

class TopopyUtils:
    
    @staticmethod
    def getFolder():
        return TOPOPY_FOLDER

    @staticmethod
    def runApp(commands, feedback):
        commandline = " ".join(commands)
        #commandline = "python " + cmd_folder + "/apps/testApp2.py"
        feedback.pushConsoleInfo("Execute topopy App")
        feedback.pushConsoleInfo(commandline)
        #exec(open('C:/Users/Usuario/AppData/Roaming/QGIS/QGIS3/profiles/default/python/plugins/topopy_tools/apps/testApp.py'.encode('utf-8')).read())
        output = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stdin=open(os.devnull), stderr=subprocess.STDOUT, universal_newlines=False).communicate()[0]
        feedback.pushConsoleInfo(output.decode("utf-8"))
