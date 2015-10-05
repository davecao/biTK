# -*- coding: utf-8 -*-
"""
    This module depends on Package Jinja2 for export HTML file
"""
import os
import sys
import glob
import shutil

from time import localtime, strftime
from biTK import TPLPATH, jinja2_ENV, VERSION
from biTK.ngs.io.pathtools import isWritable, makePath
#from jinja2 import Environment, FileSystemLoader

__all__ = ['AnalysisItems','StatsAnalysisReport']

class HTMLview(object):
    """
        variables of HTML template
    """
    def __init__(self):
        """
            Initialization
        """
        pass

class AnalysisItems(object):
    def __init__(self, name=None, description=None, data=None, 
                      hasTable=False, tableColName=[], hasImg=False,
                      img_path=None, img_name=None, 
                      tplName="section.html"):

        super(AnalysisItems, self).__init__()
        self.name = name if name else ''
        self.description = description if description else ''
        self.data = data
        self.hasTable = hasTable
        self.tableColName = tableColName
        self.hasImg = hasImg
        self.img_path = img_path if img_path else ''
        self.img_name = img_name if img_name else ''
        #self.template = tplName
        #print(jinja2_ENV.list_templates())

        self.template = jinja2_ENV.get_template(tplName)

    def __get_elements__(self):
        elems = {
            'name':self.name, 
            'description': self.description,
            'content_url': ".",
            'thumnail_url':"." + os.sep + self.img_name,
            'hasImg': self.hasImg,
            'hasTable': self.hasTable,
            'table' : self.data,
            'id': self.name.replace(' ', ''),
            'tableColName': self.tableColName
        }
        return elems

    def render(self, format='html'):
        return self.template.render(self.__get_elements__())

class StatsAnalysisReport(object):
    """docstring for StatsAnalysisReport"""
    def __init__(self, outputdir, outputfile, 
                infoBox = [],
                backupCount = 5, 
                title = None, 
                tplName = "report.html"):
        """
            Initialization
        Args:
            outputpath (str): the path of output file
            title (str) : the title of the report
            tplName (str) : file name of the jinja2 template

        Return:
             StatsAnalysisReport object
        """
        super(StatsAnalysisReport, self).__init__()
        self.title = title
        self.odir = outputdir
        self.ofile = outputfile
        self.template = jinja2_ENV.get_template(tplName)
        self.analysisItems = []
        self.backupCount = backupCount
        self.infoBox = infoBox
        self.doRollover()

    def doRollover(self):
        """ Determine if rollover should occur """
#        directories = self.odir.split(os.sep)
#        search_dir = os.path.join(directories[:-1])[0]
#        dir_name = directories[-1]
#        n = (glob.glob("{}{}*{}*".format(search_dir,os.sep,dir_name))[-1]).split('.')[-1]
#        n = int(n)
        if self.backupCount > 0:
            for i in range(self.backupCount - 1, 0, -1):
                sfn = "%s.%d" % (self.odir, i)
                dfn = "%s.%d" % (self.odir, i + 1)
                if os.path.exists(sfn):
                    #print "%s -> %s" % (sfn, dfn)
                    if os.path.exists(dfn):
                        try:
                            os.remove(dfn)
                        except:
                            shutil.rmtree(dfn, ignore_errors=True)
                    os.rename(sfn, dfn)
            dfn = self.odir + ".1"
            if os.path.exists(dfn):
                try:
                    os.remove(dfn)
                except Exception:
                    shutil.rmtree(dfn, ignore_errors=True)
            # Issue 18940: A file may not have been created if delay is True.
            if os.path.exists(self.odir):
                os.rename(self.odir, dfn)
        self.outputpath = makePath(self.odir)

    def append(self, analysisItems):
        """
            Append results 
        """
        if isinstance(analysisItems, AnalysisItems):
            self.analysisItems.append(analysisItems)
        else:
            print("TypeError: input argument is not AnalysisItems")
            sys.exit(1)

    def render(self, format='html'):
        block = ''
        li_item = ''
        for item in self.analysisItems:
            block +=item.render()
            li_item += "<li><a href='#"+item.name.replace(' ','') +"'>"+item.name+"</a></li>"
            #sd_content.append(li_item)
        #result = self.template.render({'title':self.title, 'items':block})
        with open(self.outputpath+os.sep+self.ofile, 'w') as fout:
            fout.write(self.template.render({
                'package' : "biTK ver "+ VERSION,
                'title' : self.title, 
                'rptTime' : strftime("%d %b %Y %H:%M:%S", localtime()),
                'body_content' : block,
                'sidebar_content' : li_item,
                'infoBox' : self.infoBox
                }))
