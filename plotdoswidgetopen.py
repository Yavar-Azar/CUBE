# importing various libraries 
import sys 
from PyQt5.QtWidgets import QWidget,QDialog, QApplication, QPushButton, QVBoxLayout, QComboBox , QHBoxLayout,QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas 
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar 
import matplotlib.pyplot as plt 
import random 
import numpy as np

# main window 
# which inherits QDialog 
class plotwin(QWidget): 
	
	# constructor 
    def __init__(self, parent=None): 
        super(plotwin, self).__init__(parent) 

		# a figure instance to plot on 
        self.figure = plt.figure() 

		# this is the Canvas Widget that 
		# displays the 'figure'it takes the 
		# 'figure' instance as a parameter to __init__ 
        self.canvas = FigureCanvas(self.figure) 

		# this is the Navigation widget 
		# it takes the Canvas widget and a parent 
        self.toolbar = NavigationToolbar(self.canvas, self) 
        

        self.buttondata=QPushButton('Open data')
        
		# Just some button connected to 'plot' method 
        self.button = QPushButton('Plot') 
		
		# adding action to the button 
        self.button.clicked.connect(self.plot)
        self.buttondata.clicked.connect(self.open_dialog_box)
        self.combox = QComboBox()
        self.comboy = QComboBox()
        self.datapath = "None"
        self.headlines = []


   		# creating a Vertical Box layout 
        layout = QVBoxLayout() 
		# adding tool bar to the layout 
        layout.addWidget(self.toolbar) 
		# adding canvas to the layout 
        layout.addWidget(self.canvas)

		# adding push button to the layout 
        horizlayout=QHBoxLayout()
        horizlayout.addWidget(self.buttondata)
        horizlayout.addWidget(self.combox)
        horizlayout.addWidget(self.comboy)

        layout.addLayout(horizlayout)
        layout.addWidget(self.button)
        # setting layout to the main window 
        self.setLayout(layout)


    # action called by thte push button



    def open_dialog_box(self):
        filename = QFileDialog.getOpenFileName()
        path = filename[0]
        directad = path.split('/')[:-1]
        fname=path.split('/')[-1]
        
        basename=fname.split(".")[0]
        print(path)
        print("your data file is ", fname) 
        
        f=open(path, "r")
        firstline=f.readline()
        headers=firstline.split()[1:]
        f.close()
        
        
        self.datapath = path
        list1=headers
        self.combox.addItems(list1)
        self.comboy.addItems(list1)
        for val in list1:
        	self.headlines.append(val)
        print(self.headlines)




    


    def plot(self):

        xstr = self.combox.currentText()
        ystr = self.comboy.currentText()

        xind = self.headlines.index(xstr)
        yind = self.headlines.index(ystr)
        data = np.loadtxt(self.datapath)
        # clearing old figure
        self.figure.clear()
# create an axis
        ax = self.figure.add_subplot(111) 

        # plot data
        ax.set_xlim(min(data[:,xind]), 1.1*max(data[:,xind]))
        ax.set_ylim(min(data[:,yind]), 1.1*max(data[:,yind]))
        ax.set_xlabel(xstr, fontsize=14)
        ax.set_ylabel(ystr, fontsize=12) 
        ax.plot(data[:,xind], data[:,yind], '-') 
        plt.tight_layout()
        plt.savefig("test00.png",bbox_inches='tight', dpi=300)


		# refresh canvas 
        self.canvas.draw() 

# driver code 
if __name__ == '__main__': 
	
	# creating apyqt5 application 
	app = QApplication(sys.argv) 

	# creating a window object 
	main = plotwin() 
	
	# showing the window 
	main.show() 

	# loop 
	sys.exit(app.exec_()) 

