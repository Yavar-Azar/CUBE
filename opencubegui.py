import sys

from PyQt5.QtGui import QPixmap, QIcon, QCursor, QIntValidator

from PyQt5.QtCore import (QCoreApplication, 
                          QObject,
                          QRunnable,
                          QThread,
                          QThreadPool,
                          pyqtSignal)

from PyQt5.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
    QSpacerItem,
    QSizePolicy,   
    QLineEdit,
    QFileDialog,
    QComboBox,
    QProgressBar,
    QMessageBox,
    QGridLayout,
    QGroupBox,
    QTextEdit,
    QTableWidget,
    QTableWidgetItem,
    QSpacerItem
)


from CUBE_analyzer import cube1

class cubegui(QWidget):

    def __init__(self):
        super().__init__()
        self.setFixedSize(600, 600)
        self.setWindowTitle("Cube analyze")
        self.cubeobject = None


        opencubepb = QPushButton("open cube file")

        self.label1 = QLabel()
        self.cubenamelb = QLabel()
        self.cubetextte = QTextEdit()
        self.cubetextte.setReadOnly(True)        

        opencubegroup = QGroupBox("CUBE FILE")
        analyzecube = QGroupBox("Analyze")


        grid1 = QGridLayout()
        grid1.addWidget(opencubepb)
        grid1.addWidget(self.label1, 0, 1)
        grid1.addWidget(self.cubenamelb, 0, 2)
        grid1.addWidget(self.cubetextte, 1, 0, 1, 3)

        opencubegroup.setLayout(grid1)


        selectactionlb = QLabel("Select act")
        
        self.actscb = QComboBox()
        self.actscb.addItems(["sum", "", "c"])
        self.plotxypb = QPushButton("plotxy")
        self.plotxypb.setEnabled(False)
        

        gridanalyze = QGridLayout()
        gridanalyze.addWidget(selectactionlb ,0, 0)
        gridanalyze.addWidget(self.actscb, 0, 1)
        gridanalyze.addWidget(self.plotxypb)
        
        
        
        
        analyzecube.setLayout(gridanalyze)
        
        mylayout = QVBoxLayout()
        

        mylayout.addWidget(opencubegroup)
        mylayout.addWidget(analyzecube)
        self.setLayout(mylayout)




###########################  Connections


        opencubepb.clicked.connect(self.opencube_dialog_box)




    def opencube_dialog_box(self):
        filename = QFileDialog.getOpenFileName(None, "Select a file...", "./", filter="cube files (*.cube)")
        fpath = filename[0]
        print("this is path :   "+fpath)
        directoryl = fpath.split('/')[:-1]
#       directory = self.currentdir
        fname=fpath.split('/')[-1]
        basename=fname.split(".")[0]
        print("your cif file is ", fname) 
 #       with open(path, "r") as f:
 #           print(f.readline())
 #       qeinpgen(fname, basename)
        self.base1=basename
        self.cifname=fname
        self.label1.setText('Selected cube file is')
        self.cubenamelb.setText(fname)
        self.cifpath=fpath
        self.cubeobject = cube1(fpath)

        self.extract_data(fpath)
        self.plotxypb.setEnabled(True)
        self.plotxypb.clicked.connect(self.plotxyslice)





# !!! there is a problem in nwchem cubes



    def extract_data(self, cubefile):
        datastr = "Data extracted from {} file \n\n".format(self.cifname)
        temp = self.cubeobject
        datastr = datastr + "Number of atoms is \n {} \n".format(str(temp.natoms))
        datastr = datastr + "sum of data * dv (total charge for charge data) \n {:8.4f} \n\n".format(float(temp.datasum))
        datastr = datastr + "Cell parameters are \n{} \n\n\n".format(str(temp.cell))
        datastr = datastr + "atomic positions in XYZ format is \n {} \n".format(str(temp.position))
        self.cubetextte.setText(datastr)
        
        
    def plotxyslice(self):
        temp = self.cubeobject 
        print("clicked")
        

        








app = QApplication(sys.argv)
w = cubegui()
w.show()
app.exec()