from PyQt5.QtWidgets import *
import sys
from PyQt5 import QtGui ,QtCore, QtWidgets
from PyQt5.QtCore import QRect
from functools import partial
import time
import os
from os.path import join
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Staubility_Code import *
# os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"


class PushButton:
    """
           A class used to implement Push button

           ...

           Attributes
           ----------
           button : QPushButton
               the actual push button from PyQt5 library

           Methods
           -------
           init(self)
               init the push button
           """
    def __init__(self,parent, name, location, Tooltip,function,arg1=None,arg2=None,arg3=None,arg4=None,arg5=None):
        """
        Parameters
        ----------
        parent : QGroupBox
            pointer to parent groupbox
        name : str
            The name of the button - will be displayed on the button
        location : 4dim location
            location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
        Tooltip : str
            description of the pushbutton functionality
        function : function
            link to the function that the pushbutton call to
        registers : list of register
            global registers list
        arg1,arg2,arg3,arg4,arg5:
            optiaonal input arguments to the callback function
        """
        self.button = QPushButton(parent=parent,text=name)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        if arg1==None:
            self.button.clicked.connect(partial(function))
        elif arg2==None:
            self.button.clicked.connect(partial(function, arg1))
        elif arg3==None:
            self.button.clicked.connect(partial(function, arg1,arg2))
        elif arg4 == None:
            self.button.clicked.connect(partial(function, arg1, arg2,arg3))
        elif arg5==None:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3,arg4))
        else:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3, arg4,arg5))
    def hide(self):
        self.button.hide()
    def show(self):
        self.button.show()

class FreeCheckButton:
    """
               A class used to implement Free Check button

               ...

               Attributes
               ----------
               button : QCheckButton
                   the actual check button from PyQt5 library

               Methods
               -------
               init(self)
                   init the check button
               """
    def __init__(self,parent, name, location, Tooltip,function=None,arg1=None,arg2=None,arg3=None,arg4=None,arg5=None):
        """
               Parameters
               ----------
               parent : QGroupBox
                   pointer to parent groupbox
               name : str
                   The name of the button - will be displayed near the button
               location : 4dim location
                   location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
               Tooltip : str
                   description of the checkbutton functionality
               function : function
                    link to the function that the pushbutton call to
               registers : list of register
                    global registers list
               arg1,arg2,arg3,arg4,arg5:
                    optiaonal input arguments to the callback function
               """
        self.button = QCheckBox(name, parent)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        self.button.setChecked(0)
        if function==None:
            return
        if arg1==None:
            self.button.clicked.connect(partial(function))
        elif arg2==None:
            self.button.clicked.connect(partial(function, arg1))
        elif arg3==None:
            self.button.clicked.connect(partial(function, arg1,arg2))
        elif arg4 == None:
            self.button.clicked.connect(partial(function, arg1, arg2,arg3))
        elif arg5==None:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3,arg4))
        else:
            self.button.clicked.connect(partial(function, arg1, arg2, arg3, arg4,arg5))
    def hide(self):
        self.button.hide()
    def show(self):
        self.button.show()
class FreeComboBox:
    """
      A class used to implement ComboBox -using for decimal or str values
      this comboBox will NOT change values in the registers and fields structure
      it will change other values
      ...

      Attributes
      ----------
      button : QComboBox
          the actual ComboBox from PyQt5 library

      Methods
      -------
      init(self)
          init the ComboBox
      """
    def __init__(self, parent, items, location, Tooltip,function,value=None):
        """
           Parameters
           ----------
           parent : QGroupBox
               pointer to parent groupbox
           name : str
               The name of the button - not in use
           items: list
                list of strings to display
           location : 4dim location
               location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
           Tooltip : str
               description of the ComboBox functionality
           address : link to variable
               name of the variable to get the chosen value
           function : link to function
               callback function of this button
           """
        self.button = QComboBox(parent=parent)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        self.button.addItems(items)
        if value!=None:
            self.button.setCurrentText(value)
        self.button.currentIndexChanged.connect(partial(function, self.button))
    def hide(self):
        self.button.hide()
    def show(self):
        self.button.show()
class FreeLineEdit:
    """
              A class used to implement LineEdit
              this LineEdit will NOT change values in the registers and fields structure
              it will change other values
              ...

              Attributes
              ----------
              button : QLineEdit
                  the actual LineEdit from PyQt5 library

              Methods
              -------
              init(self)
                  init the LineEdit
              """
    def __init__(self, parent, location, Tooltip,value=None,function =None):
        """
               Parameters
               ----------
               parent : QGroupBox
                   pointer to parent groupbox
               location : 4dim location
                   location of the button. suppose to get the next input:  QRect(10, 150, 120, 20)
               Tooltip : str
                   description of the ComboBox functionality
        """
        self.button = QLineEdit(parent=parent)
        self.button.setGeometry(location)
        self.button.setToolTip(Tooltip)
        if value==None:
            self.button.setText("16")
        else:
            self.button.setText(value)
        if function!= None:
            self.button.editingFinished.connect(partial(function, self.button))

    def hide(self):
        self.button.hide()
    def show(self):
        self.button.show()




class App(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "Staubility Optimizer Tool"
        self.left = 100
        self.top = 100
        self.width = 800
        self.height = 600
        self.icon = QtGui.QIcon("In\\STAUBILITY2.png")
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.setWindowIcon(self.icon)

        # # Add tabs to widget
        self.table_widget = InitApp(self)
        self.setCentralWidget(self.table_widget)

        self.show()


class InitApp(QWidget):
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)

        self.codon_methods = ["use_best_codon", "match_codon_usage", "harmonize_rca"]
        self.linker_types = ["fusion_linker", "2A"]
        self.codon_usage_opt = ["tAI", "NTE", "TDR", "fraction"]
        self.restrriction_opt = ["None"]
        for s in list(rest_dict.keys()):
            self.restrriction_opt.append(s + "_site")
        #
        self.input_target_gene = "None"
        self.output_path = "None"


        self.best_genes = [0]*10
        self.best_vals = [0]*10
        self.linkers = None
        self.number_of_optional_sites = [str(x) for x in range(16)]
        self.number_of_sites = self.number_of_optional_sites[10]
        self.miniGC = 0.3
        self.maxiGC = 0.7
        self.method = self.codon_methods[0]
        self.linker_type = self.linker_types[0]
        self.codon_usage = self.codon_usage_opt[0]
        self.restriction_enzyme = self.restrriction_opt[0]
        self.ess_gen = None
        self.final_linker = None

        self.input_disabled = False
        self.system_first_up = True
        self.finished_first_step = False
        self.finished_second_step = False
        self.finished_third_step = False


        self.main = parent
        self.layout = QVBoxLayout(self)

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tabs.resize(300, 200)

        # Initialize  tab screen
        self.gsa = QTabWidget()
        self.predictor = QTabWidget()

        # Add tabs
        self.tabs.addTab(self.gsa, "Staubility Optimizer Tool")


        # create tab
        self.init_tab()

        # set layout
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

    def browse_path(self):
        # load tx registes file
        options = QFileDialog.Options()

        fileName = QFileDialog.getExistingDirectory(self, "Select Folder For Output File", "Folder",
                                                  options=options)
        if fileName:
            self.output_path = fileName
            self.all_gsa.output_path_line.button.setText(fileName)

    def check_input_string(self):
        valid_aa = ['G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D','T']
        allowd_chars = ["A", "C", "G", "T"]
        ret = ""
        if not self.input_disabled:
            string = self.all_gsa.input_string_line.button.text().upper()
            if string == "NONE" or string == '':
                ret = "Wrong Input! Please Insert a String"
                QMessageBox.about(self, "User Warning",
                                 ret)
                return
            if string[:4]=="NONE":
                string = string[4:]
            if string[-3:] == "TAG" or string[-3:] == "TAA" or string[-3:] == "TGA":
                string = string[:-3]
                if len(string)<3:
                    ret+="Your sequence contains only stop codon\n"
            if len(string)%3!=0:
                ret+= "String must be divided by 3!\n"
            for char in string:
                if char not in allowd_chars:
                    if char==" ":
                        char = "'space'"
                    ret += "Character not allowed {} !\n".format(char)
            target_gene_aa = gene_seq_2_aa(string)
            if not (all(i in valid_aa for i in target_gene_aa)):
                ret += "There is a stop codon in the middle of your target gene.\nPlease check that the sequence you inserted is a valid reading frame\n"
            if ret != "":
                ret = "Wrong Input! Please fix the next errors and try again:\n\n"+ret
                QMessageBox.about(self, "User Warning",
                                 ret)
            else:
                self.input_target_gene = string
                self.all_gsa.input_string_line.button.setText(self.input_target_gene)
                self.input_disabled = True
                self.all_gsa.input_string_line.button.setDisabled(self.input_disabled)
                self.all_gsa.set_input.button.setText("Change")
        else:
            self.input_target_gene = "None"
            self.all_gsa.input_string_line.button.setText(self.input_target_gene)
            self.input_disabled = False
            self.all_gsa.input_string_line.button.setDisabled(self.input_disabled)
            self.all_gsa.set_input.button.setText("Insert")

    def change_linker_type(self,button):
        self.linker_type = button.currentText()

    def change_codon_usage(self,button):
        self.codon_usage = button.currentText()
        if self.codon_usage == "fraction":
            self.all_gsa.label_codon_usage_meth.show()
            self.all_gsa.codon_usage_meth_combo.show()
        else:
            self.all_gsa.label_codon_usage_meth.hide()
            self.all_gsa.codon_usage_meth_combo.hide()

    def change_codon_method(self, button):
        self.method = button.currentText()

    def change_gc_min(self,button):
        try:
            value = float(button.text())
            self.miniGC = value
        except:
            QMessageBox.about(self, "User Warning",
                              "GC Min must be a number")
            button.setText(str(self.miniGC))

    def change_gc_max(self,button):
        try:
            value = float(button.text())
            self.maxiGC = value
        except:
            QMessageBox.about(self, "User Warning",
                              "GC Max must be a number")
            button.setText(str(self.maxiGC))

    def change_restriction(self,button):
        self.restriction_enzyme = button.currentText()

    def change_efm_sites(self,button):
        self.number_of_sites = button.currentText()

    def open_user_guide(self):
        try:
            os.system("start Staubility_Enhancer_User_Guide.pdf")
        except:
            QMessageBox.about(self, "User Warning",
                              "Oh No!\nCan't Find The UserGuide file. Please locate it in the main folder and keep the original name (UserGuide.txt)")

    def write_to_log(self,text):
        log = self.all_gsa.log.button
        log.setText(text)
        time.sleep(0.5)
        QtCore.QCoreApplication.processEvents()

    def change_gene(self,button):
        self.ess_gen = button.currentText()

    def change_linker(self,button):
        final_linker = button.currentText()
        self.final_linker = self.linkers[final_linker]

    def plot_linker(self,linker):
        path = self.output_path + '\\Disorder_profile_{}.png'.format(linker)
        try:
            os.system("start "+path)
        except:
            QMessageBox.about(self, "User Warning",
                              "Oh No!\nCan't Find The UserGuide file. Please locate it in the main folder and keep the original name (UserGuide.txt)")

    def choose_genes(self):
        self.first_dialog = QtWidgets.QDialog(self)
        self.first_dialog.setWindowTitle("Predictor Results")
        self.first_dialog.setGeometry(100, 100, 500, 500)
        #
        self.first_dialog.label = QLabel("Step 1: Choose Conjugated Gene\n      from Best Candidates",self.first_dialog)
        self.first_dialog.label.setFont(QtGui.QFont('Arial', 15))
        self.first_dialog.label.setGeometry(QRect(110, 10, 300, 70))
        #
        self.first_dialog.label_name = QLabel("Conjugated Gene Name", self.first_dialog)
        self.first_dialog.label_name.setFont(QtGui.QFont('Arial', 12))
        self.first_dialog.label_name.setGeometry(QRect(40, 80, 200, 25))
        #
        self.first_dialog.label_score = QLabel("Conjugated Gene Normalized Score", self.first_dialog)
        self.first_dialog.label_score.setFont(QtGui.QFont('Arial', 12))
        self.first_dialog.label_score.setGeometry(QRect(230, 80, 250, 25))
        #
        for i in range(len(self.best_genes)):
            self.first_dialog.label_1 = QLabel("{}".format(self.best_genes[i]), self.first_dialog)
            self.first_dialog.label_1.setFont(QtGui.QFont('Arial', 10))
            self.first_dialog.label_1.setGeometry(QRect(70, 100+i*22, 100, 20))
            self.first_dialog.label_2 = QLabel("{0:3f}".format(self.best_vals[i]), self.first_dialog)
            self.first_dialog.label_2.setFont(QtGui.QFont('Arial', 10))
            self.first_dialog.label_2.setGeometry(QRect(280, 100+i*22, 140, 20))

        #
        self.first_dialog.choose_gene = QLabel("Choose Gene:", self.first_dialog)
        self.first_dialog.choose_gene.setFont(QtGui.QFont('Arial', 12))
        self.first_dialog.choose_gene.setGeometry(QRect(50, 350, 100, 20))
        self.first_dialog.gene_combo = FreeComboBox(self.first_dialog, self.best_genes, QRect(200, 350, 120, 20),
                                                      "Choose Conjugated Gene",self.change_gene)
        self.first_dialog.finish = PushButton(self.first_dialog, "Start Step 2\nLinker Candidates", QRect(200, 400, 130, 50),
                                        "Save your gene", self.start_step_2)
        self.first_dialog.finish.button.setFont(QtGui.QFont('Arial', 8))
        self.first_dialog.finish.button.setStyleSheet("background-color: green; font-weight: bold")


        self.first_dialog.exec_()

    def choose_linker(self):
        self.second_dialog = QtWidgets.QDialog(self)
        self.second_dialog.setWindowTitle("Best Linker Results")
        self.second_dialog.setGeometry(100, 100, 700, 500)
        #
        self.second_dialog.label = QLabel("Step 2: Choose Linker \n from Best Candidates", self.second_dialog)
        self.second_dialog.label.setFont(QtGui.QFont('Arial', 15))
        self.second_dialog.label.setGeometry(QRect(200, 10, 300, 70))
        #
        self.second_dialog.label_name = QLabel("Linker Name", self.second_dialog)
        self.second_dialog.label_name.setFont(QtGui.QFont('Arial', 12))
        self.second_dialog.label_name.setGeometry(QRect(20, 80, 100, 25))
        #
        self.second_dialog.label_seq = QLabel("Linker Sequence", self.second_dialog)
        self.second_dialog.label_seq.setFont(QtGui.QFont('Arial', 12))
        self.second_dialog.label_seq.setGeometry(QRect(120, 80, 250, 25))
        #
        self.second_dialog.label_score = QLabel("Linker Score", self.second_dialog)
        self.second_dialog.label_score.setFont(QtGui.QFont('Arial', 12))
        self.second_dialog.label_score.setGeometry(QRect(450, 80, 100, 25))
        #
        names = list(self.linkers.keys())
        vals = list(self.linkers.values())
        sequences = [val[0] for val in vals]
        scores = [val[1] for val in vals]
        indexes = list(np.array(scores).argsort())#[::-1])
        new_names = []
        loc = 0
        for i in indexes:
            self.second_dialog.label_1 = QLabel("{}".format(names[i]), self.second_dialog)
            self.second_dialog.label_1.setFont(QtGui.QFont('Arial', 10))
            self.second_dialog.label_1.setGeometry(QRect(20, 100 + loc * 22, 100, 20))
            self.second_dialog.label_2 = QLabel("{}".format(sequences[i]), self.second_dialog)
            self.second_dialog.label_2.setFont(QtGui.QFont('Arial', 10))
            self.second_dialog.label_2.setGeometry(QRect(120, 100 + loc * 22, 270, 20))
            self.second_dialog.label_3 = QLabel("{0:3f}".format(scores[i]), self.second_dialog)
            self.second_dialog.label_3.setFont(QtGui.QFont('Arial', 10))
            self.second_dialog.label_3.setGeometry(QRect(450, 100 + loc * 22, 200, 20))
            if self.linker_type =="fusion_linker":
                self.second_dialog.open_image =  PushButton(self.second_dialog, "Plot linker graph",
                                                  QRect(550, 100 + loc * 22, 120, 20),
                                                  "Plot linker {} graph - sequence:{}".format(names[i],sequences[i]), self.plot_linker,names[i])
            new_names.append(names[i])
            loc+=1


        #
        self.second_dialog.choose_gene = QLabel("Choose Linker:", self.second_dialog)
        self.second_dialog.choose_gene.setFont(QtGui.QFont('Arial', 12))
        self.second_dialog.choose_gene.setGeometry(QRect(50, 350, 120, 20))
        self.second_dialog.linker_combo = FreeComboBox(self.second_dialog, new_names, QRect(200, 350, 120, 20),
                                                    "Choose Linker", self.change_linker)
        self.second_dialog.finish = PushButton(self.second_dialog, "Start Step 3\nOptimization",
                                              QRect(200, 400, 130, 50),
                                              "Start optimization step", self.start_step_3)
        self.second_dialog.finish.button.setFont(QtGui.QFont('Arial', 8))
        self.second_dialog.finish.button.setStyleSheet("background-color: green; font-weight: bold")

        self.second_dialog.exec_()

    def start_step_1(self):
        if self.input_target_gene == "None":
            QMessageBox.about(self, "User Warning",
                              "Error: Invalid Input String!")
        elif self.output_path == "None":
            QMessageBox.about(self, "User Warning",
                              "Error: Invalid Output Path!")
        else:
            try:
                self.write_to_log("Starting Now - Please Wait Until The Process Will End!")
                ret,self.best_genes,self.best_vals = step_1(target_seq = self.input_target_gene, gene_data_list_path=join(getcwd(), r'In\utils\orf_coding_all.fasta.gz'),
                           genes_sequence_path=join(getcwd(), r'In\utils\seq_orf.csv'),
                           neighbors_indices_data_path=join(getcwd(), r'In\utils\max_5_ind_ORF.csv'),
                           chosen_path=join(getcwd(), r'In\utils\Final_Features_Table.csv'))
                if ret == "Success!":
                    log = "Finished the First Step!"
                    self.write_to_log(log)
                    print("\n\n"+log)
                    self.finished_first_step = True
                    self.choose_genes()
                else:
                    QMessageBox.about(self, "Message", ret)
                    log = "Failed! please try again the first step"
                    self.write_to_log(log)


            except Exception as e:
                QMessageBox.about(self, "ERROR",
                                  "Oh No!\nSomething Went Wrong, Please Try Again!\n\nThe Error is:\n{}".format(e))

    def start_step_2(self):
        self.change_gene(self.first_dialog.gene_combo.button)
        self.first_dialog.close()
        try:
            self.write_to_log("Starting Step 2 - Please Wait Until The Process Will End!")
            QMessageBox.about(self, "Message","Starting Step 2 - Finding Linker Candidates\n\nTo proceed into step 2 press OK")
            ret,self.linkers= step_2(target_seq = self.input_target_gene,genes_sequence_path= join(getcwd(), r'In\utils\seq_orf.csv'),
            essen_gene_orf = self.ess_gen,type_of_linker = self.linker_type,output_path = self.output_path,linkers = join(getcwd(), r'In\utils\linkers.csv'))
            if ret == "Success!":
                log = "Finished the Second Step!"
                self.write_to_log(log)
                print("\n\n" + log)
                self.finished_second_step = True
                self.choose_linker()
            else:
                QMessageBox.about(self, "Message", ret)
                log = "Failed! please try again the first step"
                self.write_to_log(log)
        except Exception as e:
            QMessageBox.about(self, "ERROR",
                              "Oh No!\nSomething Went Wrong, Please Try Again!\n\nThe Error is:\n{}".format(e))

    def start_step_3(self):
        self.change_linker(self.second_dialog.linker_combo.button)
        self.second_dialog.close()
        try:
            self.write_to_log("Starting Step 3 - Please Wait Until The Process Will End!")
            QMessageBox.about(self, "Message","Starting Step 3 - Optimization\nPlease wait, this process will take approximatly 10 minutes\n\nTo proceed into step 3 press OK")
            ret= step_3(target_seq = self.input_target_gene,genes_sequence_path= join(getcwd(), r'In\utils\seq_orf.csv'),essen_gene_orf= self.ess_gen,
                        linker_aa=self.final_linker[0],res_enzyme_name=self.restriction_enzyme,miniGC=self.miniGC,maxiGC=self.maxiGC,
                          opt_parameter=self.codon_usage,method=self.method,output_path =self.output_path,num_sites=int(self.number_of_sites))
            if ret == "Success!":
                log = "Finished the Last Step! - You Can Find Your Files at the Output Folder"
                self.write_to_log(log)
                print("\n\n" + log)
                self.finished_third_step_step = True
                QMessageBox.about(self, "Message", log)
            else:
                QMessageBox.about(self, "Message", ret)
                log = "Failed! please try again the first step"
                self.write_to_log(log)
        except Exception as e:
            QMessageBox.about(self, "ERROR",
                              "Oh No!\nSomething Went Wrong, Please Try Again!\n\nThe Error is:\n{}".format(e))



    def init_tab(self):
        self.all_gsa = QGroupBox("", self.gsa)
        self.all_gsa.setGeometry(QRect(0, 0, self.main.width, self.main.height))
        self.all_gsa.label = QLabel("Stability Optimizer",self.all_gsa)
        self.all_gsa.label.setFont(QtGui.QFont('Arial', 20))
        self.all_gsa.label.setGeometry(QRect(self.main.width/2-150, 20, 400, 35))
        #
        self.all_gsa.label1 = QLabel("Main Window",self.all_gsa)
        self.all_gsa.label1.setFont(QtGui.QFont('Arial', 15))
        self.all_gsa.label1.setGeometry(QRect(self.main.width/2-100, 45, 400, 35))
        #
        self.centralwidget = QtWidgets.QWidget(self.all_gsa)
        self.photo = QtWidgets.QLabel(self.centralwidget)
        self.photo.setGeometry(QtCore.QRect(600, 00, 170,170))
        self.photo.setPixmap(QtGui.QPixmap("In\\STAUBILITY2.png"))
        self.photo.setStyleSheet("border: 0")
        self.photo.setScaledContents(True)
        self.photo.setObjectName("photo")
        #

        self.all_gsa.label_input_string = QLabel("Target Gene String:", self.all_gsa)
        self.all_gsa.label_input_string.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_input_string.setGeometry(QRect(30, 100, 150, 20))
        self.all_gsa.input_string_line = FreeLineEdit(self.all_gsa, QRect(180, 100, 350, 20), "Insert String of Your Input Target Gene\nMust be divided by 3\nMust be Alphabet ACGT", self.input_target_gene)
        self.all_gsa.set_input = PushButton(self.all_gsa,"Insert",QRect(550, 100, 50, 20),"Check your string",self.check_input_string)
        #
        self.all_gsa.label_output_path = QLabel("Output Path:", self.all_gsa)
        self.all_gsa.label_output_path.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_output_path.setGeometry(QRect(30, 140,150, 20))
        self.all_gsa.output_path_line = FreeLineEdit(self.all_gsa, QRect(180, 140, 350, 20), "Insert Path To Desired Output Folder",self.output_path)
        self.all_gsa.browse_output = PushButton(self.all_gsa,"Browse",QRect(550, 140, 50, 20),"Browse Output Path",self.browse_path)
        #
        self.all_gsa.label_linker_types = QLabel("Linker Type:", self.all_gsa)
        self.all_gsa.label_linker_types.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_linker_types.setGeometry(QRect(30, 180, 200, 20))
        self.all_gsa.linker_type_combo = FreeComboBox(self.all_gsa, self.linker_types, QRect(180, 180, 100, 20), "Choose Linker Type", self.change_linker_type)
        #
        self.all_gsa.label2 = QLabel("Optimization Parameters",self.all_gsa)
        self.all_gsa.label2.setFont(QtGui.QFont('Arial', 15))
        self.all_gsa.label2.setGeometry(QRect(self.main.width/2-150, 230, 400, 35))
        #
        self.all_gsa.label_codon_usage = QLabel("Codon Usage Objective:", self.all_gsa)
        self.all_gsa.label_codon_usage.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_codon_usage.setGeometry(QRect(30, 270, 200, 20))
        self.all_gsa.codon_usage_combo = FreeComboBox(self.all_gsa, self.codon_usage_opt, QRect(300, 270, 120, 20), "Choose Codon Usage Objective", self.change_codon_usage)
        #
        self.all_gsa.label_codon_usage_meth = QLabel("Codon Usage Optimization Method:", self.all_gsa)
        self.all_gsa.label_codon_usage_meth.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_codon_usage_meth.setGeometry(QRect(30, 300, 270, 20))
        self.all_gsa.codon_usage_meth_combo = FreeComboBox(self.all_gsa, self.codon_methods, QRect(300, 300, 120, 20), "Choose Codon Usage Optimization Method", self.change_codon_method)
        self.all_gsa.label_codon_usage_meth.hide()
        self.all_gsa.codon_usage_meth_combo.hide()

        #
        self.all_gsa.label_gc_content = QLabel("GC Content:", self.all_gsa)
        self.all_gsa.label_gc_content.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_gc_content.setGeometry(QRect(30, 340, 150, 20))
        #
        self.all_gsa.label_gc_content_min = QLabel("Min:", self.all_gsa)
        self.all_gsa.label_gc_content_min.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_gc_content_min.setGeometry(QRect(180, 340, 100, 20))
        self.all_gsa.gc_min = FreeLineEdit(self.all_gsa, QRect(230, 340, 100, 20), "",
                                  str(self.miniGC), self.change_gc_min)
        #
        self.all_gsa.label_gc_content_max = QLabel("Max:", self.all_gsa)
        self.all_gsa.label_gc_content_max.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_gc_content_max.setGeometry(QRect(350, 340, 100, 20))
        self.all_gsa.gc_max = FreeLineEdit(self.all_gsa, QRect(400, 340, 100, 20), "",
                                  str(self.maxiGC), self.change_gc_max)
        #
        self.all_gsa.label_restriction_pattern = QLabel("Rest Enz Pattern To Avoid:", self.all_gsa)
        self.all_gsa.label_restriction_pattern.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_restriction_pattern.setGeometry(QRect(30, 370, 200, 20))
        self.all_gsa.restriction_combo = FreeComboBox(self.all_gsa, self.restrriction_opt, QRect(300, 370, 120, 20), "Choose Restriction Enzyme Pattern To Avoid", self.change_restriction)
        #
        self.all_gsa.label_efm_sites = QLabel("Number Of EFM Sites:", self.all_gsa)
        self.all_gsa.label_efm_sites.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.label_efm_sites.setGeometry(QRect(30, 400, 200, 20))
        self.all_gsa.efm_sites_combo = FreeComboBox(self.all_gsa, self.number_of_optional_sites, QRect(300, 400, 120, 20), "Choose Number Of EFM Sites", self.change_efm_sites,self.number_of_sites)
        #
        self.all_gsa.log_label = QLabel("Log:", self.all_gsa)
        self.all_gsa.log_label.setFont(QtGui.QFont('Arial', 12))
        self.all_gsa.log_label.setGeometry(QRect(400, 170, 200, 20))
        self.all_gsa.log = FreeLineEdit(self.all_gsa,QRect(400, 190, 350, 20), "Here you can get some updates","Waiting to start!")
        #
        self.all_gsa.user_guide = PushButton(self.all_gsa,"Open User Guide",QRect(250, 450, 100, 50),"Check Out Our User Guide",self.open_user_guide)
        self.all_gsa.user_guide.button.setFont(QtGui.QFont('Arial', 8))
        self.all_gsa.user_guide.button.setStyleSheet("background-color: turquoise; font-weight: bold")
        #
        self.all_gsa.step1 = PushButton(self.all_gsa,"Start Step 1",QRect(400, 450, 100, 50),"Calculate the features and find the best  gene",self.start_step_1)
        self.all_gsa.step1.button.setFont(QtGui.QFont('Arial', 8))
        self.all_gsa.step1.button.setStyleSheet("background-color: green; font-weight: bold")
        #
        # self.all_gsa.step2 = PushButton(self.all_gsa,"Start Step 2",QRect(400, 450, 100, 50),"Start from Step 2",self.start_step_2)
        # self.all_gsa.step2.button.setFont(QtGui.QFont('Arial', 8))
        # self.all_gsa.step2.button.setStyleSheet("background-color: green; font-weight: bold")
        # self.all_gsa.step2.hide()
        # #
        # self.all_gsa.step3 = PushButton(self.all_gsa,"Start Step 3",QRect(550, 450, 100, 50),"Start from step 3",self.start_step_3)
        # self.all_gsa.step3.button.setFont(QtGui.QFont('Arial', 8))
        # self.all_gsa.step3.button.setStyleSheet("background-color: green; font-weight: bold")
        # self.all_gsa.step3.hide()


if __name__=="__main__":
    print("##########################################################################################\n")
    print("PLEASE WAIT WHILE THE SOFTWARE TURNS ON- IT WILL TAKE A FEW SECONDS\n")
    print("##########################################################################################")
    # Run Program
    app = QApplication(sys.argv)
    # if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    #     app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    #
    # if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    #     app.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
    ex = App()
    sys.exit(app.exec())