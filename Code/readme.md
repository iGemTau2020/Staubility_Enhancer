Staubility Enhancer Code
-----------------------
This folder contains the source code of our product, the Staubility Enhancer.

Our code combines many frameworks, some of whom only work on certain operating systems, and some of whom require highly nontrivial installations and dependencies. Thus, running our code as a script is extremely tricky. This is part of the reason we built the GUI framework, to avoid having our users go through this process, and compiled it into a software.exe, so our users can run it without creating and installing the different dependencies.

Our software.exe was compiled from the code in this folder which is specified in the next section.

Folder Contents
----------------
 - Staubility_code.py - The staubility code contains the different models of the project (feature_generation,feature_selection, linker module,sequence_optimizer,conjugated_gene_predictor) inside one file with few alignments for running inside a complex program, and a main function which activates the models according to the user inputs and the program step.
 
 - Gui_Code.py -  The Gui code is our software interface code. The user interface was created with the unique framework we built for this project, based on the pyqt library. The Gui_code file contains the customer side functions which test and modify the different inputs, the user interfaces structured code and our framework.
