import sys
import json
import os
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QComboBox,
    QCheckBox,
    QDoubleSpinBox,
    QSpinBox,
    QFileDialog,
    QScrollArea,
    QGridLayout,
    QFrame,
)
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import matplotlib.pyplot as plt
from datetime import datetime


class SimulationGUI(QWidget):

    def __init__(self):
        super().__init__()
        self.params_file = "guiconfig/last_params.json"
        self.counter_file = "guiconfig/run_counter.json"
        self.icon_file = "guiconfig/question_icon.png"
        self.resize(1400, 800)  # Set the initial size of the GUI
        self.setMinimumSize(500, 500)  # Set the minimum size of the GUI
        self.initUI()

    def initUI(self):
        self.main_layout = QVBoxLayout(self)

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        self.main_layout.addWidget(scroll_area)

        content_widget = QWidget()
        self.content_layout = QHBoxLayout(content_widget)

        scroll_area.setWidget(content_widget)

        # Create the columns using QVBoxLayout
        self.column1 = QVBoxLayout()
        self.column2 = QVBoxLayout()
        self.column3 = QVBoxLayout()

        # Create vertical lines (separators) between the columns
        self.line1 = QFrame()
        self.line1.setFrameShape(QFrame.VLine)
        self.line1.setFrameShadow(QFrame.Sunken)
        
        self.line2 = QFrame()
        self.line2.setFrameShape(QFrame.VLine)
        self.line2.setFrameShadow(QFrame.Sunken)

        # Add the columns and lines to the horizontal layout
        self.content_layout.addLayout(self.column1)
        self.content_layout.addWidget(self.line1)
        self.content_layout.addLayout(self.column2)
        self.content_layout.addWidget(self.line2)
        self.content_layout.addLayout(self.column3)

        # Add your fields to the first column
        self.particleTypeInput = self.addFieldWithExplanation(
            "Particle Type:",
            QComboBox(self),
            explanation="Select the type of particle to simulate.",
            column=self.column1
        )
        self.particleTypeInput.addItems(["Pion", "Eta"])

        self.nParticlesInput = self.addFieldWithExplanation(
            "Number of Particles:",
            QSpinBox(self),
            max_value=100000000,
            default_value=8000000,
            explanation="Enter the number of particles to simulate.",
            column=self.column1
        )

        self.ptMaxInput = self.addFieldWithExplanation(
            "PT Max (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=50,
            explanation="Enter the maximum transverse momentum (PT) in GeV.",
            column=self.column1
        )

        self.ptMinInput = self.addFieldWithExplanation(
            "PT Min (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=0,
            explanation="Enter the minimum transverse momentum (PT) in GeV.",
            column=self.column1
        )

        self.weightMethodInput = self.addFieldWithExplanation(
            "Weight Method:",
            QComboBox(self),
            explanation="Select the weighting method.",
            column=self.column1
        )
        self.weightMethodInput.addItems(["EXP", "POWER", "WSHP", "HAGEDORN"])

        self.asymmCutInput = self.addFieldWithExplanation(
            "Apply Asymm Cut:",
            QCheckBox(self),
            default_value=True,
            explanation="Check to apply asymmetric cut.",
            column=self.column1
        )

        self.asymmCutValueInput = self.addFieldWithExplanation(
            "Asymm Cut Value:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=3,
            default_value=0.6,
            explanation="Enter the value for the asymmetric cut.",
            column=self.column1
        )

        self.clusterOverlapInput = self.addFieldWithExplanation(
            "Cluster Overlap:",
            QCheckBox(self),
            default_value=True,
            explanation="Check to enable cluster overlap.",
            column=self.column1
        )

        self.clusterOverlapProbInput = self.addFieldWithExplanation(
            "Cluster Overlap Probability:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=3,
            step=0.001,
            default_value=0.99,
            explanation="Enter the probability of NO cluster overlap. i.e. 0.99 means 1% overlap probability.",
            column=self.column1
        )

        # Add fields to the second column
        self.deltaRcutMaxInput = self.addFieldWithExplanation(
            "Delta R Cut Max:",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=1.1,
            explanation="Enter the maximum delta R cut value.",
            column=self.column2
        )

        self.pt1cutInput = self.addFieldWithExplanation(
            "PT1 Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=1.5,
            explanation="Enter the cut value for PT1 in GeV.",
            column=self.column2
        )

        self.pt2cutInput = self.addFieldWithExplanation(
            "PT2 Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=1.5,
            explanation="Enter the cut value for PT2 in GeV.",
            column=self.column2
        )

        self.combPtCutInput = self.addFieldWithExplanation(
            "Combined PT Cut Multiplier:",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=0,
            explanation="This value is multiplied by the sum of the photon cuts. If the pion has a pt smaller than that value, it is cut.",
            column=self.column2
        )

        self.ptMaxCutInput = self.addFieldWithExplanation(
            "PT Max Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=50,
            explanation="Enter the maximum PT cut value in GeV.",
            column=self.column2
        )

        self.PT_Max_binInput = self.addFieldWithExplanation(
            "PT Max Bin:",
            QSpinBox(self),
            max_value=100,
            default_value=20,
            explanation="Specify the maximum PT bin.",
            column=self.column2
        )

        self.MassNBinsInput = self.addFieldWithExplanation(
            "Mass N Bins:",
            QSpinBox(self),
            max_value=5000,
            default_value=1200,
            explanation="Specify the number of mass bins.",
            column=self.column2
        )

        self.binresInput = self.addFieldWithExplanation(
            "Bin Resolution:",
            QSpinBox(self),
            max_value=100,
            default_value=2,
            explanation="Specify the pT bin resolution.",
            column=self.column2
        )

        self.n_binsInput = self.addFieldWithExplanation(
            "Number of Bins:",
            QSpinBox(self),
            max_value=1000,
            default_value=80,
            explanation="Specify the number of pT bins.",
            column=self.column2
        )

        # Add fields to the third column
        self.DebugInput = self.addFieldWithExplanation(
            "Debug Mode:",
            QCheckBox(self),
            default_value=False,
            explanation="Enable or disable debug mode.",
            column=self.column3
        )

        self.etCutInput = self.addFieldWithExplanation(
            "ET Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=1.0,
            explanation="Enter the ET cut value in GeV.",
            column=self.column3
        )

        self.Apply_Eta_CutInput = self.addFieldWithExplanation(
            "Apply Eta Cut:",
            QCheckBox(self),
            default_value=False,
            explanation="Check to apply the eta cut.",
            column=self.column3
        )

        self.eta_cut_valInput = self.addFieldWithExplanation(
            "Eta Cut Value:",
            QDoubleSpinBox(self),
            max_value=10,
            decimals=2,
            default_value=0.6,
            explanation="Enter the value for the eta cut.",
            column=self.column3
        )

        self.smeared_lower_bin_limitInput = self.addFieldWithExplanation(
            "Smeared Lower Bin Limit:",
            QDoubleSpinBox(self),
            max_value=10,
            decimals=2,
            default_value=0.0,
            explanation="Specify the lower limit for (smeared) mass bins.",
            column=self.column3
        )

        self.smeared_upper_bin_limitInput = self.addFieldWithExplanation(
            "Smeared Upper Bin Limit:",
            QDoubleSpinBox(self),
            max_value=10,
            decimals=2,
            default_value=1.2,
            explanation="Specify the upper limit for (smeared) mass bins.",
            column=self.column3
        )

        self.smear_factor_sqrtEInput = self.addFieldWithExplanation(
            "Smear Factor sqrt(E):",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=4,
            default_value=0.154,
            explanation="Specify the smearing factor A/sqrt(E).",
            column=self.column3
        )

        self.nclusPtCutInput = self.addFieldWithExplanation(
            "NCLUS PT Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=0.0,
            explanation="Enter the NCLUS PT cut value in GeV.",
            column=self.column3
        )

        self.positSmearingFactorInput = self.addFieldWithExplanation(
            "Position Smearing Factor (mm):",
            QDoubleSpinBox(self),
            max_value=10,
            decimals=2,
            default_value=2.8,
            explanation="Enter the position smearing factor.",
            column=self.column3
        )

        self.saveToTreeInput = self.addFieldWithExplanation(
            "Save to Tree:",
            QCheckBox(self),
            default_value=False,
            explanation="Check to save the results to a tree.",
            column=self.column3
        )

        self.SmearFactorconstInput = self.addFieldWithExplanation(
            "Base Smearing Percent:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=3,
            step=0.001,
            default_value=0.168,
            explanation="Enter the base constant smearing percent.",
            column=self.column3
        )

        self.nStepsInput = self.addFieldWithExplanation(
            "Number of Steps:",
            QSpinBox(self),
            max_value=100,
            default_value=1,
            explanation="Enter the number of constant smearing steps for the simulation to run over. Currently the only itterative value is the constant smearing factor.",
            column=self.column3
        )

        self.stepSizeInput = self.addFieldWithExplanation(
            "Step Size:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=4,
            step=0.0001,
            default_value=0.001,
            explanation="Enter the step size for the simulation.",
            column=self.column3
        )

        self.outputFormatInput = self.addFieldWithExplanation(
            "Output Format:",
            QComboBox(self),
            explanation="Select the format for saving the output.",
            column=self.column3
        )
        self.outputFormatInput.addItems(["PDF", "PNG", "JPEG"])

        self.runButton = QPushButton("Run Simulation", self)
        self.main_layout.addWidget(self.runButton)
        # Connect the button click to the runSimulation method
        self.runButton.clicked.connect(self.runSimulation)

        self.setStyleSheet(
            """
            QLabel {
                font-size: 18px;
            }
            QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox {
                font-size: 18px;
            }
            QPushButton {
                font-size: 20px;
            }
            QScrollArea {
                background-color: white;
            }
        """
        )

        self.loadParams()

    def addFieldWithExplanation(
        self,
        label_text,
        widget,
        max_value=None,
        decimals=None,
        step=None,
        default_value=None,
        explanation=None,
        column=None,  # New parameter to specify the column
    ):
        label = QLabel(label_text, self)
        explanation_label = QLabel(self)
        explanation_pixmap = QPixmap(self.icon_file)
        explanation_pixmap = explanation_pixmap.scaled(
            20, 20, Qt.KeepAspectRatio, Qt.SmoothTransformation
        )
        explanation_label.setPixmap(explanation_pixmap)
        explanation_label.setToolTip(explanation)
        explanation_label.setCursor(Qt.PointingHandCursor)

        layout = QHBoxLayout()
        layout.addWidget(label)
        layout.addWidget(widget)
        layout.addWidget(explanation_label)

        # Set fixed width for input fields
        widget.setFixedWidth(150)
        explanation_label.setFixedWidth(20)
        explanation_label.setFixedHeight(20)

        container = QWidget()
        container.setLayout(layout)
        column.addWidget(container)  # Add the container to the specified column

        if isinstance(widget, QDoubleSpinBox):
            if max_value is not None:
                widget.setMaximum(max_value)
            if decimals is not None:
                widget.setDecimals(decimals)
            if step is not None:
                widget.setSingleStep(step)
            if default_value is not None:
                widget.setValue(default_value)

        if isinstance(widget, QSpinBox):
            if max_value is not None:
                widget.setMaximum(max_value)
            if default_value is not None:
                widget.setValue(default_value)

        if isinstance(widget, QCheckBox):
            if default_value is not None:
                widget.setChecked(default_value)

        return widget
    
    def loadParams(self):
        if not os.path.exists("guiconfig"):
            os.makedirs("guiconfig")

        if os.path.exists(self.params_file):
            with open(self.params_file, "r") as file:
                params = json.load(file)
                self.particleTypeInput.setCurrentText(params.get("particleType", "Pion"))
                self.nParticlesInput.setValue(params.get("nParticles", 8000000))
                self.ptMaxInput.setValue(params.get("ptMax", 50))
                self.ptMinInput.setValue(params.get("ptMin", 0))
                self.weightMethodInput.setCurrentText(params.get("weightMethod", "WSHP"))
                self.asymmCutInput.setChecked(params.get("applyAsymmCut", True))
                self.asymmCutValueInput.setValue(params.get("asymmCutValue", 0.6))
                self.clusterOverlapInput.setChecked(params.get("clusterOverlap", True))
                self.clusterOverlapProbInput.setValue(params.get("clusterOverlapProb", 0.99))
                self.deltaRcutMaxInput.setValue(params.get("deltaRcutMax", 1.1))
                self.pt1cutInput.setValue(params.get("pt1cut", 1.5))
                self.pt2cutInput.setValue(params.get("pt2cut", 1.5))
                self.combPtCutInput.setValue(params.get("combPtCut", 0))
                self.ptMaxCutInput.setValue(params.get("ptMaxCut", 50))
                self.PT_Max_binInput.setValue(params.get("PT_Max_bin", 20))
                self.MassNBinsInput.setValue(params.get("MassNBins", 1200))
                self.binresInput.setValue(params.get("binres", 2))
                self.n_binsInput.setValue(params.get("n_bins", 80))
                self.DebugInput.setChecked(params.get("Debug", False))
                self.etCutInput.setValue(params.get("etCut", 1.0))
                self.Apply_Eta_CutInput.setChecked(params.get("Apply_Eta_Cut", False))
                self.eta_cut_valInput.setValue(params.get("eta_cut_val", 0.6))
                self.smeared_lower_bin_limitInput.setValue(params.get("smeared_lower_bin_limit", 0.0))
                self.smeared_upper_bin_limitInput.setValue(params.get("smeared_upper_bin_limit", 1.2))
                self.smear_factor_sqrtEInput.setValue(params.get("smear_factor_sqrtE", 0.154))
                self.nclusPtCutInput.setValue(params.get("nclusPtCut", 0.0))
                self.positSmearingFactorInput.setValue(params.get("positSmearingFactor", 2.8))
                self.saveToTreeInput.setChecked(params.get("saveToTree", False))
                self.SmearFactorconstInput.setValue(params.get("SmearFactorconst", 0.168))
                self.nStepsInput.setValue(params.get("nSteps", 1))
                self.stepSizeInput.setValue(params.get("stepSize", 0.001))
                self.outputFormatInput.setCurrentText(params.get("outputFormat", "PDF"))
            
    def saveParams(self):
        if not os.path.exists("guiconfig"):
            os.makedirs("guiconfig")

        params = {
            "particleType": self.particleTypeInput.currentText(),
            "nParticles": self.nParticlesInput.value(),
            "ptMax": self.ptMaxInput.value(),
            "ptMin": self.ptMinInput.value(),
            "weightMethod": self.weightMethodInput.currentText(),
            "applyAsymmCut": self.asymmCutInput.isChecked(),
            "asymmCutValue": self.asymmCutValueInput.value(),
            "clusterOverlap": self.clusterOverlapInput.isChecked(),
            "clusterOverlapProb": self.clusterOverlapProbInput.value(),
            "deltaRcutMax": self.deltaRcutMaxInput.value(),
            "pt1cut": self.pt1cutInput.value(),
            "pt2cut": self.pt2cutInput.value(),
            "combPtCut": self.combPtCutInput.value(),
            "ptMaxCut": self.ptMaxCutInput.value(),
            "PT_Max_bin": self.PT_Max_binInput.value(),
            "MassNBins": self.MassNBinsInput.value(),
            "binres": self.binresInput.value(),
            "n_bins": self.n_binsInput.value(),
            "Debug": self.DebugInput.isChecked(),
            "etCut": self.etCutInput.value(),
            "Apply_Eta_Cut": self.Apply_Eta_CutInput.isChecked(),
            "eta_cut_val": self.eta_cut_valInput.value(),
            "smeared_lower_bin_limit": self.smeared_lower_bin_limitInput.value(),
            "smeared_upper_bin_limit": self.smeared_upper_bin_limitInput.value(),
            "smear_factor_sqrtE": self.smear_factor_sqrtEInput.value(),
            "nclusPtCut": self.nclusPtCutInput.value(),
            "positSmearingFactor": self.positSmearingFactorInput.value(),
            "saveToTree": self.saveToTreeInput.isChecked(),
            "SmearFactorconst": self.SmearFactorconstInput.value(),
            "nSteps": self.nStepsInput.value(),
            "stepSize": self.stepSizeInput.value(),
            "outputFormat": self.outputFormatInput.currentText(),
        }
        with open(self.params_file, "w") as file:
            json.dump(params, file)

    def getRunCounter(self):
        if os.path.exists(self.counter_file):
            with open(self.counter_file, "r") as file:
                counter = json.load(file).get("counter", 0)
        else:
            counter = 0
        return counter

    def incrementRunCounter(self):
        counter = self.getRunCounter() + 1
        with open(self.counter_file, "w") as file:
            json.dump({"counter": counter}, file)
        return counter

    def runSimulation(self):
        self.saveParams()
        run_number = self.incrementRunCounter()

        particleType = self.particleTypeInput.currentText().capitalize()
        nParticles = self.nParticlesInput.value()
        ptMax = self.ptMaxInput.value()
        ptMin = self.ptMinInput.value()
        weightMethod = self.weightMethodInput.currentText()
        applyAsymmCut = "true" if self.asymmCutInput.isChecked() else "false"
        asymmCutValue = self.asymmCutValueInput.value()
        clusterOverlap = "true" if self.clusterOverlapInput.isChecked() else "false"
        clusterOverlapProb = self.clusterOverlapProbInput.value()
        deltaRcutMax = self.deltaRcutMaxInput.value()
        pt1cut = self.pt1cutInput.value()
        pt2cut = self.pt2cutInput.value()
        combPtCut = self.combPtCutInput.value()
        ptMaxCut = self.ptMaxCutInput.value()
        nclusPtCut = self.nclusPtCutInput.value()
        positSmearingFactor = self.positSmearingFactorInput.value()
        saveToTree = "true" if self.saveToTreeInput.isChecked() else "false"
        SmearFactorconst = self.SmearFactorconstInput.value()
        nSteps = self.nStepsInput.value()
        stepSize = self.stepSizeInput.value()
        PT_Max_bin = self.PT_Max_binInput.value()
        MassNBins = self.MassNBinsInput.value()
        binres = self.binresInput.value()
        n_bins = self.n_binsInput.value()
        Debug = "true" if self.DebugInput.isChecked() else "false"
        etCut = self.etCutInput.value()
        Apply_Eta_Cut = "true" if self.Apply_Eta_CutInput.isChecked() else "false"
        eta_cut_val = self.eta_cut_valInput.value()
        smeared_lower_bin_limit = self.smeared_lower_bin_limitInput.value()
        smeared_upper_bin_limit = self.smeared_upper_bin_limitInput.value()
        smear_factor_sqrtE = self.smear_factor_sqrtEInput.value()
        outputFormat = self.outputFormatInput.currentText().lower()

        script = "gen_res"  # Replace this with your actual script name

        command = (
            f"./{script} -particleType={particleType} -nParticles={nParticles} -PT_Max={ptMax} "
            f"-PT_Min={ptMin} -weightMethod={weightMethod} -applyAsymmCut={applyAsymmCut} "
            f"-asymmCutValue={asymmCutValue} -clusterOverlap={clusterOverlap} "
            f"-clusterOverlapProb={clusterOverlapProb} -DeltaRcut_MAX={deltaRcutMax} "
            f"-pt1cut={pt1cut} -pt2cut={pt2cut} -comb_ptcut={combPtCut} -ptMaxCut={ptMaxCut} "
            f"-nclus_ptCut={nclusPtCut} -posit_smearingFactor={positSmearingFactor} "
            f"-saveToTree={saveToTree} -smear_factor_const={SmearFactorconst} "
            f"-smear_factor_const_step_size={stepSize} -smear_factor_const_num_steps={nSteps} "
            f"-PT_Max_bin={PT_Max_bin} -MassNBins={MassNBins} "
            f"-binres={binres} -n_bins={n_bins} "
            f"-Debug={Debug} -etCut={etCut} "
            f"-Apply_Eta_Cut={Apply_Eta_Cut} -eta_cut_val={eta_cut_val} "
            f"-smeared_lower_bin_limit={smeared_lower_bin_limit} "
            f"-smeared_upper_bin_limit={smeared_upper_bin_limit} "
            f"-smear_factor_sqrtE={smear_factor_sqrtE}"
        )

        print("Running command:", command)
        os.system(command)

        params = {
            "Particle Type": particleType,
            "Number of Particles": nParticles,
            "PT Max (GeV)": ptMax,
            "PT Min (GeV)": ptMin,
            "Weight Method": weightMethod,
            "Apply Asymm Cut": applyAsymmCut,
            "Asymm Cut Value": asymmCutValue,
            "Cluster Overlap": clusterOverlap,
            "Cluster Overlap Probability": clusterOverlapProb,
            "Delta R Cut Max": deltaRcutMax,
            "PT1 Cut (GeV)": pt1cut,
            "PT2 Cut (GeV)": pt2cut,
            "Combined PT Cut factor": combPtCut,
            "PT Max Cut (GeV)": ptMaxCut,
            "NCLUS PT Cut (GeV)": nclusPtCut,
            "Position Smearing Factor (mm)": positSmearingFactor,
            "Save to Tree": saveToTree,
            "Base Constant Smearing (%)": SmearFactorconst,
            "Number of Steps": nSteps,
            "Step Size": stepSize,
            "PT Max Bin": PT_Max_bin,
            "Mass N Bins": MassNBins,
            "Bin Resolution": binres,
            "Number of Bins": n_bins,
            "Debug": Debug,
            "ET Cut (GeV)": etCut,
            "Apply Eta Cut": Apply_Eta_Cut,
            "Eta Cut Value": eta_cut_val,
            "Smeared Lower Bin Limit": smeared_lower_bin_limit,
            "Smeared Upper Bin Limit": smeared_upper_bin_limit,
            "Smear Factor sqrt(E)": smear_factor_sqrtE,
            "Output Format": outputFormat,
            "Run Number": run_number,
            "Timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        if not os.path.exists("pioncode/canvas_pdf"):
            os.makedirs("pioncode/canvas_pdf")

        filename = f"pioncode/canvas_pdf/simulation_params_{run_number}.{outputFormat}"
        if outputFormat == "pdf":
            self.printParametersToPDF(params, filename)
        elif outputFormat in ["png", "jpeg"]:
            self.printParametersToImage(params, filename)

    def printParametersToPDF(self, params, filePath):
        if not os.path.exists("pioncode/canvas_pdf"):
            os.makedirs("pioncode/canvas_pdf")

        c = canvas.Canvas(filePath, pagesize=letter)
        width, height = letter
        c.drawString(100, height - 40, "Simulation Parameters")
        y = height - 80
        for key, value in params.items():
            c.drawString(100, y, f"{key}: {value}")
            y -= 20
        c.save()
        print(f"Parameters saved to {filePath}")

    def printParametersToImage(self, params, filePath):
        if not os.path.exists("pioncode/canvas_pdf"):
            os.makedirs("pioncode/canvas_pdf")

        fig, ax = plt.subplots(figsize=(8.5, 11))
        ax.axis("off")
        text = "\n".join([f"{key}: {value}" for key, value in params.items()])
        ax.text(
            0.5,
            0.5,
            text,
            transform=ax.transAxes,
            fontsize=12,
            va="center",
            ha="center",
            wrap=True,
        )
        plt.savefig(filePath, bbox_inches="tight")
        plt.close()
        print(f"Parameters saved to {filePath}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = SimulationGUI()
    ex.show()
    sys.exit(app.exec_())
