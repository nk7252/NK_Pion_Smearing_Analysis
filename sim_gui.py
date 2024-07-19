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
        self.resize(550, 800)  # Set the initial size of the GUI
        self.setMinimumSize(550, 800)  # Set the minimum size of the GUI
        self.initUI()

    def initUI(self):
        self.main_layout = QVBoxLayout(self)

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        self.main_layout.addWidget(scroll_area)

        content_widget = QWidget()
        self.layout = QVBoxLayout(content_widget)

        scroll_area.setWidget(content_widget)

        self.particleTypeInput = self.addFieldWithExplanation(
            "Particle Type:",
            QComboBox(self),
            explanation="Select the type of particle to simulate.",
        )
        self.particleTypeInput.addItems(["Pion", "Eta"])

        self.nParticlesInput = self.addFieldWithExplanation(
            "Number of Particles:",
            QSpinBox(self),
            max_value=100000000,
            default_value=8000000,
            explanation="Enter the number of particles to simulate.",
        )

        self.ptMaxInput = self.addFieldWithExplanation(
            "PT Max (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=50,
            explanation="Enter the maximum transverse momentum (PT) in GeV.",
        )

        self.ptMinInput = self.addFieldWithExplanation(
            "PT Min (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=0,
            explanation="Enter the minimum transverse momentum (PT) in GeV.",
        )

        self.weightMethodInput = self.addFieldWithExplanation(
            "Weight Method:",
            QComboBox(self),
            explanation="Select the weighting method.",
        )
        self.weightMethodInput.addItems(["EXP", "POWER", "WSHP", "HAGEDORN"])

        self.asymmCutInput = self.addFieldWithExplanation(
            "Apply Asymm Cut:",
            QCheckBox(self),
            default_value=True,
            explanation="Check to apply asymmetric cut.",
        )

        self.asymmCutValueInput = self.addFieldWithExplanation(
            "Asymm Cut Value:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=3,
            default_value=0.6,
            explanation="Enter the value for the asymmetric cut.",
        )

        self.clusterOverlapInput = self.addFieldWithExplanation(
            "Cluster Overlap:",
            QCheckBox(self),
            default_value=True,
            explanation="Check to enable cluster overlap.",
        )

        self.clusterOverlapProbInput = self.addFieldWithExplanation(
            "Cluster Overlap Probability:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=3,
            step=0.001,
            default_value=0.99,
            explanation="Enter the probability of NO cluster overlap. i.e. 0.99 means 1% overlap probability.",
        )

        self.deltaRcutMaxInput = self.addFieldWithExplanation(
            "Delta R Cut Max:",
            QDoubleSpinBox(self),
            max_value=10,
            decimals=2,
            default_value=1.1,
            explanation="Enter the maximum delta R cut value.",
        )

        self.pt1cutInput = self.addFieldWithExplanation(
            "PT1 Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=1.5,
            explanation="Enter the cut value for PT1 in GeV.",
        )

        self.pt2cutInput = self.addFieldWithExplanation(
            "PT2 Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=1.5,
            explanation="Enter the cut value for PT2 in GeV.",
        )

        self.combPtCutInput = self.addFieldWithExplanation(
            "Combined PT Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=0,
            explanation="Enter the combined PT cut value in GeV.",
        )

        self.ptMaxCutInput = self.addFieldWithExplanation(
            "PT Max Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=50,
            explanation="Enter the maximum PT cut value in GeV.",
        )

        self.nclusPtCutInput = self.addFieldWithExplanation(
            "NCLUS PT Cut (GeV):",
            QDoubleSpinBox(self),
            max_value=100,
            decimals=2,
            default_value=0.0,
            explanation="Enter the NCLUS PT cut value in GeV.",
        )

        self.positSmearingFactorInput = self.addFieldWithExplanation(
            "Position Smearing Factor (mm):",
            QDoubleSpinBox(self),
            max_value=10,
            decimals=2,
            default_value=2.8,
            explanation="Enter the position smearing factor.",
        )

        self.saveToTreeInput = self.addFieldWithExplanation(
            "Save to Tree:",
            QCheckBox(self),
            default_value=False,
            explanation="Check to save the results to a tree.",
        )

        self.baseSmearPercentInput = self.addFieldWithExplanation(
            "Base Smearing Percent:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=3,
            step=0.001,
            default_value=0.168,
            explanation="Enter the base smearing percent.",
        )

        self.nStepsInput = self.addFieldWithExplanation(
            "Number of Steps:",
            QSpinBox(self),
            max_value=100,
            default_value=1,
            explanation="Enter the number of steps for the simulation.",
        )

        self.stepSizeInput = self.addFieldWithExplanation(
            "Step Size:",
            QDoubleSpinBox(self),
            max_value=1,
            decimals=4,
            step=0.0001,
            default_value=0.001,
            explanation="Enter the step size for the simulation.",
        )

        self.outputFormatInput = self.addFieldWithExplanation(
            "Output Format:",
            QComboBox(self),
            explanation="Select the format for saving the output.",
        )
        self.outputFormatInput.addItems(["PDF", "PNG", "JPEG"])

        self.runButton = QPushButton("Run Simulation", self)
        self.runButton.clicked.connect(self.runSimulation)
        self.layout.addWidget(self.runButton)

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
        self.layout.addWidget(container)

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
                self.particleTypeInput.setCurrentText(
                    params.get("particleType", "Pion")
                )
                self.nParticlesInput.setValue(params.get("nParticles", 8000000))
                self.ptMaxInput.setValue(params.get("ptMax", 50))
                self.ptMinInput.setValue(params.get("ptMin", 0))
                self.weightMethodInput.setCurrentText(
                    params.get("weightMethod", "WSHP")
                )
                self.asymmCutInput.setChecked(params.get("applyAsymmCut", True))
                self.asymmCutValueInput.setValue(params.get("asymmCutValue", 0.6))
                self.clusterOverlapInput.setChecked(params.get("clusterOverlap", True))
                self.clusterOverlapProbInput.setValue(
                    params.get("clusterOverlapProb", 0.99)
                )
                self.deltaRcutMaxInput.setValue(params.get("deltaRcutMax", 1.1))
                self.pt1cutInput.setValue(params.get("pt1cut", 1.5))
                self.pt2cutInput.setValue(params.get("pt2cut", 1.5))
                self.combPtCutInput.setValue(params.get("combPtCut", 0))
                self.ptMaxCutInput.setValue(params.get("ptMaxCut", 50))
                self.nclusPtCutInput.setValue(params.get("nclusPtCut", 0.0))
                self.positSmearingFactorInput.setValue(
                    params.get("positSmearingFactor", 2.8)
                )
                self.saveToTreeInput.setChecked(params.get("saveToTree", False))
                self.baseSmearPercentInput.setValue(
                    params.get("baseSmearPercent", 0.168)
                )
                self.nStepsInput.setValue(params.get("nSteps", 25))
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
            "nclusPtCut": self.nclusPtCutInput.value(),
            "positSmearingFactor": self.positSmearingFactorInput.value(),
            "saveToTree": self.saveToTreeInput.isChecked(),
            "baseSmearPercent": self.baseSmearPercentInput.value(),
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
        baseSmearPercent = self.baseSmearPercentInput.value()
        nSteps = self.nStepsInput.value()
        stepSize = self.stepSizeInput.value()
        outputFormat = self.outputFormatInput.currentText().lower()

        script = "gen_res"

        command = (
            f"./{script} -particleType={particleType} -nParticles={nParticles} -PT_Max={ptMax} "
            f"-PT_Min={ptMin} -weightMethod={weightMethod} -applyAsymmCut={applyAsymmCut} "
            f"-asymmCutValue={asymmCutValue} -clusterOverlap={clusterOverlap} "
            f"-clusterOverlapProb={clusterOverlapProb} -DeltaRcut_MAX={deltaRcutMax} "
            f"-pt1cut={pt1cut} -pt2cut={pt2cut} -comb_ptcut={combPtCut} -ptMaxCut={ptMaxCut} "
            f"-nclus_ptCut={nclusPtCut} -posit_smearingFactor={positSmearingFactor} "
            f"-saveToTree={saveToTree} -smear_factor_basevalue={baseSmearPercent} "
            f"-smear_factor_step={stepSize} -smear_factor_steps={nSteps}"
        )

        print("Running command:", command)

        os.system(command)

        params = {
            "Particle Type": self.particleTypeInput.currentText(),
            "Number of Particles": self.nParticlesInput.value(),
            "PT Max (GeV)": self.ptMaxInput.value(),
            "PT Min (GeV)": self.ptMinInput.value(),
            "Weight Method": self.weightMethodInput.currentText(),
            "Apply Asymm Cut": self.asymmCutInput.isChecked(),
            "Asymm Cut Value": self.asymmCutValueInput.value(),
            "Cluster Overlap": self.clusterOverlapInput.isChecked(),
            "Cluster Overlap Probability": self.clusterOverlapProbInput.value(),
            "Delta R Cut Max": self.deltaRcutMaxInput.value(),
            "PT1 Cut (GeV)": self.pt1cutInput.value(),
            "PT2 Cut (GeV)": self.pt2cutInput.value(),
            "Combined PT Cut factor": self.combPtCutInput.value(),
            "PT Max Cut (GeV)": self.ptMaxCutInput.value(),
            "NCLUS PT Cut (GeV)": self.nclusPtCutInput.value(),
            "Position Smearing Factor (mm)": self.positSmearingFactorInput.value(),
            "Save to Tree": self.saveToTreeInput.isChecked(),
            "Base Smearing Percent": self.baseSmearPercentInput.value(),
            "Number of Steps": self.nStepsInput.value(),
            "Step Size": self.stepSizeInput.value(),
            "Output Format": self.outputFormatInput.currentText(),
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
