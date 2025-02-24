import os, csv
import vtk, qt, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
import sitkUtils

import numpy as np
import SimpleITK as sitk
import pydicom
import math
from scipy import ndimage, stats, spatial
from mahotas.features.shape import roundness

try:
    import matplotlib

except ModuleNotFoundError:
    slicer.util.pip_install("matplotlib")
    import matplotlib

import matplotlib.pyplot as plt
import matplotlib.colors as cls

matplotlib.use("Agg")


from QAHybridLib import Auxiliar_functions as af


# QAHybrid
#


class QAHybrid(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "QAHybrid"
        self.parent.categories = ["Quantification"]
        self.parent.dependencies = []
        self.parent.contributors = [
            "Tobias Fechter (Department of Radiation Oncology, University Medical Centre Freiburg), Carmen Salvador Ribes (Biomedical Imaging Research Group, La Fe Health Research Institute), Montserrat Carles Farina (Biomedical Imaging Research Group, La Fe Health Research Institute)"
        ]
        self.parent.helpText = """
This extension enables the analysis of different quality  parameters for PET, CT and MR images of the Quality-Assurance section within the Hybrid-phantom: quantification, resolution, registration, distortion and radiomics.
"""
        self.parent.acknowledgementText = """
"""


#
# QAHybridWidget
#


class QAHybridWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)  # needed for parameter node observation
        self.logic = None
        self._parameterNode = None
        self._updatingGUIFromParameterNode = False
        self._updatePetGuiValues = False
        self._updateInsertTypeValues = False
        self.timer = None

    def onQuantSliceCTButton(self):
        self.logic.reportProgress(10, "CT Slice detection")
        centerSliceIdx = QuantificationRoutine.getCenterSlicerForCtInsert(self.ui.ctSelector.currentNode())
        self.logic.reportProgress(80, "CT Slice detection")
        self.ui.firstSliceCTSpinBox.value = centerSliceIdx
        self.logic.reportProgress(100, "CT Slice detection")

    def onQuantSliceT1Button(self):
        pass

    def onQuantSliceT2Button(self):
        pass

    def onQuantSlicePETButton(self):
        self.logic.reportProgress(10, "PET Slice detection")
        centerSliceIdx = QuantificationRoutine.getCenterSlicerForPetInsert(self.ui.petSelector.currentNode())
        self.logic.reportProgress(80, "PET Slice detection")
        self.ui.firstSlicePETSpinBox.value = centerSliceIdx
        self.logic.reportProgress(100, "PET Slice detection")

    def onDetectDistortionSlicesCTButton(self):
        self.logic.reportProgress(10, "Distortion Insert CT detection")
        startSlice, endSlice, panalSlice = self.logic.distortionStartEndSliceCylinderDetection(
            self.ui.ctSelector.currentNode(), "CT", self.resourcePath("RadiomicsData"), shrinkImage=True
        )
        self.logic.reportProgress(80, "Distortion Insert CT detection")
        self.ui.distortionFirstSliceCTSpinBox.value = startSlice
        self.ui.distortionLastSliceCTSpinBox.value = endSlice
        self.ui.panalSliceCTSpinBox.value = panalSlice
        self.logic.reportProgress(100, "Distortion Insert CT detection")

    def onDetectDistortionSlicesT1Button(self):
        self.logic.reportProgress(10, "Distortion Insert T1 detection")
        startSlice, endSlice, panalSlice = self.logic.distortionStartEndSliceCylinderDetection(
            self.ui.t1Selector.currentNode(), "T1", self.resourcePath("RadiomicsData")
        )
        self.logic.reportProgress(80, "Distortion Insert T1 detection")
        self.ui.distortionFirstSliceT1SpinBox.value = startSlice
        self.ui.distortionLastSliceT1SpinBox.value = endSlice
        self.ui.panalSliceT1SpinBox.value = panalSlice
        self.logic.reportProgress(100, "Distortion Insert T1 detection")

    def onDetectDistortionSlicesT2Button(self):
        self.logic.reportProgress(10, "Distortion Insert T2 detection")
        startSlice, endSlice, panalSlice = self.logic.distortionStartEndSliceCylinderDetection(
            self.ui.t2Selector.currentNode(), "T2", self.resourcePath("RadiomicsData"), shrinkImage=True
        )
        self.logic.reportProgress(80, "Distortion Insert T2 detection")
        self.ui.distortionFirstSliceT2SpinBox.value = startSlice
        self.ui.distortionLastSliceT2SpinBox.value = endSlice
        self.ui.panalSliceT2SpinBox.value = panalSlice
        self.logic.reportProgress(100, "Distortion Insert T2 detection")

    def onDetectResolutionSliceCTButton(self):
        self.logic.reportProgress(10, "Resolution Insert CT detection")
        self.logic.resolutionSliceDetection(
            self.ui.ctSelector.currentNode(), "CT", self.resourcePath("RadiomicsData"), shrinkImage=True
        )
        self.logic.reportProgress(80, "Resolution Insert CT detection")
        self.logic.reportProgress(100, "Resolution Insert CT detection")
        self.updateGUIFromParameterNode()

    def onDetectResolutionSlicePETButton(self):
        self.logic.reportProgress(10, "Resolution Insert PET detection")
        self.logic.resolutionSliceDetection(
            self.ui.petSelector.currentNode(),
            "PET",
            self.resourcePath("RadiomicsData"),
        )
        self.logic.reportProgress(80, "Resolution Insert PET detection")
        self.logic.reportProgress(100, "Resolution Insert PET detection")
        self.updateGUIFromParameterNode()

    def onDetectResolutionSliceT1Button(self):
        self.logic.reportProgress(10, "Resolution Insert T1 detection")
        self.logic.resolutionSliceDetection(
            self.ui.t1Selector.currentNode(),
            "T1",
            self.resourcePath("RadiomicsData"),
        )
        self.logic.reportProgress(80, "Resolution Insert T1 detection")
        self.logic.reportProgress(100, "Resolution Insert T1 detection")
        self.updateGUIFromParameterNode()

    def onDetectResolutionSliceT2Button(self):
        self.logic.reportProgress(10, "Resolution Insert T2 detection")
        self.logic.resolutionSliceDetection(
            self.ui.t2Selector.currentNode(), "T2", self.resourcePath("RadiomicsData"), shrinkImage=True
        )
        self.logic.reportProgress(80, "Resolution Insert T2 detection")
        self.logic.reportProgress(100, "Resolution Insert T2 detection")
        self.updateGUIFromParameterNode()

    def setup(self):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/QAHybrid.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = QAHybridLogic()

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
        # (in the selected parameter node).
        self.ui.petSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.ctSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.t1Selector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.t2Selector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)

        self.ui.analysisTypeSelector.addItems(list(self.logic.getSupportedAnalysesMethods().keys()))
        self.ui.analysisTypeSelector.connect("checkedIndexesChanged()", self.updateParameterNodeFromGUI)

        self.ui.insertTypeSelector.addItems(list(self.logic.getSupportedInserts().keys()))
        self.ui.insertTypeSelector.connect("currentIndexChanged(int)", self.updateParameterNodeFromGUI)

        self.ui.tracerSelector.addItems(list(self.logic.getSupportedTracers()))
        self.ui.tracerSelector.connect("currentIndexChanged(int)", self.updateParameterNodeFromGUI)

        # Buttons
        self.ui.applyButton.connect("clicked(bool)", self.onApplyButton)

        self.ui.detectQuantSlicePETButton.connect("clicked(bool)", self.onQuantSlicePETButton)
        self.ui.detectQuantSliceCTButton.connect("clicked(bool)", self.onQuantSliceCTButton)
        self.ui.detectQuantSliceT1Button.connect("clicked(bool)", self.onQuantSliceT1Button)
        self.ui.detectQuantSliceT2Button.connect("clicked(bool)", self.onQuantSliceT2Button)

        self.ui.detectResolutionSliceCTButton.connect("clicked(bool)", self.onDetectResolutionSliceCTButton)
        self.ui.detectResolutionSlicePETButton.connect("clicked(bool)", self.onDetectResolutionSlicePETButton)
        self.ui.detectResolutionSliceT1Button.connect("clicked(bool)", self.onDetectResolutionSliceT1Button)
        self.ui.detectResolutionSliceT2Button.connect("clicked(bool)", self.onDetectResolutionSliceT2Button)

        self.ui.distortionInsertCTButton.connect("clicked(bool)", self.onDetectDistortionSlicesCTButton)
        self.ui.distortionInsertT1Button.connect("clicked(bool)", self.onDetectDistortionSlicesT1Button)
        self.ui.distortionInsertT2Button.connect("clicked(bool)", self.onDetectDistortionSlicesT2Button)

        # input numbers
        self.ui.rcVolume.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
        self.ui.decayFactor.connect("valueChanged(double)", self.updateParameterNodeFromGUI)

        self.ui.initialActivity.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.timeActivity.connect("valueChanged(int)", self.updateParameterNodeFromGUI)

        self.ui.firstSlicePETSpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.firstSliceCTSpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.firstSliceT1SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.firstSliceT2SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)

        self.ui.distortionFirstSliceCTSpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.distortionLastSliceCTSpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.panalSliceCTSpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)

        self.ui.distortionFirstSliceT1SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.distortionLastSliceT1SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.panalSliceT1SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.distortionFirstSliceT2SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.distortionLastSliceT2SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)
        self.ui.panalSliceT2SpinBox.connect("valueChanged(int)", self.updateParameterNodeFromGUI)

        # input fields positions
        self.ui.p1.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p2.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p3.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p4.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p5.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p6.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p7.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p8.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p9.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p10.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.p11.connect("valueChanged(int)", self.updateParameterNodeFromGUIWithTimer)
        self.ui.l1.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l2.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l3.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l4.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l5.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l6.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l7.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l8.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l9.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l10.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)
        self.ui.l11.textEdited.connect(self.updateParameterNodeFromGUIWithTimer)

        self.ui.defaultRadiomicSegmentPathButton.directory = self.resourcePath("RadiomicsData")

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

    def cleanup(self):
        """
        Called when the application closes and the module widget is destroyed.
        """
        self.removeObservers()

    def enter(self):
        """
        Called each time the user opens this module.
        """
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self):
        """
        Called each time the user opens a different module.
        """
        # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
        self.removeObserver(
            self._parameterNode,
            vtk.vtkCommand.ModifiedEvent,
            self.updateGUIFromParameterNode,
        )

    def onSceneStartClose(self, caller, event):
        """
        Called just before the scene is closed.
        """
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event):
        """
        Called just after the scene is closed.
        """
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self):
        """
        Ensure parameter node exists and observed.
        """
        # Parameter node stores all user choices in parameter values, node selections, etc.
        # so that when the scene is saved and reloaded, these settings are restored.

        self.setParameterNode(self.logic.getParameterNode())

        # Select default input nodes if nothing is selected yet to save a few clicks for the user
        # if not self._parameterNode.GetNodeReference("PetInputVolume"):
        #   firstVolumeNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLScalarVolumeNode")
        #   if firstVolumeNode:
        #     self._parameterNode.SetNodeReferenceID("PetInputVolume", firstVolumeNode.GetID())

    def setParameterNode(self, inputParameterNode):
        """
        Set and observe parameter node.
        Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
        """

        if inputParameterNode:
            self.logic.setDefaultParameters(inputParameterNode)

        # Unobserve previously selected parameter node and add an observer to the newly selected.
        # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
        # those are reflected immediately in the GUI.
        if self._parameterNode is not None:
            self.removeObserver(
                self._parameterNode,
                vtk.vtkCommand.ModifiedEvent,
                self.updateGUIFromParameterNode,
            )
        self._parameterNode = inputParameterNode
        if self._parameterNode is not None:
            self.addObserver(
                self._parameterNode,
                vtk.vtkCommand.ModifiedEvent,
                self.updateGUIFromParameterNode,
            )

        # Initial GUI update
        self.updateGUIFromParameterNode()

    def updateGUIFromParameterNode(self, caller=None, event=None):
        """
        This method is called whenever parameter node is changed.
        The module GUI is updated to show the current state of the parameter node.
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
        self._updatingGUIFromParameterNode = True

        self.ui.warningTextBrowser.clear()

        if self._parameterNode.GetParameter("insertTypeSelector"):
            self.ui.insertTypeSelector.setCurrentText(self._parameterNode.GetParameter("insertTypeSelector"))

        if self._parameterNode.GetParameter("decayFactor"):
            self.ui.decayFactor.setValue(float(self._parameterNode.GetParameter("decayFactor")))

        if self._parameterNode.GetParameter("rcVolume"):
            self.ui.rcVolume.setValue(float(self._parameterNode.GetParameter("rcVolume")))

        if self._parameterNode.GetParameter("tracerSelector"):
            self.ui.tracerSelector.setCurrentText(self._parameterNode.GetParameter("tracerSelector"))

        if self._parameterNode.GetParameter("initialActivity"):
            self.ui.initialActivity.setValue(int(self._parameterNode.GetParameter("initialActivity")))

        if self._parameterNode.GetParameter("timeActivity"):
            self.ui.timeActivity.setValue(int(self._parameterNode.GetParameter("timeActivity")))

        if self._updateInsertTypeValues:
            self._updateInsertTypeValues = False
            self.logic.setInsertTypePositions(
                self._parameterNode.GetParameter("insertTypeSelector"), self._parameterNode
            )

        for insertPosKey in self.logic.materialCharacteristicsPosition:
            if self._parameterNode.GetParameter("p" + str(insertPosKey)):
                uiElement = getattr(self.ui, "p" + str(insertPosKey))
                uiElement.setValue(int(self._parameterNode.GetParameter("p" + str(insertPosKey))))
            if self._parameterNode.GetParameter("l" + str(insertPosKey)):
                uiElement = getattr(self.ui, "l" + str(insertPosKey))
                uiElement.setText(self._parameterNode.GetParameter("l" + str(insertPosKey)))

        # Update node selectors and sliders
        if not self._checkPropertiesOfSelectedPetNode(self._parameterNode.GetNodeReference("PetInputVolume")):
            self.ui.petSelector.setCurrentNode(self._parameterNode.GetNodeReference("PetInputVolume"))

            if self._updatePetGuiValues:
                self._updatePetGuiValues = False

                if self._parameterNode.GetNodeReference("PetInputVolume") is None:
                    self.ui.firstSlicePETSpinBox.setEnabled(False)
                    self.ui.detectQuantSlicePETButton.setEnabled(False)
                    self.ui.detectResolutionSlicePETButton.setEnabled(False)
                else:
                    self.ui.firstSlicePETSpinBox.setEnabled(True)
                    self.ui.detectResolutionSlicePETButton.setEnabled(True)
                    self.ui.detectQuantSlicePETButton.setEnabled(True)
                    numberOfSlices = (
                        self._parameterNode.GetNodeReference("PetInputVolume").GetImageData().GetDimensions()[2]
                    )
                    self.ui.firstSlicePETSpinBox.maximum = numberOfSlices - 1

                time_activity = self.logic.getDicomTagOrSequence(
                    self._parameterNode.GetNodeReference("PetInputVolume"), 0x0054, 0x0016, 0, 0x0018, 0x1072
                )
                if time_activity is not None:
                    self.ui.timeActivity.setValue(int(time_activity.value))

                radionuclideTotalDose = self.logic.getDicomTagOrSequence(
                    self._parameterNode.GetNodeReference("PetInputVolume"), 0x0054, 0x0016, 0, 0x0018, 0x1074
                )
                if radionuclideTotalDose is not None:
                    self.ui.initialActivity.setValue(int(radionuclideTotalDose.value))

                radiopharmaceutical = self.logic.getDicomTagOrSequence(
                    self._parameterNode.GetNodeReference("PetInputVolume"), 0x0054, 0x0016, 0, 0x0018, 0x0031
                )
                if radiopharmaceutical is not None:
                    if "fluor" in radiopharmaceutical.value.lower() and "18" in radiopharmaceutical.value.lower():
                        self.ui.tracerSelector.setCurrentText(self.logic.getSupportedTracers()[0])
                        self.ui.decayFactor.setValue(self.logic.decayFactor)

        else:
            self.ui.firstSlicePETSpinBox.setEnabled(False)
            self.ui.petSelector.setCurrentNode(None)
            self.ui.warningTextBrowser.setText(
                "Warning: selected PET Volume does not contain all necessary Information. Hint: Load PET Volume with Dicom-Browser."
            )

        if self._checkPropertiesOfSelectedCtNode(self._parameterNode.GetNodeReference("CtInputVolume")):
            self.ui.firstSliceCTSpinBox.setEnabled(True)
            self.ui.detectQuantSliceCTButton.setEnabled(True)

            self.ui.detectResolutionSliceCTButton.setEnabled(True)

            numberOfSlices = self._parameterNode.GetNodeReference("CtInputVolume").GetImageData().GetDimensions()[2]
            self.ui.firstSliceCTSpinBox.maximum = numberOfSlices - 1
            self.ui.ctSelector.setCurrentNode(self._parameterNode.GetNodeReference("CtInputVolume"))

            self.ui.distortionFirstSliceCTSpinBox.setEnabled(True)
            self.ui.distortionFirstSliceCTSpinBox.maximum = numberOfSlices - 1

            self.ui.distortionLastSliceCTSpinBox.setEnabled(True)
            self.ui.distortionLastSliceCTSpinBox.maximum = numberOfSlices - 1

            self.ui.panalSliceCTSpinBox.setEnabled(True)
            self.ui.panalSliceCTSpinBox.maximum = numberOfSlices - 1

            self.ui.distortionInsertCTButton.setEnabled(True)

        else:
            self.ui.ctSelector.setCurrentNode(None)
            self.ui.firstSliceCTSpinBox.setEnabled(False)
            self.ui.detectQuantSliceCTButton.setEnabled(False)
            self.ui.detectResolutionSliceCTButton.setEnabled(False)

            self.ui.distortionFirstSliceCTSpinBox.setEnabled(False)
            self.ui.distortionLastSliceCTSpinBox.setEnabled(False)

            self.ui.panalSliceCTSpinBox.setEnabled(False)
            self.ui.distortionInsertCTButton.setEnabled(False)

            self.ui.warningTextBrowser.setText(
                "Warning: selected CT Volume does not contain all necessary Information. Hint: Load CT Volume with Dicom-Browser."
            )

        if self._checkPropertiesOfSelectedMRNode(self._parameterNode.GetNodeReference("T1InputVolume")):
            self.ui.firstSliceT1SpinBox.setEnabled(True)
            self.ui.detectQuantSliceT1Button.setEnabled(False)
            numberOfSlices = self._parameterNode.GetNodeReference("T1InputVolume").GetImageData().GetDimensions()[2]
            self.ui.firstSliceT1SpinBox.maximum = numberOfSlices - 1

            self.ui.detectResolutionSliceT1Button.setEnabled(True)

            self.ui.distortionFirstSliceT1SpinBox.setEnabled(True)
            self.ui.distortionFirstSliceT1SpinBox.maximum = numberOfSlices - 1

            self.ui.distortionLastSliceT1SpinBox.setEnabled(True)
            self.ui.distortionLastSliceT1SpinBox.maximum = numberOfSlices - 1

            self.ui.panalSliceT1SpinBox.setEnabled(True)
            self.ui.panalSliceT1SpinBox.maximum = numberOfSlices - 1

            self.ui.distortionInsertT1Button.setEnabled(True)

            self.ui.t1Selector.setCurrentNode(self._parameterNode.GetNodeReference("T1InputVolume"))
        else:
            self.ui.firstSliceT1SpinBox.setEnabled(False)
            self.ui.detectQuantSliceT1Button.setEnabled(False)
            self.ui.detectResolutionSliceT1Button.setEnabled(False)

            self.ui.distortionFirstSliceT1SpinBox.setEnabled(False)
            self.ui.distortionInsertT1Button.setEnabled(False)

            self.ui.panalSliceT1SpinBox.setEnabled(False)

            self.ui.distortionLastSliceT1SpinBox.setEnabled(False)

            self.ui.t1Selector.setCurrentNode(None)
            self.ui.warningTextBrowser.setText(
                "Warning: selected T1 Volume does not contain all necessary Information. Hint: Load T1 Volume with Dicom-Browser."
            )

        if self._checkPropertiesOfSelectedMRNode(self._parameterNode.GetNodeReference("T2InputVolume")):
            self.ui.firstSliceT2SpinBox.setEnabled(True)
            self.ui.detectQuantSliceT2Button.setEnabled(False)
            numberOfSlices = self._parameterNode.GetNodeReference("T2InputVolume").GetImageData().GetDimensions()[2]
            self.ui.firstSliceT2SpinBox.maximum = numberOfSlices - 1
            self.ui.t2Selector.setCurrentNode(self._parameterNode.GetNodeReference("T2InputVolume"))

            self.ui.detectResolutionSliceT2Button.setEnabled(True)
            self.ui.distortionLastSliceT2SpinBox.setEnabled(True)
            self.ui.distortionLastSliceT2SpinBox.maximum = numberOfSlices - 1
            self.ui.distortionFirstSliceT2SpinBox.setEnabled(True)
            self.ui.distortionFirstSliceT2SpinBox.maximum = numberOfSlices - 1
            self.ui.distortionInsertT2Button.setEnabled(True)

            self.ui.panalSliceT2SpinBox.setEnabled(True)
            self.ui.panalSliceT2SpinBox.maximum = numberOfSlices - 1
        else:
            self.ui.firstSliceT2SpinBox.setEnabled(False)
            self.ui.detectQuantSliceT2Button.setEnabled(False)
            self.ui.detectResolutionSliceT2Button.setEnabled(False)
            self.ui.distortionFirstSliceT2SpinBox.setEnabled(False)
            self.ui.distortionLastSliceT1SpinBox.setEnabled(False)
            self.ui.distortionInsertT2Button.setEnabled(False)
            self.ui.distortionLastSliceT2SpinBox.setEnabled(False)
            self.ui.panalSliceT2SpinBox.setEnabled(False)
            self.ui.t2Selector.setCurrentNode(None)
            self.ui.warningTextBrowser.setText(
                "Warning: selected T2 Volume does not contain all necessary Information. Hint: Load T2 Volume with Dicom-Browser."
            )

        # if not self._checkFrameOfReferenceUID(
        #     self._parameterNode.GetNodeReference("PetInputVolume"),
        #     self._parameterNode.GetNodeReference("CtInputVolume"),
        # ):
        #     self.ui.warningTextBrowser.setText(
        #         "Warning: selected CT and PET Volumes do not share the same frame of reference."
        #     )

        idx = 0
        methodChecked = False
        selectedIdxs = []
        for analysesMethod in self.logic.getSupportedAnalysesMethods().keys():
            if (
                self._parameterNode.GetParameter(analysesMethod)
                and self._parameterNode.GetParameter(analysesMethod) == "1"
            ):
                self.ui.analysisTypeSelector.model().item(idx).setCheckState(2)
                methodChecked = True
                selectedIdxs.append(idx)
            else:
                self.ui.analysisTypeSelector.model().item(idx).setCheckState(0)
            idx = idx + 1
        # self._parameterNode.SetParameter("AnalysisType", self.ui.analysisTypeSelector.checkedIndexes())

        if 0 in selectedIdxs:
            self.ui.additionalQuantificationParameters.collapsed = False
            self.ui.InsertPositionParameters.collapsed = False

            self.ui.firstSliceCTSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            if self.ui.ctSelector.currentNode() is not None:
                if self.ui.firstSliceCTSpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: CT start slice for quantifiction analysis not set.")
                    self.ui.firstSliceCTSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")

            self.ui.firstSlicePETSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            if self.ui.petSelector.currentNode() is not None:
                if self.ui.firstSlicePETSpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: PET start slice for quantifiction analysis not set.")
                    self.ui.firstSlicePETSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")

            self.ui.firstSliceT1SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            if self.ui.t1Selector.currentNode() is not None:
                if self.ui.firstSliceT1SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T1 start slice for quantifiction analysis not set.")
                    self.ui.firstSliceT1SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")

            self.ui.firstSliceT2SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            if self.ui.t2Selector.currentNode() is not None:
                if self.ui.firstSliceT2SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T2 start slice for quantifiction analysis not set.")
                    self.ui.firstSliceT2SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")

        else:
            self.ui.additionalQuantificationParameters.collapsed = True
            self.ui.InsertPositionParameters.collapsed = True

        self.ui.detectResolutionSliceCTButton.setStyleSheet("QPushButton { }")
        self.ui.detectResolutionSlicePETButton.setStyleSheet("QPushButton { }")
        self.ui.detectResolutionSliceT1Button.setStyleSheet("QPushButton { }")
        self.ui.detectResolutionSliceT2Button.setStyleSheet("QPushButton { }")
        if 1 in selectedIdxs:
            self.ui.ResolutionParameters.collapsed = False

            if self.ui.ctSelector.currentNode() is not None:
                nodeName = f"ResolutionTriangle_CT_1"
                markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
                if markupNodeList.GetNumberOfItems() == 0:
                    self.ui.detectResolutionSliceCTButton.setStyleSheet(
                        "QPushButton { border: 2px solid red; background-color: none;}"
                    )
                    self.ui.warningTextBrowser.setText("Warning: CT insert for resolution analysis not set.")

            if self.ui.petSelector.currentNode() is not None:
                nodeName = f"ResolutionTriangle_PET_1"
                markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
                if markupNodeList.GetNumberOfItems() == 0:
                    self.ui.detectResolutionSlicePETButton.setStyleSheet(
                        "QPushButton { border: 2px solid red; background-color: none;}"
                    )
                    self.ui.warningTextBrowser.setText("Warning: PET insert for resolution analysis not set.")

            if self.ui.t1Selector.currentNode() is not None:
                nodeName = f"ResolutionTriangle_T1_1"
                markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
                if markupNodeList.GetNumberOfItems() == 0:
                    self.ui.warningTextBrowser.setText("Warning: T1 insert for resolution analysis not set.")
                    self.ui.detectResolutionSliceT1Button.setStyleSheet(
                        "QPushButton { border: 2px solid red; background-color: none;}"
                    )

            if self.ui.t2Selector.currentNode() is not None:
                nodeName = f"ResolutionTriangle_T2_1"
                markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
                if markupNodeList.GetNumberOfItems() == 0:
                    self.ui.warningTextBrowser.setText("Warning: T2 insert for resolution analysis not set.")
                    self.ui.detectResolutionSliceT2Button.setStyleSheet(
                        "QPushButton { border: 2px solid red; background-color: none;}"
                    )
        else:
            self.ui.ResolutionParameters.collapsed = True

        if 3 in selectedIdxs:
            self.ui.DistortionParameters.collapsed = False

            self.ui.distortionFirstSliceCTSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            self.ui.distortionLastSliceCTSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            self.ui.panalSliceCTSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")

            self.ui.distortionLastSliceT1SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            self.ui.distortionFirstSliceT1SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            self.ui.panalSliceT1SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")

            self.ui.distortionFirstSliceT2SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            self.ui.distortionLastSliceT2SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")
            self.ui.panalSliceT2SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: white; }")

            if self.ui.ctSelector.currentNode() is not None:
                if self.ui.distortionFirstSliceCTSpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: CT slices for distortion analysis not set.")
                    self.ui.distortionFirstSliceCTSpinBox.lineEdit().setStyleSheet(
                        "QLineEdit { background-color: red; }"
                    )

                if self.ui.distortionLastSliceCTSpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: CT slices for distortion analysis not set.")
                    self.ui.distortionLastSliceCTSpinBox.lineEdit().setStyleSheet(
                        "QLineEdit { background-color: red; }"
                    )
                if self.ui.panalSliceCTSpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: CT slices for distortion analysis not set.")
                    self.ui.panalSliceCTSpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")

            if self.ui.t1Selector.currentNode() is not None:
                if self.ui.distortionFirstSliceT1SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T1 slices for distortion analysis not set.")
                    self.ui.distortionFirstSliceT1SpinBox.lineEdit().setStyleSheet(
                        "QLineEdit { background-color: red; }"
                    )

                if self.ui.distortionLastSliceT1SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T1 slices for distortion analysis not set.")
                    self.ui.distortionLastSliceT1SpinBox.lineEdit().setStyleSheet(
                        "QLineEdit { background-color: red; }"
                    )

                if self.ui.panalSliceT1SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T1 slices for distortion analysis not set.")
                    self.ui.panalSliceT1SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")

            if self.ui.t2Selector.currentNode() is not None:
                if self.ui.distortionFirstSliceT2SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T2 slices for distortion analysis not set.")
                    self.ui.distortionFirstSliceT2SpinBox.lineEdit().setStyleSheet(
                        "QLineEdit { background-color: red; }"
                    )

                if self.ui.distortionLastSliceT2SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T2 slices for distortion analysis not set.")
                    self.ui.distortionLastSliceT2SpinBox.lineEdit().setStyleSheet(
                        "QLineEdit { background-color: red; }"
                    )

                if self.ui.panalSliceT2SpinBox.value < 0:
                    self.ui.warningTextBrowser.setText("Warning: T2 slices for distortion analysis not set.")
                    self.ui.panalSliceT2SpinBox.lineEdit().setStyleSheet("QLineEdit { background-color: red; }")
        else:
            self.ui.DistortionParameters.collapsed = True

        if 4 in selectedIdxs:
            self.ui.RadiomicsParameters.collapsed = False
        else:
            self.ui.RadiomicsParameters.collapsed = True

        # Update buttons states and tooltips
        # if (
        #     methodChecked
        #     and self._parameterNode.GetNodeReference("PetInputVolume")
        #     and self._parameterNode.GetNodeReference("CtInputVolume")
        # ) or (
        #     3 in selectedIdxs
        #     and (
        #         self._parameterNode.GetNodeReference("T1InputVolume")
        #         or (self._parameterNode.GetNodeReference("T2InputVolume"))
        #     )
        # ):
        #     self.ui.applyButton.toolTip = "Compute output volume"
        #     self.ui.applyButton.enabled = True
        # else:
        #     self.ui.applyButton.toolTip = "Select input and output volume nodes"
        #     self.ui.applyButton.enabled = False
        #
        if methodChecked:
            self.ui.applyButton.toolTip = "Compute output volume"
            self.ui.applyButton.enabled = True

        # All the GUI updates are done
        self._updatingGUIFromParameterNode = False

    def _checkPropertiesOfSelectedPetNode(self, volumeNode):
        returnValue = False
        if volumeNode and volumeNode.IsA("vtkMRMLScalarVolumeNode"):
            returnValue = True
            modality = self.logic.getDicomTagOrSequence(volumeNode, 0x0008, 0x0060)
            if modality is not None and modality.value == "PT":
                returnValue = False

            acquisitionTimePET = self.logic.getDicomTagOrSequence(
                self._parameterNode.GetNodeReference("PetInputVolume"), 0x0008, 0x0032
            )
            if acquisitionTimePET is not None:
                returnValue = False

        return returnValue

    def _checkFrameOfReferenceUID(self, volumeNode, otherVolumeNode):
        if volumeNode and volumeNode.IsA("vtkMRMLScalarVolumeNode"):
            if otherVolumeNode and otherVolumeNode.IsA("vtkMRMLScalarVolumeNode"):
                frameOfRefUID = self.logic.getDicomTagOrSequence(volumeNode, 0x0020, 0x0052)
                otherFrameOfRefUID = self.logic.getDicomTagOrSequence(otherVolumeNode, 0x0020, 0x0052)
                if frameOfRefUID == otherFrameOfRefUID:
                    return True
        return False

    def _checkPropertiesOfSelectedCtNode(self, volumeNode):
        if volumeNode and volumeNode.IsA("vtkMRMLScalarVolumeNode"):
            modality = self.logic.getDicomTagOrSequence(volumeNode, 0x0008, 0x0060)
            if modality is not None and modality.value == "CT":
                return True
        return False

    def _checkPropertiesOfSelectedMRNode(self, volumeNode):
        if volumeNode and volumeNode.IsA("vtkMRMLScalarVolumeNode"):
            modality = self.logic.getDicomTagOrSequence(volumeNode, 0x0008, 0x0060)
            if modality is not None and modality.value == "MR":
                return True
        return False

    def timerEvent(self, event=None):
        self.timer = None
        self.updateParameterNodeFromGUI()

    def updateParameterNodeFromGUIWithTimer(self, caller=None, event=None):
        if self.timer is not None:
            # self.timer.killTimer(self.timer.timerId())
            self.timer.stop()
        self.timer = qt.QTimer()
        self.timer.timeout.connect(self.timerEvent)
        self.timer.start(2000)

    def updateParameterNodeFromGUI(self, caller=None, event=None):
        """
        This method is called when the user makes any change in the GUI.
        The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        wasModified = self._parameterNode.StartModify()
        if self._parameterNode.GetNodeReferenceID("PetInputVolume") != self.ui.petSelector.currentNodeID:
            self._updatePetGuiValues = True
        else:
            self._updatePetGuiValues = False
        self._parameterNode.SetNodeReferenceID("PetInputVolume", self.ui.petSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("CtInputVolume", self.ui.ctSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("T1InputVolume", self.ui.t1Selector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("T2InputVolume", self.ui.t2Selector.currentNodeID)

        self._parameterNode.SetParameter("decayFactor", str(self.ui.decayFactor.value))
        self._parameterNode.SetParameter("rcVolume", str(self.ui.rcVolume.value))

        idx = 0
        for analysesMethod in self.logic.getSupportedAnalysesMethods().keys():
            if self.ui.analysisTypeSelector.model().item(idx).checkState() == 2:
                self._parameterNode.SetParameter(analysesMethod, "1")
            else:
                self._parameterNode.SetParameter(analysesMethod, "0")
            idx = idx + 1

        if self._parameterNode.GetParameter("insertTypeSelector") != str(self.ui.insertTypeSelector.currentText):
            self._updateInsertTypeValues = True
        else:
            self._updateInsertTypeValues = False
        self._parameterNode.SetParameter("insertTypeSelector", str(self.ui.insertTypeSelector.currentText))

        self._parameterNode.SetParameter("tracerSelector", str(self.ui.tracerSelector.currentText))

        self._parameterNode.SetParameter("initialActivity", str(self.ui.initialActivity.value))

        self._parameterNode.SetParameter("timeActivity", str(self.ui.timeActivity.value))

        self._parameterNode.SetParameter("firstSliceCT", str(self.ui.firstSliceCTSpinBox.value))
        self._parameterNode.SetParameter("firstSlicePET", str(self.ui.firstSlicePETSpinBox.value))
        self._parameterNode.SetParameter("firstSliceT1", str(self.ui.firstSliceT1SpinBox.value))
        self._parameterNode.SetParameter("firstSliceT2", str(self.ui.firstSliceT2SpinBox.value))

        self._parameterNode.SetParameter("distortionFirstSliceCT", str(self.ui.distortionFirstSliceCTSpinBox.value))
        self._parameterNode.SetParameter("distortionLastSliceCT", str(self.ui.distortionLastSliceCTSpinBox.value))
        self._parameterNode.SetParameter("panalSliceCT", str(self.ui.panalSliceCTSpinBox.value))
        self._parameterNode.SetParameter("panalSliceT1", str(self.ui.panalSliceT1SpinBox.value))
        self._parameterNode.SetParameter("panalSliceT2", str(self.ui.panalSliceT2SpinBox.value))
        self._parameterNode.SetParameter("distortionFirstSliceT1", str(self.ui.distortionFirstSliceT1SpinBox.value))
        self._parameterNode.SetParameter("distortionLastSliceT1", str(self.ui.distortionLastSliceT1SpinBox.value))
        self._parameterNode.SetParameter("distortionFirstSliceT2", str(self.ui.distortionFirstSliceT2SpinBox.value))
        self._parameterNode.SetParameter("distortionLastSliceT2", str(self.ui.distortionLastSliceT2SpinBox.value))

        for insertPosKey in self.logic.materialCharacteristicsPosition:
            uiElement = getattr(self.ui, "p" + str(insertPosKey))
            self._parameterNode.SetParameter("p" + str(insertPosKey), str(uiElement.value))
            uiElement = getattr(self.ui, "l" + str(insertPosKey))
            self._parameterNode.SetParameter("l" + str(insertPosKey), str(uiElement.text))

        self._parameterNode.EndModify(wasModified)

    def onApplyButton(self):
        """
        Run processing when user clicks "Apply" button.
        """
        try:
            checkedIndics = self.ui.analysisTypeSelector.checkedIndexes()
            checkedMethods = [checkedIdx.data() for checkedIdx in checkedIndics]

            instertPositions = {}
            for insertPosKey in self.logic.materialCharacteristicsPosition:
                uiElement = getattr(self.ui, "p" + str(insertPosKey))
                position = uiElement.value
                uiElement = getattr(self.ui, "l" + str(insertPosKey))
                insertName = str(uiElement.text)
                instertPositions[position] = insertName

            self.logic.setOutputDirectory(self.ui.outputPahtDirectoryButton.directory)
            self.logic.setShowPlots(self.ui.showPlotsCheckBox.checked)

            # Compute output
            self.logic.process(
                self.ui.petSelector.currentNode(),
                self.ui.ctSelector.currentNode(),
                self.ui.t1Selector.currentNode(),
                self.ui.t2Selector.currentNode(),
                self.ui.petSegmentSelector.currentNode(),
                self.ui.ctSegmentSelector.currentNode(),
                self.ui.t1SegmentSelector.currentNode(),
                self.ui.t2SegmentSelector.currentNode(),
                self.ui.rcVolume.value,
                self.ui.decayFactor.value,
                self.ui.insertTypeSelector.currentText,
                self.ui.tracerSelector.currentText,
                self.ui.initialActivity.value,
                self.ui.timeActivity.value,
                instertPositions,
                checkedMethods,
                self.ui.defaultRadiomicSegmentPathButton.directory,
                self.ui.firstSlicePETSpinBox.value,
                self.ui.firstSliceCTSpinBox.value,
                self.ui.firstSliceT1SpinBox.value,
                self.ui.firstSliceT2SpinBox.value,
                self.ui.distortionFirstSliceCTSpinBox.value,
                self.ui.distortionLastSliceCTSpinBox.value,
                self.ui.panalSliceCTSpinBox.value,
                self.ui.panalSliceT1SpinBox.value,
                self.ui.panalSliceT2SpinBox.value,
                self.ui.distortionFirstSliceT1SpinBox.value,
                self.ui.distortionLastSliceT1SpinBox.value,
                self.ui.distortionFirstSliceT2SpinBox.value,
                self.ui.distortionLastSliceT2SpinBox.value,
                self.ui.distortionInsertFilledCheckBox.checked,
            )

            # Thread_obj = Thread(target=self.logic.process, args=(self.ui.petSelector.currentNode(), self.ui.ctSelector.currentNode(),checkedMethods))
            # Thread_obj.start()
            # # Compute inverted output (if needed)
            # if self.ui.invertedOutputSelector.currentNode():
            #   # If additional output volume is selected then result with inverted threshold is written there
            #   self.logic.process(self.ui.petSelector.currentNode(), self.ui.invertedOutputSelector.currentNode(),
            #     self.ui.imageThresholdSliderWidget.value, not self.ui.invertOutputCheckBox.checked, showResult=False)

        except Exception as e:
            slicer.util.errorDisplay("Failed to compute results: " + str(e))
            import traceback

            traceback.print_exc()


#
# QAHybridLogic
#


class QAHybridLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)

        self.petVolume = None
        self.ctVolume = None
        self.t1Volume = None
        self.t2Volume = None

        self.petSegment = None
        self.ctSegment = None
        self.t1Segment = None
        self.t2Segment = None

        self.resourcePath = None

        self.cylinderNames = None
        self.updatedPoints = set()

        self.decayFactor = 6586.200165
        self.rcVolume = 1200.0
        self.insertType = "Calibration Curve"
        self.tracerType = "18F-FDG"
        self.timeActivity = 143500
        self.initialActivity = 69190000
        self.insertPositions = None

        self.quantificationStartSlicePET = -1
        self.quantificationStartSliceCT = -1
        self.quantificationStartSliceT1 = -1
        self.quantificationStartSliceT2 = -1

        self.outputDirectory = "."
        self.showOutputPlots = False

        self.theor_diameter_pet = ["D10", "D28a", "D28b", "D28c", "D28d", "D20", "D12"]

        self.callibrationCurvePositions = {
            1: "D10",
            2: "D28a",
            3: "D28b",
            4: "D28c",
            5: "D28d",
            6: "D20",
            7: "D12",
            8: "Trabecular bone",
            9: "Dense bone",
            10: "Muscle",
            11: "Lung (exhale)",
        }

        self.materialCharacteristicsPosition = {
            1: "D10",
            2: "D28a",
            3: "D28b",
            4: "D28c",
            5: "D28d",
            6: "D20",
            7: "D12",
            8: "DAP0102",
            9: "SolidWater",
            10: "HW04",
            11: "DAP0203",
        }
        self.insertDensities = {
            "Muscle": 1.06,
            "Lung (inhale)": 0.20,
            "Liver": 1.07,
            "Breast": 0.99,
            "Water": 1.004,
            "Adipose": 0.97,
            "Trabecular bone": 1.16,
            "Lung (exhale)": 0.50,
            "Dense bone": 1.61,
        }

        self.pltPlots = {}

    def setShowPlots(self, showOutputPlots):
        self.showOutputPlots = showOutputPlots

    def setOutputDirectory(self, outputDirectory):
        self.outputDirectory = outputDirectory

    def setDefaultParameters(self, parameterNode):
        """
        Initialize parameter node with default settings.
        """
        if not parameterNode.GetParameter("rcVolume"):
            parameterNode.SetParameter("rcVolume", str(self.rcVolume))
        if not parameterNode.GetParameter("decayFactor"):
            parameterNode.SetParameter("decayFactor", str(self.decayFactor))
        if not parameterNode.GetParameter("timeActivity"):
            parameterNode.SetParameter("timeActivity", str(self.timeActivity))
        if not parameterNode.GetParameter("initialActivity"):
            parameterNode.SetParameter("initialActivity", str(self.initialActivity))
        if not parameterNode.GetParameter("vv"):
            parameterNode.SetParameter("insertType", self.insertType)
        if not parameterNode.GetParameter("tracerType"):
            parameterNode.SetParameter("tracerType", self.tracerType)

        self.setInsertTypePositions(self.insertType, parameterNode)

    def setInsertTypePositions(self, insertType, parameterNode):
        if insertType == "Calibration Curve":
            insertPositions = self.callibrationCurvePositions
        else:
            insertPositions = self.materialCharacteristicsPosition
        for insertPosKey in insertPositions.keys():
            parameterNode.SetParameter("p" + str(insertPosKey), str(insertPosKey))
            parameterNode.SetParameter("l" + str(insertPosKey), insertPositions[insertPosKey])

    def getDicomTag(self, slicerNode, dicomTag):
        shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
        structSeriesUID = shNode.GetItemUID(shNode.GetItemByDataNode(slicerNode), "DICOM")
        instUids = slicer.dicomDatabase.instancesForSeries(structSeriesUID)
        return slicer.dicomDatabase.instanceValue(instUids[0], dicomTag)

    # #
    # # if we want to access sequences we can use the arguemnts in *args
    # # the first one should be the element number of the sequence; then another group and elemt value can follow
    # # e.g. group = 0054, element= 0016, *args = 0, 0054, 0300, 0, 0008, 0100
    # #
    # #
    def getDicomTagOrSequence(self, slicerNode, group, element, *args):
        dt = None
        if slicerNode is not None:
            shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
            structSeriesUID = shNode.GetItemUID(shNode.GetItemByDataNode(slicerNode), "DICOM")
            fileList = slicer.dicomDatabase.filesForSeries(structSeriesUID)
            dt = self.getDicomTagOrSequenceForFile(fileList[0], group, element, *args)
        return dt

    def getAllDicomTagsOrSequences(self, slicerNode, group, element, *args):
        dts = []
        if slicerNode is not None:
            shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
            structSeriesUID = shNode.GetItemUID(shNode.GetItemByDataNode(slicerNode), "DICOM")
            fileList = slicer.dicomDatabase.filesForSeries(structSeriesUID)
            for fileName in fileList:
                dt = self.getDicomTagOrSequenceForFile(fileName, group, element, *args)
                dts.append(dt)
        return dts

    def getDicomTagOrSequenceForFile(self, file, group, element, *args):
        dt = None
        if file is not None:
            ds = pydicom.dcmread(file)
            if (group, element) in ds:
                dt = ds[group, element]
                nextItemIsSequence = True
                i = 0
                while dt is not None:
                    if nextItemIsSequence and i < len(args):
                        if args[i] < len(dt.value):
                            dt = dt[args[i]]
                            i = i + 1
                        else:
                            dt = None
                    elif not nextItemIsSequence and i + 1 < len(args):
                        if (args[i], args[i + 1]) in dt:
                            dt = dt[args[i], args[i + 1]]
                            i = i + 2
                        else:
                            dt = None
                    else:
                        break
                    nextItemIsSequence = not nextItemIsSequence
        return dt

    def getSupportedAnalysesMethods(self):
        analyisMethods = {
            "Quantification": self.doQuantification,
            "Resolution": self.doResolution,
            "Registration": self.doRegistration,
            "Distortion": self.doDistortion,
            "Radiomics": self.doRadiomics,
        }
        return analyisMethods

    def getSupportedInserts(self):
        supportedInserts = {
            "Calibration Curve": self.insertTypeCalibrationCurve,
            "Materials Characterization": self.insertTypeMaterialsCharacterization,
        }
        return supportedInserts

    def getSupportedTracers(self):
        return ("18F-FDG", "other")

    def getMarkupNodeList(self, markupNodeListName, parentNodeID=None):
        markupNodeList = slicer.mrmlScene.GetNodesByName(markupNodeListName)
        if markupNodeList.GetNumberOfItems() > 0:
            pointListNode = markupNodeList.GetItemAsObject(0)
        else:
            pointListNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode", markupNodeListName)

            if parentNodeID:
                shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
                shNode.CreateItem(parentNodeID, pointListNode)
        return pointListNode

    def getTableNode(self, tableName, parentNodeID=None):
        tableNodeList = slicer.mrmlScene.GetNodesByName(tableName)
        if tableNodeList.GetNumberOfItems() > 0:
            tableNode = tableNodeList.GetItemAsObject(0)
        else:
            tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", tableName)
            if parentNodeID:
                shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
                shNode.CreateItem(parentNodeID, tableNode)
        return tableNode

    def getQuantificationTableNode(self, modality=None, parentNodeID=None):
        return self.getTableNode("Quantification" + modality + "_" + str(parentNodeID), parentNodeID)

    def addSitkVolume(self, sitkVolume, volumeName, volumeType, parentNodeID):
        nodeList = slicer.mrmlScene.GetNodesByName(volumeName)
        if nodeList.GetNumberOfItems() > 0:
            volumeNode = nodeList.GetItemAsObject(0)
        else:
            volumeNode = slicer.mrmlScene.AddNewNodeByClass(volumeType, volumeName)
            if parentNodeID:
                shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
                shNode.CreateItem(parentNodeID, volumeNode)
        sitkUtils.PushVolumeToSlicer(sitkVolume, volumeNode)
        return volumeNode

    def getSegmentationNode(self, nodeName, parentNodeID=None):
        nodeList = slicer.mrmlScene.GetNodesByName(nodeName)
        if nodeList.GetNumberOfItems() > 0:
            node = nodeList.GetItemAsObject(0)
            slicer.mrmlScene.RemoveNode(node)
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", nodeName)
        if parentNodeID:
            shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
            shNode.CreateItem(parentNodeID, node)
        return node

    def addRadiomicSegmentations(self, radiomicSegmentation, modality, parentNodeID):
        segmentationNodeName = "Radiomics" + modality + "Segmentation" + str(parentNodeID)
        segmentationNode = self.getSegmentationNode(segmentationNodeName, parentNodeID)

        for contourName in radiomicSegmentation.keys():
            segmentedVolume = radiomicSegmentation[contourName]
            tempVolumeNode = self.addSitkVolume(segmentedVolume, contourName, "vtkMRMLLabelMapVolumeNode", parentNodeID)
            slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(tempVolumeNode, segmentationNode)
            slicer.mrmlScene.RemoveNode(tempVolumeNode)

    def getRadiomicsTableNode(self, modality, parentNodeID=None):
        return self.getTableNode("Radiomics" + modality + "Table" + str(parentNodeID), parentNodeID)

    def createRadiomicsTable(self, radiomicFeatures, modality, parentNodeID):
        nameCol = vtk.vtkStringArray()
        nameCol.SetName("Segmentation")

        rfFeatureColumns = {}

        radiomicsTable = self.getRadiomicsTableNode(modality, parentNodeID)
        while radiomicsTable.RemoveRow(0):
            pass

        for contourName in radiomicFeatures.keys():
            nameCol.InsertNextValue(contourName)
            radiomicFeaturesForContour = radiomicFeatures[contourName]
            for radiomicFeatureName in radiomicFeaturesForContour.keys():
                if radiomicFeatureName in rfFeatureColumns.keys():
                    rfCol = rfFeatureColumns[radiomicFeatureName]
                else:
                    rfCol = vtk.vtkDoubleArray()
                    rfCol.SetName(radiomicFeatureName)
                    rfFeatureColumns[radiomicFeatureName] = rfCol

                rfCol.InsertNextValue(radiomicFeaturesForContour[radiomicFeatureName])

        radiomicsTable.RemoveAllColumns()
        radiomicsTable.AddColumn(nameCol)

        for radiomicFeatureName in rfFeatureColumns.keys():
            radiomicsTable.AddColumn(rfFeatureColumns[radiomicFeatureName])

    def calculateActivityDecay(self):
        activity_decay = 0.0
        if self.petVolume is not None:
            acquisitionTimesPet = self.getAllDicomTagsOrSequences(self.petVolume, 0x0008, 0x0032)
            minAcquisitionTime = None
            for acTime in acquisitionTimesPet:
                acTime = acTime.value.split(".")[0]
                if minAcquisitionTime is None or minAcquisitionTime > acTime:
                    minAcquisitionTime = acTime
            acquisitionTimePET = af.time_converter(minAcquisitionTime)

            timeActivity = af.time_converter(str(self.timeActivity))

            radionuclideTotalDose = self.initialActivity

            diff_time = abs(acquisitionTimePET - timeActivity)
            activity_decay = (radionuclideTotalDose / self.rcVolume) * math.exp(
                -diff_time * ((math.log(2)) / self.decayFactor)
            )

        return activity_decay

    def createQuantificationTablePET(
        self,
        parentNodeID,
        names,
        diameters,
        meansPT,
        devsPT,
        rcs,
    ):
        nameCol = vtk.vtkStringArray()
        nameCol.SetName("Position")
        insertCol = vtk.vtkStringArray()
        insertCol.SetName("Insert")
        diameterCol = vtk.vtkDoubleArray()
        diameterCol.SetName("Diameter")
        meanCol = vtk.vtkDoubleArray()
        meanCol.SetName("Mean( Bq/m)")
        devCol = vtk.vtkDoubleArray()
        devCol.SetName("Standard Deviation (Bq/ml)")
        rcCol = vtk.vtkDoubleArray()
        rcCol.SetName("RC_Activity")

        insertNames = self.theor_diameter_pet

        diameters, insertNames, names, meansPT, devsPT, rcs = self.sort_lists(
            diameters, insertNames, names, meansPT, devsPT, rcs, descendingOrder=False
        )

        qunatificationTableNode = self.getQuantificationTableNode("PET", parentNodeID)
        while qunatificationTableNode.RemoveRow(0):
            pass

        for i in range(0, len(names)):
            nameCol.InsertNextValue(str(names[i]))
            insertCol.InsertNextValue(str(insertNames[i]))
            meanCol.InsertNextValue(meansPT[i])
            devCol.InsertNextValue(devsPT[i])
            rcCol.InsertNextValue(rcs[i])
            diameterCol.InsertNextValue(diameters[i])

        qunatificationTableNode.RemoveAllColumns()
        qunatificationTableNode.AddColumn(nameCol)
        qunatificationTableNode.AddColumn(insertCol)
        qunatificationTableNode.AddColumn(diameterCol)
        qunatificationTableNode.AddColumn(meanCol)
        qunatificationTableNode.AddColumn(devCol)
        qunatificationTableNode.AddColumn(rcCol)

    def sort_lists(self, reference, *lists, descendingOrder: bool):
        zipped = list(zip(reference, *lists))
        sorted_zipped = sorted(zipped, reverse=descendingOrder)
        return [list(t) for t in zip(*sorted_zipped)]

    def createQuantificationTableMR(self, modality, parentNodeID, names, meansCT, devsCT):
        nameCol = vtk.vtkStringArray()
        nameCol.SetName("Position")

        volumeColCt = vtk.vtkStringArray()
        volumeColCt.SetName("Inserts")
        meanColCt = vtk.vtkDoubleArray()
        meanColCt.SetName("Mean value")
        devColCt = vtk.vtkDoubleArray()
        devColCt.SetName("Mean dev.")

        diameterCT = list(self.insertPositions.values())

        qunatificationTableNode = self.getQuantificationTableNode("MR" + modality, parentNodeID)
        while qunatificationTableNode.RemoveRow(0):
            pass

        for i in range(0, len(names)):
            nameCol.InsertNextValue(str(names[i]))
            meanColCt.InsertNextValue(meansCT[i])
            devColCt.InsertNextValue(devsCT[i])
            volumeColCt.InsertNextValue(diameterCT[i])

        qunatificationTableNode.RemoveAllColumns()
        qunatificationTableNode.AddColumn(nameCol)
        qunatificationTableNode.AddColumn(volumeColCt)
        qunatificationTableNode.AddColumn(meanColCt)
        qunatificationTableNode.AddColumn(devColCt)

    def createQuantificationTableCT(self, parentNodeID, names, meansCT, devsCT):
        nameCol = vtk.vtkStringArray()
        nameCol.SetName("Position")

        volumeColCt = vtk.vtkStringArray()
        volumeColCt.SetName("Inserts")
        meanColCt = vtk.vtkDoubleArray()
        meanColCt.SetName("Mean value (HU)")
        devColCt = vtk.vtkDoubleArray()
        devColCt.SetName("Mean dev. (HU)")
        densityColCt = vtk.vtkDoubleArray()
        densityColCt.SetName("Density (g/cm3")

        densityValues = []
        for idx, _ in enumerate(meansCT):
            idxP1 = idx + 1
            densityFound = False
            if idxP1 in self.insertPositions.keys():
                insert_name = self.insertPositions[idxP1]
                if insert_name in self.insertDensities:
                    densityValues.append(self.insertDensities[insert_name])
                    densityFound = True
            if not densityFound:
                densityValues.append(-1.0)

        diameterCT = list(self.insertPositions.values())

        densityValues, names, meansCT, devsCT, diameterCT = self.sort_lists(
            densityValues, names, meansCT, devsCT, diameterCT, descendingOrder=False
        )

        qunatificationTableNode = self.getQuantificationTableNode("CT", parentNodeID)
        while qunatificationTableNode.RemoveRow(0):
            pass

        for i in range(0, len(names)):
            nameCol.InsertNextValue(str(names[i]))
            meanColCt.InsertNextValue(meansCT[i])
            devColCt.InsertNextValue(devsCT[i])
            volumeColCt.InsertNextValue(diameterCT[i])
            densityColCt.InsertNextValue(densityValues[i])

        qunatificationTableNode.RemoveAllColumns()
        qunatificationTableNode.AddColumn(nameCol)
        qunatificationTableNode.AddColumn(volumeColCt)
        qunatificationTableNode.AddColumn(meanColCt)
        qunatificationTableNode.AddColumn(devColCt)
        qunatificationTableNode.AddColumn(densityColCt)

    def registrationPlot(self, diceCoefficients, DIF_CORTES, modality):
        diceVals = list(diceCoefficients.values())
        zCoordinates = list(diceCoefficients.keys())
        mean_dice = np.mean(diceVals)
        std_dice = np.std(diceVals)
        min_dice = np.min(diceVals)
        max_dice = np.max(diceVals)

        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        ax1.plot(zCoordinates, diceVals, color=cls.to_rgba("C7", 0.7))
        ax1.axhline(y=mean_dice, color=cls.to_rgba("C9", 0.8))
        ax1.axhline(y=0.93, color=cls.to_rgba("C10", 0.8), linestyle="-")
        ax1.fill_between(zCoordinates, diceVals, color=cls.to_rgba("C7", 0.3))
        ax1.axhline(y=1, color="k", linestyle="dashed")
        ax1.set_xlim([min(zCoordinates), max(zCoordinates)])
        ax1.set_ylim([0, 1.05])
        ax1.set_xlabel("Slice", fontsize=15)
        ax1.set_ylabel("DICE", fontsize=15)

        if mean_dice >= 0.91:
            if DIF_CORTES > 2:
                color = ("black", "green", "black")
                ax1.legend(
                    [
                        f"DSC coefficient: The {modality} system is not properly co-registered",
                        f"Mean DSC= {round(mean_dice, 2)} \u00B1 {round(std_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm",
                        f"DSC threshold value: 0.91",
                    ],
                    labelcolor=color,
                    title=f"Co-registration properties",
                )
            else:
                color = ("black", "green", "black")
                ax1.legend(
                    [
                        f"DSC coefficient: The {modality} system is well co-registered",
                        f"Mean DSC= {round(mean_dice, 2)} \u00B1 {round(std_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm",
                        f"DSC threshold value: 0.91",
                    ],
                    labelcolor=color,
                    title=f"Co-registration properties",
                )
        else:
            color = ("black", "red", "black")
            ax1.legend(
                [
                    f"DSC coefficient: The {modality} system is not properly co-registered",
                    f"Mean DSC= {round(mean_dice, 2)} \u00B1 {round(std_dice, 2)} and $\Delta$z={round(DIF_CORTES, 2)} mm",
                    f"DSC threshold value: 0.91",
                ],
                labelcolor=color,
                title=f"Co-registration properties",
            )

        loc_min = []
        loc_max = []
        for i in range(0, len(diceVals)):
            if diceVals[i] == min_dice:
                loc_min.append(i)
            if diceVals[i] == max_dice:
                loc_max.append(i)

        if len(loc_min) == 1:
            ax1.plot(loc_min[0] + min(zCoordinates), min_dice, marker=".", color="red")
            ax1.text(loc_min[0] + min(zCoordinates), min_dice - 0.03 * min_dice, round(min_dice, 4))

        if len(loc_min) != 1 and len(loc_min) != 0:
            ax1.plot(loc_min[0] + min(zCoordinates), min_dice, marker=".", color="red")
            ax1.text(loc_min[0] + min(zCoordinates), min_dice - 0.03 * min_dice, round(min_dice, 4))
            for i in range(1, len(loc_min)):
                ax1.plot(loc_min[i] + min(zCoordinates), min_dice, marker=".", color="red")

        if len(loc_max) == 1:
            ax1.plot(loc_max[0] + min(zCoordinates), max_dice, marker=".", color="green")
            ax1.text(loc_max[0] + min(zCoordinates), max_dice - 0.03 * max_dice, round(max_dice, 4))

        if len(loc_max) != 1 and len(loc_max) != 0:
            ax1.plot(loc_max[0] + min(zCoordinates), max_dice, marker=".", color="green")
            ax1.text(loc_max[0] + min(zCoordinates), max_dice - 0.03 * max_dice, round(max_dice, 4))
            for i in range(1, len(loc_max)):
                ax1.plot(loc_max[i] + min(zCoordinates), max_dice, marker=".", color="green")

        filename = os.path.join(self.outputDirectory, f"Coregistration_DSC_PET{modality}.png")
        plt.savefig(filename, bbox_inches="tight")
        if self.showOutputPlots:
            pm = qt.QPixmap(filename)
            imageWidget = qt.QLabel()
            imageWidget.setPixmap(pm)
            imageWidget.setScaledContents(True)
            imageWidget.show()
            self.pltPlots[f"Coregistration_DSC_PET{modality}"] = imageWidget
        else:
            plt.close()

    def resolutionPlot2(self, resultValuesForPlt2, modality):
        valid_pixels, x_fit, y_fit, centroid1, centroid2, th_cyl = resultValuesForPlt2
        plt.figure(figsize=(10, 6))
        plt.hist(valid_pixels, bins=100, density=True, alpha=0.6, color="g", label="Datos")
        plt.plot(x_fit, y_fit, label="Ajuste Doble Gaussiana", color="blue")
        plt.axvline(centroid1, color="red", linestyle="dashed", linewidth=1, label=f"Centroide Agua: {centroid1:.2f}")
        plt.axvline(
            centroid2, color="orange", linestyle="dashed", linewidth=1, label=f"Centroide Cilindros: {centroid2:.2f}"
        )
        plt.axvline(
            th_cyl, color="orange", linestyle="dotted", linewidth=1, label=f"Límite Inferior Cilindros: {th_cyl:.2f}"
        )
        plt.legend()
        plt.xlabel("Intensidad de píxel")
        plt.ylabel("Densidad")
        plt.title("Ajuste de Doble Gaussiana")
        # plt.show()

        filename = os.path.join(self.outputDirectory, f"Resolution{modality}GaussianPlot.png")
        plt.savefig(filename, bbox_inches="tight")
        if self.showOutputPlots:
            pm = qt.QPixmap(filename)
            imageWidget = qt.QLabel()
            imageWidget.setPixmap(pm)
            imageWidget.setScaledContents(True)
            imageWidget.show()
            self.pltPlots[f"Resolution{modality}GaussianPlot"] = imageWidget
        else:
            plt.close()

    def resolutionPlot(
        self,
        plotImageValuesList,
        normalized_profile,
        starting_coord,
        plotThreshHolesValuesList,
        new_end_coord,
        unit,
        modality,
    ):
        list_diameters = [5, 7.5, 9, 11, 12, 15]
        for section in range(0, len(plotImageValuesList)):
            fig, ax = plt.subplots(nrows=1, ncols=2, num=f"RES{modality}SEC{section}", clear=True)
            fig.suptitle(f"{modality} Section {section+1}")
            ax[0].plot(normalized_profile[section])
            ax[0].set_title(f"Resolution {modality} Section {section} (D={list_diameters[section]}mm)")
            ax[0].set_ylabel("Normalized Intensity")

            ax[1].imshow(plotImageValuesList[section], cmap="gray")
            ax[1].imshow(plotThreshHolesValuesList[section], cmap="Blues", alpha=0.3)
            # ax[1].arrow(
            #     starting_coord[section][2],  # x start
            #     starting_coord[section][1],  # y start
            #     new_end_coord[section][2] - starting_coord[section][2],  # dx
            #     new_end_coord[section][1] - starting_coord[section][1],  # dy
            #     color="green",
            #     width=1.0,  # 0.001,  # arrow width
            #     head_width=0.2,  # arrow head width
            #     head_length=0.3,  # arrow head length
            #     length_includes_head=True,  # include head in arrow length
            # )
            ax[1].annotate(
                "",
                xy=(new_end_coord[section][2], new_end_coord[section][1]),  # arrow tip
                xytext=(starting_coord[section][2], starting_coord[section][1]),  # arrow base
                arrowprops=dict(
                    color="red",
                    width=1,
                    headwidth=5,  # increase this for a bigger head
                    headlength=5,  # increase this for a longer head
                ),
            )
            ax[1].grid(False)
            ax[1].set_title("Intensity profile")

            filename = os.path.join(self.outputDirectory, f"Resolution{modality}Plot{section}.png")
            plt.savefig(filename, bbox_inches="tight")
            if self.showOutputPlots:
                pm = qt.QPixmap(filename)
                imageWidget = qt.QLabel()
                imageWidget.setPixmap(pm)
                imageWidget.setScaledContents(True)
                imageWidget.show()
                self.pltPlots[f"Resolution{modality}Plot{section}"] = imageWidget
            else:
                plt.close(fig)

    def quantificationPlotPET(self, rc_ac):
        data_points = list(zip(self.theor_diameter_pet, rc_ac))
        data_points_sorted = sorted(data_points, key=lambda point: point[1])
        theor_diameter_pet_sorted, rc_ac_sorted = zip(*data_points_sorted)
        _ = plt.figure(num="RC", clear=True)
        plt.plot(theor_diameter_pet_sorted, rc_ac_sorted, "o", label="Activity concentration RC")
        plt.xlabel("Insert diameter (mm)")
        plt.ylabel("Recovery Coefficient")
        plt.ylim([0.3, 1.3])
        plt.title("Quantification PET")
        plt.legend()

        filename = os.path.join(self.outputDirectory, "Quantification_PET_RC.png")
        plt.savefig(filename, bbox_inches="tight")
        if self.showOutputPlots:
            pm = qt.QPixmap(filename)
            imageWidget = qt.QLabel()
            imageWidget.setPixmap(pm)
            imageWidget.setScaledContents(True)
            imageWidget.show()
            self.pltPlots["Quantification_PET_RC"] = imageWidget
        else:
            plt.close()

    def quantificationPlotCT(self, intensity_ct):
        insert_position = []
        density_list = []
        for insertPosKey in self.insertPositions.keys():
            insert_name = self.insertPositions[insertPosKey]
            if insert_name in self.insertDensities:
                insert_position.append(insertPosKey)
                density_list.append(self.insertDensities[insert_name])

        if len(insert_position) > 0:
            intensity_ct_list = np.array([intensity_ct[indice - 1] for indice in insert_position])

            regression = np.polyfit(density_list, intensity_ct_list, 1)
            x = np.linspace(np.min(density_list), np.max(density_list), 2000)
            y = regression[0] * x + regression[1]

            slope, intercept, r_value, _, _ = stats.linregress(density_list, intensity_ct_list)
            R2 = r_value**2

            plt.figure(num="HUDENS", clear=True)
            plt.plot(density_list, intensity_ct_list, "bo")
            plt.plot(x, y, "g--", alpha=0.5, label="Linear regression")
            if intercept > 0:
                plt.title(f"Linear regression: {slope:.0f}x + {intercept:.0f}. R square value = {R2:.4f}")
            else:
                plt.title(f"Linear regression: {slope:.0f}x - {-intercept:.0f}. R square value = {R2:.4f}")
            plt.xlabel("Density (g/cm³)")
            plt.ylabel("Hounsfield Units")

            filename = os.path.join(self.outputDirectory, "Quantification_CT_CalibrationCurve.png")
            if os.path.exists(filename):
                os.remove(filename)

            plt.savefig(filename, bbox_inches="tight")
            if self.showOutputPlots:
                pm = qt.QPixmap(filename)
                imageWidget = qt.QLabel()
                imageWidget.setPixmap(pm)
                imageWidget.setScaledContents(True)
                imageWidget.show()
                self.pltPlots["Quantification_CT_CalibrationCurve"] = imageWidget
            else:
                plt.close()

    def radiomicsCallback(self, radiomicValues, radiomicSegmentations):
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        for modality in radiomicValues.keys():
            radiomics = radiomicValues[modality]
            radiomicSegmentation = radiomicSegmentations[modality]
            if modality == "PET":
                parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.petVolume))
            elif modality == "CT":
                parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.ctVolume))
            elif modality == "T1":
                parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.t1Volume))
            elif modality == "T2":
                parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.t2Volume))

            self.createRadiomicsTable(radiomics, modality, parentNodeID)
            if len(radiomicSegmentation) > 0:
                self.addRadiomicSegmentations(radiomicSegmentation, modality, parentNodeID)

    def quantificationCallback(
        self,
        names_ct,
        intensity_ct,
        err_intensity_ct,
        names_pet,
        diameter_pet,
        intensity_pet,
        err_intensity_pet,
        rc_ac,
        names_t1,
        intensity_t1,
        err_intensity_t1,
        names_t2,
        intensity_t2,
        err_intensity_t2,
    ):
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()

        if self.ctVolume is not None and names_ct is not None:
            parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.ctVolume))
            self.createQuantificationTableCT(parentNodeID, names_ct, intensity_ct, err_intensity_ct)
            self.quantificationPlotCT(intensity_ct)

        if self.petVolume is not None and names_pet is not None:
            parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.petVolume))
            self.createQuantificationTablePET(
                parentNodeID,
                names_pet,
                diameter_pet,
                intensity_pet,
                err_intensity_pet,
                rc_ac,
            )
            self.quantificationPlotPET(rc_ac)

        if self.t1Volume is not None and names_t1 is not None:
            parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.t1Volume))
            self.createQuantificationTableMR("T1", parentNodeID, names_t1, intensity_t1, err_intensity_t1)

        if self.t2Volume is not None and names_t2 is not None:
            parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(self.t2Volume))
            self.createQuantificationTableMR("T2", parentNodeID, names_t2, intensity_t2, err_intensity_t2)

    def getResolutionTableNode(self, modality=None, parentNodeID=None):
        return self.getTableNode("Resolution" + modality + "_" + str(parentNodeID), parentNodeID)

    def createDistortionTablePanalWalls(
        self,
        parentNodeID,
        modality,
        distances_coord_r_mm,
        distances_coord_l_mm,
        distances_coord_dsu_mm,
        distances_coord_dsd_mm,
        distances_coord_dmu_mm,
        distances_coord_dmd_mm,
    ):
        rightCol = vtk.util.numpy_support.numpy_to_vtk(distances_coord_r_mm)
        leftCol = vtk.util.numpy_support.numpy_to_vtk(distances_coord_l_mm)
        upRightCol = vtk.util.numpy_support.numpy_to_vtk(distances_coord_dsu_mm)
        downLeftCol = vtk.util.numpy_support.numpy_to_vtk(distances_coord_dsd_mm)
        upLeftCol = vtk.util.numpy_support.numpy_to_vtk(distances_coord_dmu_mm)
        downRightCol = vtk.util.numpy_support.numpy_to_vtk(distances_coord_dmd_mm)

        rightCol.SetName("Right (mm)")
        leftCol.SetName("Left (mm)")
        upRightCol.SetName("Diag. up-right (mm)")
        downLeftCol.SetName("Diag. down-left (mm)")
        upLeftCol.SetName("Diag. up-left (mm)")
        downRightCol.SetName("Diag. down-right (mm)")

        distortionTableNode = self.getTableNode(
            "DistortionTablePanalWalls_" + modality + "_" + str(parentNodeID), parentNodeID
        )
        while distortionTableNode.RemoveRow(0):
            pass

        distortionTableNode.RemoveAllColumns()
        distortionTableNode.AddColumn(rightCol)
        distortionTableNode.AddColumn(leftCol)
        distortionTableNode.AddColumn(upRightCol)
        distortionTableNode.AddColumn(downLeftCol)
        distortionTableNode.AddColumn(upLeftCol)
        distortionTableNode.AddColumn(downRightCol)

    def createDistortionTableAxialDistances(self, parentNodeID, modality, axial_distance):
        axialDistanceCol = vtk.vtkDoubleArray()

        axialDistanceCol.SetName("Axial distances (mm)")

        distortionTableNode = self.getTableNode(
            "DistortionTableAxialDistances_" + modality + "_" + str(parentNodeID), parentNodeID
        )
        while distortionTableNode.RemoveRow(0):
            pass

        for i in range(0, len(axial_distance)):
            axialDistanceCol.InsertNextValue(axial_distance[i])

        distortionTableNode.RemoveAllColumns()
        distortionTableNode.AddColumn(axialDistanceCol)

    def createDistortionTableMeshSides(
        self,
        parentNodeID,
        modality,
        z_slices_filtered,
        right_sides,
        left_sides,
        up_sides,
        down_sides,
        centers_differences,
    ):
        if len(right_sides) != 0 and len(left_sides) != 0 and len(up_sides) != 0 and len(down_sides) != 0:
            centerDevX = [coord[0] for coord in centers_differences]
            centerDevY = [coord[1] for coord in centers_differences]

            sliceCol = vtk.vtkIntArray()
            sliceCol.SetName("Slice")

            rightCol = vtk.vtkDoubleArray()
            leftCol = vtk.vtkDoubleArray()
            upCol = vtk.vtkDoubleArray()
            downCol = vtk.vtkDoubleArray()
            meshCenterDevYCol = vtk.vtkDoubleArray()
            meshCenterDevXCol = vtk.vtkDoubleArray()

            rightCol.SetName("Right side (mm)")
            leftCol.SetName("Left side (mm)")
            upCol.SetName("Up side (mm)")
            downCol.SetName("Down side (mm)")
            meshCenterDevYCol.SetName("Mesh y-axis center deviation")
            meshCenterDevXCol.SetName("Mesh x-axis center deviation")

            distortionTableNode = self.getTableNode(
                "DistortionTableMeshSides_" + modality + "_" + str(parentNodeID), parentNodeID
            )
            while distortionTableNode.RemoveRow(0):
                pass

            for i in range(0, len(right_sides)):
                sliceCol.InsertNextValue(z_slices_filtered[i])
                rightCol.InsertNextValue(right_sides[i])
                leftCol.InsertNextValue(left_sides[i])
                upCol.InsertNextValue(up_sides[i])
                downCol.InsertNextValue(down_sides[i])
                meshCenterDevYCol.InsertNextValue(centerDevY[i])
                meshCenterDevXCol.InsertNextValue(centerDevX[i])

            distortionTableNode.RemoveAllColumns()
            distortionTableNode.AddColumn(sliceCol)
            distortionTableNode.AddColumn(rightCol)
            distortionTableNode.AddColumn(leftCol)
            distortionTableNode.AddColumn(upCol)
            distortionTableNode.AddColumn(downCol)
            distortionTableNode.AddColumn(meshCenterDevYCol)
            distortionTableNode.AddColumn(meshCenterDevXCol)

    def createResolutionTable(self, parentNodeID, sections, contrastCT, holesCT, rcHoles, modality):
        nameCol = vtk.vtkStringArray()
        nameCol.SetName("Section")
        contrastCtCol = vtk.vtkDoubleArray()
        contrastCtCol.SetName("Contrast")
        holesCtCol = vtk.vtkDoubleArray()
        holesCtCol.SetName("Holes (" + modality + "")
        rcHolesCtCol = vtk.vtkDoubleArray()
        rcHolesCtCol.SetName("RC Holes")
        diameterCol = vtk.vtkDoubleArray()
        diameterCol.SetName("Diameter (mm)")
        sectionDiameterMapping = {"1": 8.0, "2": 9.5, "3": 10.0, "4": 12.0, "5": 14.0, "6": 16.0}

        qunatificationTableNode = self.getResolutionTableNode(modality, parentNodeID)
        while qunatificationTableNode.RemoveRow(0):
            pass

        for i in range(0, len(sections)):
            nameCol.InsertNextValue(str(sections[i]))
            contrastCtCol.InsertNextValue(contrastCT[i])
            holesCtCol.InsertNextValue(holesCT[i])
            rcHolesCtCol.InsertNextValue(rcHoles[i])
            diameterCol.InsertNextValue(sectionDiameterMapping[sections[i]])

        qunatificationTableNode.RemoveAllColumns()
        qunatificationTableNode.AddColumn(nameCol)
        qunatificationTableNode.AddColumn(contrastCtCol)
        qunatificationTableNode.AddColumn(holesCtCol)
        qunatificationTableNode.AddColumn(rcHolesCtCol)
        qunatificationTableNode.AddColumn(diameterCol)

    def addResolutionResult(self, resolutionResult, volume, modalityName, unit):
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(volume))

        self.createResolutionTable(
            parentNodeID,
            resolutionResult[3],
            resolutionResult[0],
            resolutionResult[1],
            resolutionResult[2],
            modalityName,
        )
        self.resolutionPlot(
            resolutionResult[7],
            resolutionResult[4],
            resolutionResult[5],
            resolutionResult[8],
            resolutionResult[6],
            unit,
            modalityName,
        )
        if resolutionResult[10] is not None:
            self.resolutionPlot2(resolutionResult[10], modalityName)

        if resolutionResult[9] is not None:
            self.addSitkVolume(
                resolutionResult[9], "ResolutionSegmentedTriangles", "vtkMRMLLabelMapVolumeNode", parentNodeID
            )

    def resolutionCallback(self, resultCT, resultPET, resultT1, resultT2):
        if resultCT:
            self.addResolutionResult(resultCT, self.ctVolume, "CT", "HU")
        if resultPET:
            self.addResolutionResult(resultPET, self.petVolume, "PET", "Bq/mL")
        if resultT2:
            self.addResolutionResult(resultT2, self.t2Volume, "T2", "ms")
        if resultT1:
            self.addResolutionResult(resultT1, self.t1Volume, "T1", "ms")

    def addDistortionResultPlotPanalWalls(
        self,
        sitkImage,
        modality,
        distances_coord_r_mm,
        distances_coord_l_mm,
        distances_coord_dsu_mm,
        distances_coord_dsd_mm,
        distances_coord_dmu_mm,
        distances_coord_dmd_mm,
    ):
        walls_fig = plt.figure()
        walls_colors = ["#FF4444", "#55CC55", "#FFCC33", "#3366FF", "#CC6699", "#6666CC"]
        walls_labels = ["Right", "Left", "Diag. up-right", "Diag. down-left", "Diag. up-left", "Diag. down-right"]
        markers = ["o", "s", "D", "v", "h", "p"]
        ax = walls_fig.add_subplot(1, 1, 1)

        # Iterate over each type of side and its respective values
        for i, walls_list in enumerate(
            [
                distances_coord_r_mm,
                distances_coord_l_mm,
                distances_coord_dsu_mm,
                distances_coord_dsd_mm,
                distances_coord_dmu_mm,
                distances_coord_dmd_mm,
            ]
        ):
            walls_list_filtered = [value for value in walls_list if value is not None]
            mean_value = np.mean(walls_list_filtered)
            mean_std = np.std(walls_list_filtered) / np.sqrt(len(walls_list_filtered))

            ax.scatter(
                range(len(walls_list)), walls_list, color=walls_colors[i], label=walls_labels[i], marker=markers[i]
            )
            plt.axhline(
                y=mean_value,
                color=walls_colors[i],
                label=f"Mean {walls_labels[i]} side: {round(mean_value, 2)} \u00B1 {round(mean_std, 2)} mm",
            )

        ax.set_ylabel("Measured distance (mm)")
        ax.set_title(
            f"Distances between panal insert walls \n Pixel dimension: {round(sitkImage.GetSpacing()[0], 2)}x{round(sitkImage.GetSpacing()[1], 2)} mm\u00B2"
        )
        plt.ylim([23, 30])
        ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_xticks([])

        filename = os.path.join(self.outputDirectory, f"Distortion_{modality}_PanalWalls.png")
        walls_fig.savefig(filename, bbox_inches="tight")
        if self.showOutputPlots:
            pm = qt.QPixmap(filename)
            imageWidget = qt.QLabel()
            imageWidget.setPixmap(pm)
            imageWidget.setScaledContents(True)
            imageWidget.show()
            self.pltPlots[f"Distortion_{modality}_PanalWalls"] = imageWidget
        else:
            plt.close(walls_fig)

    def addDistortionResultPlotSideValues(
        self, sitkImage, modality, right_sides, left_sides, up_sides, down_sides, z_slices_filtered
    ):
        # Plots sides values
        sides_fig = plt.figure()
        sides_colors = ["red", "blue", "green", "orange"]
        sides_labels = ["Right", "Left", "Up", "Down"]
        markers = ["o", "s", "D", "v"]
        ax = sides_fig.add_subplot(1, 1, 1)
        # Iterate over each type of side and its respective values
        for i, side_list in enumerate([right_sides, left_sides, up_sides, down_sides]):
            side_list_filtered = [value for value in side_list if value is not None]
            mean_value = np.mean(side_list_filtered)
            mean_std = np.std(side_list_filtered) / np.sqrt(len(side_list_filtered))
            ax.scatter(z_slices_filtered, side_list, color=sides_colors[i], label=sides_labels[i], marker=markers[i])
            plt.axhline(
                y=mean_value,
                color=sides_colors[i],
                label=f"Mean {sides_labels[i]} side: {round(mean_value, 2)} \u00B1 {round(mean_std, 2)} mm",
            )

        ax.set_ylabel("Measured distance (mm)")
        ax.set_xlabel("Slice")
        ax.set_title(
            f"Distances along transaxial plane of the cylindrical insert \n Pixel dimension: {round(sitkImage.GetSpacing()[0], 2)}x{round(sitkImage.GetSpacing()[1], 2)} mm\u00B2"
        )
        ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.ylim([6, 14])

        filename = os.path.join(self.outputDirectory, f"Distortion_{modality}_MeshSides.png")
        plt.savefig(filename, bbox_inches="tight")
        if self.showOutputPlots:
            pm = qt.QPixmap(filename)
            imageWidget = qt.QLabel()
            imageWidget.setPixmap(pm)
            imageWidget.setScaledContents(True)
            imageWidget.show()
            self.pltPlots[f"Distortion_{modality}_MeshSides"] = imageWidget
        else:
            plt.close()

    def addDistortionResultPlotCentralDeviation(self, sitkImage, modality, centers_differences, z_slices_filtered):
        # Plots (y,x) mesh central deviation
        if len(centers_differences) > 0:
            central_fig = plt.figure()
            central_colors = ["blue", "green"]
            central_labels = ["y-axis", "x-axis"]
            markers = ["o", "s"]
            ax = central_fig.add_subplot(1, 1, 1)

            ax.scatter(
                z_slices_filtered,
                [coord[0] for coord in centers_differences],
                color=central_colors[0],
                label=central_labels[0],
                marker=markers[0],
            )
            ax.scatter(
                z_slices_filtered,
                [coord[1] for coord in centers_differences],
                color=central_colors[1],
                label=central_labels[1],
                marker=markers[1],
            )

            max_y = max(coord[0] for coord in centers_differences)
            min_y = min(coord[0] for coord in centers_differences)
            max_x = max(coord[1] for coord in centers_differences)
            min_x = min(coord[1] for coord in centers_differences)
            upper_limit = max(max_y, max_x)
            lower_limit = min(min_y, min_x)
            margin = 10
            upper_limit += margin
            lower_limit -= margin

            ax.set_ylabel("Mesh centres deviation")
            ax.set_xlabel("Slice")
            ax.set_title(
                f"Mesh centres deviation \n Pixel dimension: {round(sitkImage.GetSpacing()[0], 2)}x{round(sitkImage.GetSpacing()[1], 2)} mm\u00B2"
            )
            ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc="upper left")
            ax.set_ylim([lower_limit, upper_limit])

            filename = os.path.join(self.outputDirectory, f"Distortion_{modality}_MeshCentres.png")
            central_fig.savefig(filename, bbox_inches="tight")
            if self.showOutputPlots:
                pm = qt.QPixmap(filename)
                imageWidget = qt.QLabel()
                imageWidget.setPixmap(pm)
                imageWidget.setScaledContents(True)
                imageWidget.show()
                self.pltPlots[f"Distortion_{modality}_MeshCentres"] = imageWidget
            else:
                plt.close(central_fig)

    def addDistortionResultPlotAxialDistances(
        self,
        sitkImage,
        axial_distance,
        modality,
    ):
        mean_axial_distance = np.mean(axial_distance)
        max_axial_distance = np.max(axial_distance)
        min_axial_distance = np.min(axial_distance)
        std_value = np.std(axial_distance)
        # Plots axial distances
        plt.figure()
        plt.scatter(range(len(axial_distance)), axial_distance, marker="o")
        plt.axhline(y=mean_axial_distance, color="green", label="Mean distance")
        plt.ylim([min_axial_distance - 2, max_axial_distance + 3])
        plt.legend(
            [
                f"Distances between axial planes",
                f"Mean distance in axial axis: {round(mean_axial_distance, 2)} \u00B1 {round(std_value, 2)} mm",
            ]
        )
        plt.ylabel("Measured distance (mm)")
        plt.title(
            f"Distances along superior-inferior direction of cylindrical insert \n Slice spacing: {round(sitkImage.GetSpacing()[2], 2)} mm"
        )
        plt.xticks([])

        filename = os.path.join(self.outputDirectory, f"Distortion_{modality}_AxialDistances.png")
        plt.savefig(filename, bbox_inches="tight")
        if self.showOutputPlots:
            pm = qt.QPixmap(filename)
            imageWidget = qt.QLabel()
            imageWidget.setPixmap(pm)
            imageWidget.setScaledContents(True)
            imageWidget.show()
            self.pltPlots[f"Distortion_{modality}_AxialDistances"] = imageWidget
        else:
            plt.close()

    def registrationCallback(self, resultCT, resultT1, resultT2):
        if resultCT:
            self.registrationPlot(resultCT[0], resultCT[1], "CT")
        if resultT2:
            self.registrationPlot(
                resultT2[0],
                resultT2[1],
                "T2",
            )
        if resultT1:
            self.registrationPlot(resultT1[0], resultT1[1], "T1")

    def addDistortionResults(self, volume, modality, results):
        sitkImage = sitkUtils.PullVolumeFromSlicer(volume)
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        parentNodeID = shNode.GetItemParent(shNode.GetItemByDataNode(volume))
        axial_distance, meshResult, panalInsertResult = results
        right_sides, left_sides, up_sides, down_sides, z_slices_filtered, _, centers_differences = meshResult

        self.addDistortionResultPlotAxialDistances(sitkImage, axial_distance, modality)
        self.addDistortionResultPlotCentralDeviation(sitkImage, modality, centers_differences, z_slices_filtered)
        self.addDistortionResultPlotPanalWalls(sitkImage, modality, *panalInsertResult)
        self.addDistortionResultPlotSideValues(
            sitkImage, modality, right_sides, left_sides, up_sides, down_sides, z_slices_filtered
        )

        self.createDistortionTableAxialDistances(parentNodeID, modality, axial_distance)
        self.createDistortionTableMeshSides(
            parentNodeID,
            modality,
            z_slices_filtered,
            right_sides,
            left_sides,
            up_sides,
            down_sides,
            centers_differences,
        )
        self.createDistortionTablePanalWalls(parentNodeID, modality, *panalInsertResult)

    def distortionCallBack(self, ctResult, t1Result, t2Result):
        if ctResult is not None:
            self.addDistortionResults(self.ctVolume, "CT", ctResult)

        if t1Result is not None:
            self.addDistortionResults(self.t1Volume, "T1", t1Result)

        if t2Result is not None:
            self.addDistortionResults(self.t2Volume, "T2", t2Result)

    def reportProgress(self, progressVal, callingMethod):
        if not hasattr(slicer, "progressWindow") or slicer.progressWindow.wasCanceled:
            slicer.progressWindow = slicer.util.createProgressDialog()
        slicer.progressWindow.show()
        slicer.progressWindow.activateWindow()
        slicer.progressWindow.setValue(progressVal)
        slicer.progressWindow.setWindowTitle(f"{callingMethod}...")
        slicer.app.processEvents()

    def doRadiomics(self):
        radiomicsThread = RadiomicsRoutine(
            self.ctVolume,
            self.petVolume,
            self.t1Volume,
            self.t2Volume,
            self.ctSegment,
            self.petSegment,
            self.t1Segment,
            self.t2Segment,
            callBackFunction=self.radiomicsCallback,
            statusUpdateFunction=self.reportProgress,
            resourcePath=self.resourcePath,
        )
        qt.QTimer.singleShot(0, radiomicsThread.run)

    def distortionStartEndSliceCylinderDetection(self, currentNode, modality, resourcePath, shrinkImage=False):
        sitkCtVolumeOrig = sitkUtils.PullVolumeFromSlicer(currentNode)
        referenceImageFileName = os.path.join(resourcePath, modality, f"{modality}.nrrd")
        referenceSliceIdxsFileName = os.path.join(resourcePath, modality, "distortionSlices.csv")
        im_ref = sitk.ReadImage(referenceImageFileName, sitk.sitkFloat32)  # reference image
        sitkCtVolume = sitk.Cast(sitkCtVolumeOrig, im_ref.GetPixelID())

        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(sitkCtVolume)
        resampler.SetInterpolator(sitk.sitkLinear)
        resampler.SetDefaultPixelValue(0)
        resampled_moving_image = resampler.Execute(im_ref)

        if shrinkImage:
            sitkCtVolume = sitk.Shrink(sitkCtVolume, [2, 2, 2])
            resampled_moving_image = sitk.Shrink(resampled_moving_image, [2, 2, 2])
        transformation = af.register(
            sitkCtVolume,
            resampled_moving_image,
            0.01,
            shrinkFactors=[2, 1],
            smoothingSigmas=[1, 0],
        )

        transformation = transformation.GetInverse()
        with open(referenceSliceIdxsFileName, "r") as file:
            csv_reader = csv.reader(file, delimiter=",")
            cylinderStart, cylinderEnd, panalSlice = next(csv_reader)
            cylinderStart = int(cylinderStart)
            cylinderEnd = int(cylinderEnd)
            panalSlice = int(panalSlice)
            cylinerStartWorld = im_ref.TransformIndexToPhysicalPoint((0, 0, cylinderStart))
            cylinerStartWorld = transformation.TransformPoint(cylinerStartWorld)
            _, _, cylinderStart = sitkCtVolumeOrig.TransformPhysicalPointToIndex(cylinerStartWorld)

            cylinderEndWorld = im_ref.TransformIndexToPhysicalPoint((0, 0, cylinderEnd))
            cylinderEndWorld = transformation.TransformPoint(cylinderEndWorld)
            _, _, cylinderEnd = sitkCtVolumeOrig.TransformPhysicalPointToIndex(cylinderEndWorld)

            panalSliceWorld = im_ref.TransformIndexToPhysicalPoint((0, 0, panalSlice))
            panalSliceWorld = transformation.TransformPoint(panalSliceWorld)
            _, _, panalSlice = sitkCtVolumeOrig.TransformPhysicalPointToIndex(panalSliceWorld)

        if modality == "T1" or modality == "T2":
            referenceDistortionWallPointsFileName = os.path.join(
                resourcePath, modality, f"{modality}DistortionWallPointsRight.mrk.json"
            )
            nodeName = f"DistortionWallPoints_{modality}_Right"
            self.loadMarkupFileList(referenceDistortionWallPointsFileName, nodeName, currentNode, transformation)

            referenceDistortionWallPointsFileName = os.path.join(
                resourcePath, modality, f"{modality}DistortionWallPointsLeft.mrk.json"
            )
            nodeName = f"DistortionWallPoints_{modality}_Left"
            self.loadMarkupFileList(referenceDistortionWallPointsFileName, nodeName, currentNode, transformation)

            referenceDistortionWallPointsFileName = os.path.join(
                resourcePath, modality, f"{modality}DistortionWallPointsDiagonalUpRight.mrk.json"
            )
            nodeName = f"DistortionWallPoints_{modality}_DiagonalUpRight"
            self.loadMarkupFileList(referenceDistortionWallPointsFileName, nodeName, currentNode, transformation)

            referenceDistortionWallPointsFileName = os.path.join(
                resourcePath, modality, f"{modality}DistortionWallPointsDiagonalUpLeft.mrk.json"
            )
            nodeName = f"DistortionWallPoints_{modality}_DiagonalUpLeft"
            self.loadMarkupFileList(referenceDistortionWallPointsFileName, nodeName, currentNode, transformation)

            referenceDistortionWallPointsFileName = os.path.join(
                resourcePath, modality, f"{modality}DistortionWallPointsDiagonalDownRight.mrk.json"
            )
            nodeName = f"DistortionWallPoints_{modality}_DiagonalDownRight"
            self.loadMarkupFileList(referenceDistortionWallPointsFileName, nodeName, currentNode, transformation)

            referenceDistortionWallPointsFileName = os.path.join(
                resourcePath, modality, f"{modality}DistortionWallPointsDiagonalDownLeft.mrk.json"
            )
            nodeName = f"DistortionWallPoints_{modality}_DiagonalDownLeft"
            self.loadMarkupFileList(referenceDistortionWallPointsFileName, nodeName, currentNode, transformation)

        return cylinderStart, cylinderEnd, panalSlice

    def loadMarkupFileList(self, path, nodeName, parentNode, transformation):
        markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
        for node in markupNodeList:
            slicer.mrmlScene.RemoveNode(node)
        markupNode = slicer.util.loadMarkups(path)
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        shNode.CreateItem(shNode.GetItemByDataNode(parentNode), markupNode)
        markupNode.SetName(nodeName)
        meanX = 0.0
        meanY = 0.0
        meanZ = 0.0
        for controlPoindIdx in range(0, markupNode.GetNumberOfControlPoints()):
            position = [0.0, 0.0, 0.0]
            markupNode.GetNthControlPointPositionWorld(controlPoindIdx, position)
            position[0] = position[0] * -1
            position[1] = position[1] * -1
            position = transformation.TransformPoint(position)
            position = [position[0] * -1, position[1] * -1, position[2]]
            meanX += position[0]
            meanY += position[1]
            meanZ += position[2]
            markupNode.SetNthControlPointPositionWorld(controlPoindIdx, position)

        slicer.modules.markups.logic().JumpSlicesToLocation(meanX / 4.0, meanY / 4.0, meanZ / 4.0, True)

    def resolutionSliceDetection(self, currentNode, modality, resourcePath, shrinkImage=False):
        sitkCtVolume = sitkUtils.PullVolumeFromSlicer(currentNode)
        referenceImageFileName = os.path.join(resourcePath, modality, f"{modality}.nrrd")
        im_ref = sitk.ReadImage(referenceImageFileName, sitk.sitkFloat32)  # reference image
        sitkCtVolume = sitk.Cast(sitkCtVolume, im_ref.GetPixelID())

        # movingImgVolumeNode = self.addSitkVolume(
        #     im_ref, "ReferenceImageResolutionSliceDetection", "vtkMRMLScalarVolumeNode", None
        # )

        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(sitkCtVolume)
        resampler.SetInterpolator(sitk.sitkLinear)
        resampler.SetDefaultPixelValue(0)
        resampled_moving_image = resampler.Execute(im_ref)

        # transformation = af.register(
        #     sitk.Shrink(sitkCtVolume, [2, 2, 2]),
        #     sitk.Shrink(resampled_moving_image, [2, 2, 2]),
        #     0.01,
        #     metric="ssm",
        #     shrinkFactors=[2, 1],
        #     smoothingSigmas=[1, 0],
        # )
        if shrinkImage:
            sitkCtVolume = sitk.Shrink(sitkCtVolume, [2, 2, 2])
            resampled_moving_image = sitk.Shrink(resampled_moving_image, [2, 2, 2])
        transformation = af.register(
            sitkCtVolume,
            resampled_moving_image,
            0.01,
            shrinkFactors=[2, 1],
            smoothingSigmas=[1, 0],
        )

        # transformedMovingVolumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode")
        # transformNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode")

        # parameters = {}
        # parameters["fixedVolume"] = currentNode.GetID()
        # parameters["movingVolume"] = movingImgVolumeNode.GetID()
        # parameters["outputVolume"] = transformedMovingVolumeNode.GetID()
        # parameters["linearTransform"] = transformNode.GetID()
        # parameters["useRigid"] = True  # options include: "useRigid", "useAffine", "useBSpline"
        # parameters["initializeTransformMode"] = "useGeometryAlign"
        # parameters["samplingPercentage"] = 0.02
        # cliBrainsFitRigidNode = slicer.cli.run(slicer.modules.brainsfit, None, parameters, wait_for_completion=True)

        transformation = transformation.GetInverse()
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        meanX = 0.0
        meanY = 0.0
        meanZ = 0.0
        for triangleIdx in range(1, 7):
            referenceTriangleFileName = os.path.join(
                resourcePath, modality, f"{modality}Triangle{triangleIdx}.mrk.json"
            )

            nodeName = f"ResolutionTriangle_{modality}_{triangleIdx}"
            markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
            for node in markupNodeList:
                slicer.mrmlScene.RemoveNode(node)
            markupNode = slicer.util.loadMarkups(referenceTriangleFileName)
            shNode.CreateItem(shNode.GetItemByDataNode(currentNode), markupNode)
            markupNode.SetName(nodeName)

            for controlPoindIdx in range(0, markupNode.GetNumberOfControlPoints()):
                position = [0.0, 0.0, 0.0]
                markupNode.GetNthControlPointPositionWorld(controlPoindIdx, position)
                position[0] = position[0] * -1
                position[1] = position[1] * -1
                position = transformation.TransformPoint(position)
                position = [position[0] * -1, position[1] * -1, position[2]]
                meanX += position[0]
                meanY += position[1]
                meanZ += position[2]
                markupNode.SetNthControlPointPositionWorld(controlPoindIdx, position)

        slicer.modules.markups.logic().JumpSlicesToLocation(meanX / 18.0, meanY / 18.0, meanZ / 18.0, True)

    def doQuantification(self):
        activityDecay = self.calculateActivityDecay()

        quantThread = QuantificationRoutine(
            self.ctVolume,
            self.petVolume,
            self.t1Volume,
            self.t2Volume,
            activityDecay,
            self.insertType,
            self.quantificationStartSlicePET,
            self.quantificationStartSliceCT,
            self.quantificationStartSliceT1,
            self.quantificationStartSliceT2,
            callBackFunction=self.quantificationCallback,
            statusUpdateFunction=self.reportProgress,
        )
        qt.QTimer.singleShot(0, quantThread.run)

    def insertTypeCalibrationCurve(self):
        return None

    def insertTypeMaterialsCharacterization(self):
        return None

    def getResolutionInsertPointPositions(self, modality, imgVolume):
        triangles = []
        if imgVolume is not None:
            sitkVolume = sitkUtils.PullVolumeFromSlicer(imgVolume)
            z_coord = None
            for triangleIdx in range(1, 7):
                nodeName = f"ResolutionTriangle_{modality}_{triangleIdx}"
                markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
                for node in markupNodeList:
                    if node.GetNumberOfControlPoints() == 3:
                        nodePositions = []
                        for controlPoindIdx in range(0, node.GetNumberOfControlPoints()):
                            position = [0.0, 0.0, 0.0]
                            node.GetNthControlPointPositionWorld(controlPoindIdx, position)
                            position[0] = position[0] * -1
                            position[1] = position[1] * -1
                            x, y, z_coord = sitkVolume.TransformPhysicalPointToIndex(position)
                            nodePositions.append([y, x])
                        triangles.append(nodePositions)
        if len(triangles) == 6:
            return z_coord, triangles
        else:
            return None, None

    def getNodePoints(self, nodeName, sitkVolume):
        points = []
        markupNodeList = slicer.mrmlScene.GetNodesByName(nodeName)
        for node in markupNodeList:
            nodePositions = []
            for controlPoindIdx in range(0, node.GetNumberOfControlPoints()):
                position = [0.0, 0.0, 0.0]
                node.GetNthControlPointPositionWorld(controlPoindIdx, position)
                position[0] = position[0] * -1
                position[1] = position[1] * -1
                x, y, z_coord = sitkVolume.TransformPhysicalPointToIndex(position)
                points.append([y, x])
        return points

    def getDistortionPanalInsertPointPositions(self, modality, imgVolume):
        pointsRight = pointsLeft = pointsUpRight = pointsUpLeft = pointsDownRight = pointsDownLeft = None
        if imgVolume is not None:
            sitkVolume = sitkUtils.PullVolumeFromSlicer(imgVolume)
            nodeName = f"DistortionWallPoints_{modality}_Right"
            pointsRight = self.getNodePoints(nodeName, sitkVolume)

            nodeName = f"DistortionWallPoints_{modality}_Left"
            pointsLeft = self.getNodePoints(nodeName, sitkVolume)

            nodeName = f"DistortionWallPoints_{modality}_DiagonalUpRight"
            pointsUpRight = self.getNodePoints(nodeName, sitkVolume)

            nodeName = f"DistortionWallPoints_{modality}_DiagonalUpLeft"
            pointsUpLeft = self.getNodePoints(nodeName, sitkVolume)

            nodeName = f"DistortionWallPoints_{modality}_DiagonalDownRight"
            pointsDownRight = self.getNodePoints(nodeName, sitkVolume)

            nodeName = f"DistortionWallPoints_{modality}_DiagonalDownLeft"
            pointsDownLeft = self.getNodePoints(nodeName, sitkVolume)
        return pointsRight, pointsLeft, pointsUpRight, pointsUpLeft, pointsDownRight, pointsDownLeft

    def doResolution(self):
        resolutionTrianglesCT = self.getResolutionInsertPointPositions("CT", self.ctVolume)
        resolutionTrianglesPET = self.getResolutionInsertPointPositions("PET", self.petVolume)
        resolutionTrianglesT1 = self.getResolutionInsertPointPositions("T1", self.t1Volume)
        resolutionTrianglesT2 = self.getResolutionInsertPointPositions("T2", self.t2Volume)
        quantThread = ResolutionRoutine(
            self.ctVolume,
            self.petVolume,
            self.t1Volume,
            self.t2Volume,
            resolutionTrianglesCT,
            resolutionTrianglesPET,
            resolutionTrianglesT1,
            resolutionTrianglesT2,
            callBackFunction=self.resolutionCallback,
            statusUpdateFunction=self.reportProgress,
        )
        qt.QTimer.singleShot(0, quantThread.run)

    def doRegistration(self):
        registrationThread = RegistrationRoutine(
            self.ctVolume,
            self.petVolume,
            self.t1Volume,
            self.t2Volume,
            callBackFunction=self.registrationCallback,
            statusUpdateFunction=self.reportProgress,
        )
        qt.QTimer.singleShot(0, registrationThread.run)

    def doDistortion(self):
        pointsT1 = self.getDistortionPanalInsertPointPositions("T1", self.t1Volume)
        pointsT2 = self.getDistortionPanalInsertPointPositions("T2", self.t2Volume)
        distortionThread = DistortionRoutine(
            self.ctVolume,
            self.t1Volume,
            self.t2Volume,
            self.distortionFirstSliceCT,
            self.distortionLastSliceCT,
            self.distortionFirstSliceT1,
            self.distortionLastSliceT1,
            self.distortionFirstSliceT2,
            self.distortionLastSliceT2,
            self.panalSliceCT,
            self.panalSliceT1,
            self.panalSliceT2,
            pointsT1,
            pointsT2,
            self.distortionInsertFilled,
            callBackFunction=self.distortionCallBack,
            statusUpdateFunction=self.reportProgress,
        )
        qt.QTimer.singleShot(0, distortionThread.run)

    def process(
        self,
        petVolume,
        ctVolume,
        t1Volume,
        t2Volume,
        petSegment,
        ctSegment,
        t1Segment,
        t2Segment,
        rcVolume,
        decayFactor,
        selectedInsertType,
        selectedTracerType,
        initalActivity,
        timeActivity,
        insertPositions,
        checkedMethods,
        resourcePath,
        quantificationStartSlicePET,
        quantificationStartSliceCT,
        quantificationStartSliceT1,
        quantificationStartSliceT2,
        distortionFirstSliceCT,
        distortionLastSliceCT,
        panalSliceCT,
        panalSliceT1,
        panalSliceT2,
        distortionFirstSliceT1,
        distortionLastSliceT1,
        distortionFirstSliceT2,
        distortionLastSliceT2,
        distortionInsertFilled,
    ):
        analyisMethods = self.getSupportedAnalysesMethods()
        self.rcVolume = rcVolume
        self.insertType = selectedInsertType
        self.decayFactor = decayFactor
        self.tracerType = selectedTracerType
        self.initialActivity = initalActivity
        self.timeActivity = timeActivity
        self.insertPositions = insertPositions
        self.resourcePath = resourcePath

        self.quantificationStartSlicePET = quantificationStartSlicePET
        self.quantificationStartSliceCT = quantificationStartSliceCT
        self.quantificationStartSliceT1 = quantificationStartSliceT1
        self.quantificationStartSliceT2 = quantificationStartSliceT2

        self.distortionFirstSliceCT = distortionFirstSliceCT
        self.distortionLastSliceCT = distortionLastSliceCT
        self.distortionFirstSliceT1 = distortionFirstSliceT1
        self.distortionLastSliceT1 = distortionLastSliceT1
        self.distortionFirstSliceT2 = distortionFirstSliceT2
        self.distortionLastSliceT2 = distortionLastSliceT2
        self.panalSliceCT = panalSliceCT
        self.panalSliceT1 = panalSliceT1
        self.panalSliceT2 = panalSliceT2
        self.distortionInsertFilled = distortionInsertFilled

        self.petVolume = petVolume
        self.ctVolume = ctVolume
        self.t1Volume = t1Volume
        self.t2Volume = t2Volume

        self.petSegment = petSegment
        self.ctSegment = ctSegment
        self.t1Segment = t1Segment
        self.t2Segment = t2Segment

        for checkedMethod in checkedMethods:
            methodToCall = analyisMethods[checkedMethod]
            methodToCall()


class SegmentationRoutine(object):
    def __init__(self, ctVolume, petVolume, t1Volume, t2Volume, callBackFunction=None, statusUpdateFunction=None):
        super().__init__()
        self.ctVolume = ctVolume
        self.petVolume = petVolume
        self.t1Volume = t1Volume
        self.t2Volume = t2Volume
        self.callbackFunction = callBackFunction
        self.statusUpdateFunction = statusUpdateFunction

    def getSitkImageFromNumpyArray(self, numpyArray, refrenceImage):
        sitkImage = sitk.GetImageFromArray(numpyArray)
        sitkImage.SetDirection(refrenceImage.GetDirection())
        sitkImage.SetOrigin(refrenceImage.GetOrigin())
        sitkImage.SetSpacing(refrenceImage.GetSpacing())
        return sitkImage

    def increaseStatusUpdateFunction(self, value: int, calingMethond: str):
        if self.statusUpdateFunction is not None:
            self.statusUpdateFunction(value, calingMethond)


class DistortionRoutine(SegmentationRoutine):
    def __init__(
        self,
        ctVolume,
        t1Volume,
        t2Volume,
        startSliceCT,
        endSliceCT,
        startSliceT1,
        endSliceT1,
        startSliceT2,
        endSliceT2,
        panalSliceCT,
        panalSliceT1,
        panalSliceT2,
        pointsT1,
        pointsT2,
        insertFilled,
        callBackFunction=None,
        statusUpdateFunction=None,
    ):
        super().__init__(ctVolume, None, t1Volume, t2Volume, callBackFunction, statusUpdateFunction)
        self.startSliceCT = startSliceCT
        self.endSliceCT = endSliceCT
        self.startSliceT1 = startSliceT1
        self.endSliceT1 = endSliceT1
        self.startSliceT2 = startSliceT2
        self.endSliceT2 = endSliceT2
        self.pointsT1 = pointsT1
        self.pointsT2 = pointsT2
        self.insertFilled = insertFilled
        self.panalSliceCT = panalSliceCT
        self.panalSliceT1 = panalSliceT1
        self.panalSliceT2 = panalSliceT2

    def getMeshVertices(self, thresh_cyl_sitk, z_slices, erosionIterations=3):
        # OBTAINS MESH VERTICES (xy coordinates)
        thresh_cyl = sitk.GetArrayFromImage(thresh_cyl_sitk)
        fused_vertices = []
        for z in z_slices:
            z_thresh = thresh_cyl[z, :, :]  # image slice
            # Proceeds if there is any non-null pixel
            if np.any(z_thresh != 0):
                # Erosion of the segmentation to remove the wall thickness of the squares
                eroded_segmentation = ndimage.morphology.binary_erosion(z_thresh, iterations=erosionIterations)
                # plt.imshow(eroded_segmentation, cmap='gray')
                # plt.show()

                # Get coordinates of pixels that are part of the eroded segmentation
                y_coords, x_coords = np.nonzero(eroded_segmentation)

                if len(y_coords) == 0 or len(x_coords) == 0:
                    continue
                else:
                    # Get the coordinates of the vertices of each square of the mesh
                    vertices = [(z, y, x) for y, x in zip(y_coords, x_coords)]

                    # Calculate distance matrix between vertices
                    dist_matrix = spatial.distance_matrix(vertices, vertices)

                    n_vertices = len(vertices)
                    visited_vertices = [False] * n_vertices

                    for i in range(n_vertices):
                        if visited_vertices[i]:
                            continue

                        # Identify the vertices near this vertex
                        close_vertices = np.where(dist_matrix[i] < 3)[0]

                        # Take the average value of the coordinates of nearby vertices
                        fused_vertex = np.mean(np.array([vertices[j] for j in close_vertices]), axis=0)

                        # Check if the merged coordinate is already present in fused_vertices or if any of the contiguous coordinates are present
                        contiguous_vertices = [
                            vertex for vertex in fused_vertices if np.linalg.norm(vertex - fused_vertex) < 2
                        ]
                        if not contiguous_vertices:
                            fused_vertices.append(fused_vertex)  # Add merged coordinate to the list
                        else:
                            averaged_vertex = np.mean(
                                np.concatenate((contiguous_vertices, np.array([fused_vertex]))), axis=0
                            )
                            fused_vertices = [
                                v
                                for v in fused_vertices
                                if not any(np.array_equal(v, vertex) for vertex in contiguous_vertices)
                            ]  # Delete contiguous_vertices from fused_vertices
                            fused_vertices.append(averaged_vertex)  # Add the average to the list fused_vertices

                        # Mark visited vertices as merged
                        close_vertices = np.array(close_vertices).astype(int)
                        visited_vertices = np.array(visited_vertices).astype(int)
                        visited_vertices[close_vertices] = True
        vertices_array_fused = np.array(fused_vertices).astype(int)
        mesh = np.zeros_like(thresh_cyl, dtype=int)
        mesh[vertices_array_fused[:, 0], vertices_array_fused[:, 1], vertices_array_fused[:, 2]] = 1
        mesh_sitk = sitk.GetImageFromArray(mesh)
        mesh_sitk.SetDirection(thresh_cyl_sitk.GetDirection())
        mesh_sitk.SetSpacing(thresh_cyl_sitk.GetSpacing())
        mesh_sitk.SetOrigin(thresh_cyl_sitk.GetOrigin())
        return mesh_sitk

    def meshAnalysisMR(self, mesh_sitk, z_slices):
        mesh_point = sitk.GetArrayFromImage(mesh_sitk)

        right_sides = []
        left_sides = []
        up_sides = []
        down_sides = []
        z_slices_filtered = []
        mesh_centers = []
        centers_differences = []
        for i in z_slices:
            plane = mesh_point[i, :, :]
            y_nonzero, x_nonzero = np.where(plane == 1)
            if len(y_nonzero) == 9:
                z_slices_filtered.append(i)
                nonzero_coords = np.where(plane == 1)
                combined_nonzero_coords = list(zip(nonzero_coords[0], nonzero_coords[1]))
                y_nonzero_mean, x_nonzero_mean = np.mean(y_nonzero), np.mean(x_nonzero)

                distances = [
                    np.sqrt((coord[0] - y_nonzero_mean) ** 2 + (coord[1] - x_nonzero_mean) ** 2)
                    for coord in combined_nonzero_coords
                ]
                y_central, x_central = combined_nonzero_coords[np.argmin(distances)]
                mesh_centers.append([y_central, x_central])

                # Calculates deviation of (y,x) from first
                first_y, first_x = mesh_centers[0]
                difference_yx = [y_central - first_y, x_central - first_x]
                centers_differences.append(difference_yx)

                # Choose the neighbours closest to the central point
                neigh_coords = [
                    coord
                    for coord in combined_nonzero_coords
                    if np.sqrt(
                        ((coord[0] - y_central) * mesh_sitk.GetSpacing()[1]) ** 2
                        + ((coord[1] - x_central) * mesh_sitk.GetSpacing()[0]) ** 2
                    )
                    <= 12
                ]

                # Calculate horizontal sides (right, left) and vertical sides (up, down)
                for n_coord in neigh_coords:
                    y_diff = n_coord[0] - y_central
                    x_diff = n_coord[1] - x_central

                    # Right
                    if (
                        x_diff > 0
                        and abs(x_diff) >= (3 / mesh_sitk.GetSpacing()[0])
                        and abs(y_diff) <= (3 / mesh_sitk.GetSpacing()[1])
                    ):
                        x_diff = abs(x_diff * mesh_sitk.GetSpacing()[0])
                        right_sides.append(x_diff)
                    # Left
                    elif (
                        x_diff < 0
                        and abs(x_diff) >= (3 / mesh_sitk.GetSpacing()[0])
                        and abs(y_diff) <= (3 / mesh_sitk.GetSpacing()[1])
                    ):
                        x_diff = abs(x_diff * mesh_sitk.GetSpacing()[0])
                        left_sides.append(x_diff)
                    # Up
                    if (
                        y_diff < 0
                        and abs(y_diff) >= (3 / mesh_sitk.GetSpacing()[1])
                        and abs(x_diff) <= (3 / mesh_sitk.GetSpacing()[0])
                    ):
                        y_diff = abs(y_diff * mesh_sitk.GetSpacing()[1])
                        up_sides.append(y_diff)
                    # Down
                    elif (
                        y_diff > 0
                        and abs(y_diff) >= (3 / mesh_sitk.GetSpacing()[1])
                        and abs(x_diff) <= (3 / mesh_sitk.GetSpacing()[0])
                    ):
                        y_diff = abs(y_diff * mesh_sitk.GetSpacing()[1])
                        down_sides.append(y_diff)

                all_sides = [right_sides, left_sides, up_sides, down_sides]
                max_size = max(len(side_list) for side_list in all_sides)
                if len(right_sides) < max_size:
                    right_sides.extend([None] * (max_size - len(right_sides)))
                if len(left_sides) < max_size:
                    left_sides.extend([None] * (max_size - len(left_sides)))
                if len(up_sides) < max_size:
                    up_sides.extend([None] * (max_size - len(up_sides)))
                if len(down_sides) < max_size:
                    down_sides.extend([None] * (max_size - len(down_sides)))

        return right_sides, left_sides, up_sides, down_sides, z_slices_filtered, mesh_centers, centers_differences

    def meshAnalysisCT(self, mesh_sitk, z_slices):
        mesh_point = sitk.GetArrayFromImage(mesh_sitk)

        right_sides = []
        left_sides = []
        up_sides = []
        down_sides = []
        z_slices_filtered = []
        mesh_centers = []
        centers_differences = []
        for i in z_slices:
            plane = mesh_point[i, :, :]
            y_nonzero, x_nonzero = np.where(plane == 1)
            if len(y_nonzero) >= 7 and len(y_nonzero) <= 9:
                z_slices_filtered.append(i)
                nonzero_coords = np.where(plane == 1)
                combined_nonzero_coords = list(zip(nonzero_coords[0], nonzero_coords[1]))
                y_nonzero_mean, x_nonzero_mean = np.mean(y_nonzero), np.mean(x_nonzero)

                distances = [
                    np.sqrt((coord[0] - y_nonzero_mean) ** 2 + (coord[1] - x_nonzero_mean) ** 2)
                    for coord in combined_nonzero_coords
                ]
                y_central, x_central = combined_nonzero_coords[np.argmin(distances)]
                mesh_centers.append([y_central, x_central])

                # Calculates deviation of (y,x) from first
                first_y, first_x = mesh_centers[0]
                difference_yx = [y_central - first_y, x_central - first_x]
                centers_differences.append(difference_yx)

                # Choose the neighbours closest to the central point
                neigh_coords = [
                    coord
                    for coord in combined_nonzero_coords
                    if np.sqrt(
                        ((coord[0] - y_central) * mesh_sitk.GetSpacing()[1]) ** 2
                        + ((coord[1] - x_central) * mesh_sitk.GetSpacing()[0]) ** 2
                    )
                    <= 16
                ]

                # Calculate horizontal sides (x-axis) and vertical sides (y-axis):
                for n_coord in neigh_coords:
                    y_diff = n_coord[0] - y_central
                    x_diff = n_coord[1] - x_central

                    # Horizontal sides
                    if abs(y_diff) <= (3 / mesh_sitk.GetSpacing()[1]) and abs(x_diff) >= (
                        3 / mesh_sitk.GetSpacing()[0]
                    ):
                        if x_diff > 0:
                            right_sides.append(abs(x_diff * mesh_sitk.GetSpacing()[0]))
                        elif x_diff < 0:
                            left_sides.append(abs(x_diff * mesh_sitk.GetSpacing()[0]))
                    # Vertical sides
                    if abs(x_diff) <= (3 / mesh_sitk.GetSpacing()[0]) and abs(y_diff) >= (
                        3 / mesh_sitk.GetSpacing()[1]
                    ):
                        if y_diff < 0:
                            up_sides.append(abs(y_diff * mesh_sitk.GetSpacing()[1]))
                        elif y_diff > 0:
                            down_sides.append(abs(y_diff * mesh_sitk.GetSpacing()[1]))

                all_sides = [right_sides, left_sides, up_sides, down_sides]
                max_size = max(len(side_list) for side_list in all_sides)
                if len(right_sides) < max_size:
                    right_sides.extend([None] * (max_size - len(right_sides)))
                if len(left_sides) < max_size:
                    left_sides.extend([None] * (max_size - len(left_sides)))
                if len(up_sides) < max_size:
                    up_sides.extend([None] * (max_size - len(up_sides)))
                if len(down_sides) < max_size:
                    down_sides.extend([None] * (max_size - len(down_sides)))

        return right_sides, left_sides, up_sides, down_sides, z_slices_filtered, mesh_centers, centers_differences

    def cylindricalInsertAnalysis(self, thresh_cyl, threshold):
        z_planes = []  # list with cyl insert slices (float)
        z_slices = []  # list with cyl insert slices (int)
        # Iterate along dimension z to identify contiguous insert's slices
        for z in range(thresh_cyl.shape[0] - 1):
            current_slice = thresh_cyl[z, :, :]
            next_slice = thresh_cyl[z + 1, :, :]
            # Check for non-zero values in both slices that exceed the threshold
            if np.count_nonzero(current_slice) > threshold and np.count_nonzero(next_slice) > threshold:
                # Combines both slices in binary form
                combined_slice = current_slice + next_slice
                new_slice = np.where(combined_slice >= 1, 1, combined_slice)
                # Calculate the average of z
                averaged_z = int((z + (z + 1)) / 2)
                z_planes.append(((z + (z + 1)) / 2))
                z_slices.append(averaged_z)
                # Update the segmentation values
                thresh_cyl[z, :, :] = 0
                thresh_cyl[z + 1, :, :] = 0
                thresh_cyl[averaged_z, :, :] = new_slice
            elif np.count_nonzero(current_slice) < threshold and np.count_nonzero(next_slice) < threshold:
                thresh_cyl[z, :, :] = 0
                thresh_cyl[z + 1, :, :] = 0
            elif np.count_nonzero(current_slice) > threshold and np.count_nonzero(next_slice) < threshold:
                z_planes.append(z)
                z_slices.append(z)
        distance_between_planes = np.diff(z_planes)
        return (thresh_cyl, distance_between_planes, z_slices)

    def panalInsertMRAnalysis(self, sitkMRImage, userDefinedPoints):
        distances = []
        for panal_points_dir in userDefinedPoints:
            distance_dir = []
            for i in range(len(panal_points_dir) - 1):
                y1, x1 = panal_points_dir[i]
                y2, x2 = panal_points_dir[i + 1]
                distance = math.sqrt(
                    ((y2 - y1) * sitkMRImage.GetSpacing()[1]) ** 2 + ((x2 - x1) * sitkMRImage.GetSpacing()[0]) ** 2
                )
                distance_dir.append(distance)
            distances.append(distance_dir)
        return distances

    def panalInsertCTAnalysis(self, sitkImage, center_phantom):
        thresh_panal_sitk = af.ct_th_panalinsert(sitkImage, self.panalSliceCT, center_phantom, sitkImage.GetSpacing())
        thresh_panal = sitk.GetArrayFromImage(thresh_panal_sitk)
        z_thresh_panal = thresh_panal[self.panalSliceCT, :, :]  # axial slice with the panal segmentation

        # Distortion analysis
        _, x_panal = z_thresh_panal.shape
        y_ini = round(center_phantom[0])  # panal center y coordinate
        x_ini = round(center_phantom[1])  # panal center x coordinate

        # First, we check distances in the left-right direction
        pixel_right = []
        coord_right = []
        pixel_left = []
        coord_left = []
        # We move to the right
        for x in range(x_ini, x_panal, 1):
            pixel_right.append(z_thresh_panal[y_ini, x])
            coord_right.append([y_ini, x])
        # We move to the left
        for x in range(x_ini, -1, -1):
            pixel_left.append(z_thresh_panal[y_ini, x])
            coord_left.append([y_ini, x])

        # Finally, we check distances in the two diagonal direction
        # Secondary diagonal
        pixel_ds_up = []
        coord_ds_up = []
        pixel_ds_down = []
        coord_ds_down = []
        prev_x = x_ini
        prev_y = y_ini
        for _ in range(0, 125):
            new_x = round(prev_x + 2 * math.cos(math.radians(60)))
            new_y = round(prev_y - 2 * math.sin(math.radians(60)))
            pixel_ds_up.append(z_thresh_panal[new_y, new_x])
            coord_ds_up.append([new_y, new_x])
            prev_x = new_x
            prev_y = new_y

        prev_x = x_ini
        prev_y = y_ini
        for _ in range(0, 125):
            new_x = round(prev_x - 2 * math.cos(math.radians(60)))
            new_y = round(prev_y + 2 * math.sin(math.radians(60)))
            pixel_ds_down.append(z_thresh_panal[new_y, new_x])
            coord_ds_down.append([new_y, new_x])
            prev_x = new_x
            prev_y = new_y

        # Main diagonal
        pixel_dm_up = []
        coord_dm_up = []
        pixel_dm_down = []
        coord_dm_down = []
        prev_x = x_ini
        prev_y = y_ini
        for _ in range(0, 125):
            new_x = round(prev_x + 2 * math.cos(math.radians(60)))
            new_y = round(prev_y + 2 * math.sin(math.radians(60)))
            pixel_dm_down.append(z_thresh_panal[new_y, new_x])
            coord_dm_down.append([new_y, new_x])
            prev_x = new_x
            prev_y = new_y

        prev_x = x_ini
        prev_y = y_ini
        for _ in range(0, 125):
            new_x = round(prev_x - 2 * math.cos(math.radians(60)))
            new_y = round(prev_y - 2 * math.sin(math.radians(60)))
            pixel_dm_up.append(z_thresh_panal[new_y, new_x])
            coord_dm_up.append([new_y, new_x])
            prev_x = new_x
            prev_y = new_y

        # Once we have the lists with the pixel values in different directions, we obtain the distances of the
        # walls from the center of the panal in each direction (in mm)
        # Distances in list position between walls and panal origin
        wall_coord_r = af.panal_walls(pixel_right, coord_right)
        distances_coord_r_mm = af.calculate_mm_distances(
            wall_coord_r, sitkImage.GetSpacing()[1], sitkImage.GetSpacing()[0]
        )

        wall_coord_l = af.panal_walls(pixel_left, coord_left)
        distances_coord_l_mm = af.calculate_mm_distances(
            wall_coord_l, sitkImage.GetSpacing()[1], sitkImage.GetSpacing()[0]
        )

        wall_coord_dsu = af.panal_walls(pixel_ds_up, coord_ds_up)
        distances_coord_dsu_mm = af.calculate_mm_distances(
            wall_coord_dsu, sitkImage.GetSpacing()[1], sitkImage.GetSpacing()[0]
        )

        wall_coord_dsd = af.panal_walls(pixel_ds_down, coord_ds_down)
        distances_coord_dsd_mm = af.calculate_mm_distances(
            wall_coord_dsd, sitkImage.GetSpacing()[1], sitkImage.GetSpacing()[0]
        )

        wall_coord_dmu = af.panal_walls(pixel_dm_up, coord_dm_up)
        distances_coord_dmu_mm = af.calculate_mm_distances(
            wall_coord_dmu, sitkImage.GetSpacing()[1], sitkImage.GetSpacing()[0]
        )

        wall_coord_dmd = af.panal_walls(pixel_dm_down, coord_dm_down)
        distances_coord_dmd_mm = af.calculate_mm_distances(
            wall_coord_dmd, sitkImage.GetSpacing()[1], sitkImage.GetSpacing()[0]
        )

        return (
            distances_coord_r_mm,
            distances_coord_l_mm,
            distances_coord_dsu_mm,
            distances_coord_dsd_mm,
            distances_coord_dmu_mm,
            distances_coord_dmd_mm,
        )

    def runT2Analysis(self, sitkMRImage):
        center_phantom = af.phantom_center(sitkMRImage, self.panalSliceT2, sitkMRImage.GetSpacing()[2], "MR")
        thresh_cyl_sitk = af.th_cylinsert(
            sitkMRImage,
            self.startSliceT2,
            self.endSliceT2,
            sitkMRImage.GetSpacing(),
            center_phantom,
            True,
            "T2",
        )
        thresh_cyl = sitk.GetArrayFromImage(thresh_cyl_sitk)

        thresh_cyl, distance_between_planes, z_slices = self.cylindricalInsertAnalysis(thresh_cyl, 50)

        thresh_cyl_sitk = sitk.GetImageFromArray(thresh_cyl)
        thresh_cyl_sitk.CopyInformation(sitkMRImage)

        axial_distance = distance_between_planes * sitkMRImage.GetSpacing()[2]

        meshSitk = self.getMeshVertices(thresh_cyl_sitk, z_slices)

        meshResult = self.meshAnalysisMR(meshSitk, z_slices)

        panalInsertResult = self.panalInsertMRAnalysis(sitkMRImage, self.pointsT2)

        return (axial_distance, meshResult, panalInsertResult)

    def runT1Analysis(self, sitkMRImage):
        center_phantom = af.phantom_center(sitkMRImage, self.panalSliceT1, sitkMRImage.GetSpacing()[2], "MR")
        thresh_cyl_sitk = af.th_cylinsert(
            sitkMRImage,
            self.startSliceT1,
            self.endSliceT1,
            sitkMRImage.GetSpacing(),
            center_phantom,
            True,
            "T1",
        )
        thresh_cyl = sitk.GetArrayFromImage(thresh_cyl_sitk)

        thresh_cyl, distance_between_planes, z_slices = self.cylindricalInsertAnalysis(thresh_cyl, 50)

        thresh_cyl_sitk = sitk.GetImageFromArray(thresh_cyl)
        thresh_cyl_sitk.CopyInformation(sitkMRImage)

        axial_distance = distance_between_planes * sitkMRImage.GetSpacing()[2]

        meshSitk = self.getMeshVertices(thresh_cyl_sitk, z_slices, 2)

        meshResult = self.meshAnalysisMR(meshSitk, z_slices)

        panalInsertResult = self.panalInsertMRAnalysis(sitkMRImage, self.pointsT1)

        return (axial_distance, meshResult, panalInsertResult)

    def runCTAnalysis(self, sitkCtImage):
        center_phantom = af.phantom_center(sitkCtImage, self.panalSliceCT, sitkCtImage.GetSpacing()[2], "CT")
        thresh_cyl_sitk = af.th_cylinsert(
            sitkCtImage,
            self.startSliceCT,
            self.endSliceCT,
            sitkCtImage.GetSpacing(),
            center_phantom,
            self.insertFilled,
            "CT",
        )
        thresh_cyl = sitk.GetArrayFromImage(thresh_cyl_sitk)

        thresh_cyl, distance_between_planes, z_slices = self.cylindricalInsertAnalysis(thresh_cyl, 300)

        thresh_cyl_sitk = sitk.GetImageFromArray(thresh_cyl)
        thresh_cyl_sitk.CopyInformation(sitkCtImage)

        axial_distance = distance_between_planes * sitkCtImage.GetSpacing()[2]

        meshSitk = self.getMeshVertices(thresh_cyl_sitk, z_slices)
        meshResult = self.meshAnalysisCT(meshSitk, z_slices)

        panalInsertResult = self.panalInsertCTAnalysis(sitkCtImage, center_phantom)

        return axial_distance, meshResult, panalInsertResult

    def run(self):
        ctResult = t1Result = t2Result = None
        self.increaseStatusUpdateFunction(10, "Distortion")
        if self.ctVolume:
            sitkCtImage = sitkUtils.PullVolumeFromSlicer(self.ctVolume)
            ctResult = self.runCTAnalysis(sitkCtImage)
            self.increaseStatusUpdateFunction(40, "Distortion")
        if self.t1Volume:
            sitkT1Image = sitkUtils.PullVolumeFromSlicer(self.t1Volume)
            t1Result = self.runT1Analysis(sitkT1Image)
            self.increaseStatusUpdateFunction(60, "Distortion")
        if self.t2Volume:
            sitkT2Image = sitkUtils.PullVolumeFromSlicer(self.t2Volume)
            t2Result = self.runT2Analysis(sitkT2Image)
            self.increaseStatusUpdateFunction(80, "Distortion")
        self.callbackFunction(ctResult, t1Result, t2Result)
        self.increaseStatusUpdateFunction(100, "Distortion")


class RegistrationRoutine(SegmentationRoutine):
    def runImageAnalysis(self, sitkPetImage, sitkOtherImage, modality):
        voxelvol_pet = sitkPetImage.GetSpacing()[0] * sitkPetImage.GetSpacing()[1] * sitkPetImage.GetSpacing()[2]
        voi_pt, ref_voi_pt = af.pet_th(sitkPetImage, voxelvol_pet)
        threshold = af.ct_mr_th(sitkOtherImage, modality)
        ref_voi = af.resampling(sitkOtherImage, ref_voi_pt)

        # Intersection Threshold_CT/Threshold_MR and Ref_VOI_CT/Ref_VOI_MR: VOI_CT/VOI_MR
        th = sitk.GetArrayFromImage(threshold)
        th = th.astype(bool)
        rf_voi = sitk.GetArrayFromImage(ref_voi)
        rf_voi = rf_voi.astype(bool)

        th[~rf_voi] = 0  # intersection
        th = th.astype(int)

        voi = sitk.GetImageFromArray(th)
        voi.SetDirection(threshold.GetDirection())
        voi.SetSpacing(threshold.GetSpacing())
        voi.SetOrigin(threshold.GetOrigin())

        # Intersection VOI_PET and VOI_CT_resampled/VOI_MR_resampled
        # Resample VOI_CT/VOI_MR to PET image dimensions
        voi_resamp = af.resampling(sitkPetImage, voi)

        # Intersection
        voi_pet = sitk.GetArrayFromImage(voi_pt)
        voi_resampled = sitk.GetArrayFromImage(voi_resamp)
        intersection = voi_pet * voi_resampled  # intersection

        # DSC calculation #
        # Slices number calculation for DSC analysis
        li_pet = 0
        li_other = 0
        lf_pet = len(voi_pet[:, 0, 0]) - 1
        lf_other = len(voi_resampled[:, 0, 0]) - 1

        # Choose initial and final slice for PET image
        for s in range(0, len(voi_pet[:, 0, 0])):
            if np.mean(voi_pet[s, :, :]) == 0:
                li_pet = li_pet + 1
            else:
                break
        for s in range(li_pet, len(voi_pet[:, 0, 0])):
            if np.mean(voi_pet[s, :, :]) != 0:
                lf_pet = s
            else:
                break

        # Choose initial and final slice for CT image
        for s in range(0, len(voi_resampled[:, 0, 0])):
            if np.mean(voi_resampled[s, :, :]) == 0:
                li_other = li_other + 1
            else:
                break

        for s in range(li_other, len(voi_resampled[:, 0, 0])):
            if np.mean(voi_resampled[s, :, :]) != 0:
                lf_other = s
            else:
                break

        # Choose inferior and superior slice for DSC calculation
        if li_pet < li_other:
            li = li_pet
        else:
            li = li_other
        if lf_pet > li_other:
            lf = lf_pet
        else:
            lf = lf_other

        # DSC calculation
        diceCoefficients = {}
        for j in range(li, lf + 1):
            slice = voi_resampled[j, :, :]
            pet_slice = voi_pet[j, :, :]
            intersect_slice = intersection[j, :, :]
            area_ct = np.sum(slice == 1)
            area_pet = np.sum(pet_slice == 1)
            area_intersect = np.sum(intersect_slice == 1)

            dsccoeff = (2 * area_intersect) / (area_ct + area_pet)
            diceCoefficients[j] = dsccoeff

        # Get VOI_PET and VOI_CT_resampled/VOI_MR_resampled to calculate distance between CT/MR and PET
        # Load position of origin voxel
        PixelPosition_other = voi_resamp.GetOrigin()
        PixelPosition_PET = voi_pt.GetOrigin()
        PixelSpacing_other = voi_resamp.GetSpacing()
        PixelSpacing_PET = voi_pt.GetSpacing()

        # Difference between first VOI_CT/VOI_MR and VOI_PET slices in mm
        pos_other = (
            PixelPosition_other[2] + PixelSpacing_other[2] * li_other
        )  # Position in mm of the first VOI_CT_resampled/VOI_MR_resampled slice
        pos_PET = PixelPosition_PET[2] + PixelSpacing_PET[2] * li_pet  # Position in mm of the first VOI_PET slice
        DIF_CORTES = pos_PET - pos_other

        return diceCoefficients, DIF_CORTES

    def run(self):
        self.increaseStatusUpdateFunction(10, "Registration")
        if self.petVolume:
            ctResult = t1Result = t2Result = None
            sitkPetImage = sitkUtils.PullVolumeFromSlicer(self.petVolume)
            if self.ctVolume:
                sitkCTImage = sitkUtils.PullVolumeFromSlicer(self.ctVolume)
                ctResult = self.runImageAnalysis(sitkPetImage, sitkCTImage, "CT")

            self.increaseStatusUpdateFunction(60, "Registration")
            if self.t1Volume:
                sitkT1Image = sitkUtils.PullVolumeFromSlicer(self.t1Volume)
                t1Result = self.runImageAnalysis(sitkPetImage, sitkT1Image, "T1")
            self.increaseStatusUpdateFunction(80, "Registration")
            if self.t2Volume:
                sitkT2Image = sitkUtils.PullVolumeFromSlicer(self.t2Volume)
                t2Result = self.runImageAnalysis(sitkPetImage, sitkT2Image, "T2")
            self.increaseStatusUpdateFunction(90, "Registration")
            self.callbackFunction(ctResult, t1Result, t2Result)
        self.increaseStatusUpdateFunction(100, "Registration")


class ResolutionRoutine(SegmentationRoutine):
    def __init__(
        self,
        ctVolume,
        petVolume,
        t1Volume,
        t2Volume,
        resolutionTrianglesCT,
        resolutionTrianglesPET,
        resolutionTrianglesT1,
        resolutionTrianglesT2,
        callBackFunction,
        statusUpdateFunction,
    ):
        SegmentationRoutine.__init__(
            self, ctVolume, petVolume, t1Volume, t2Volume, callBackFunction, statusUpdateFunction
        )
        self.startCoordCT, self.resolutionTrianglesCT = resolutionTrianglesCT
        self.startCoordPET, self.resolutionTrianglesPET = resolutionTrianglesPET
        self.startCoordT1, self.resolutionTrianglesT1 = resolutionTrianglesT1
        self.startCoordT2, self.resolutionTrianglesT2 = resolutionTrianglesT2

    def runPETAnalysis(self, sitkPETImage):
        (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        ) = af.resol_pet(
            sitkPETImage,
            self.startCoordPET,
            self.resolutionTrianglesPET,
            sitkPETImage.GetSpacing(),
        )
        finalSegmentation = self.getSitkImageFromNumpyArray(finalSegmentation, sitkPETImage)
        return (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        )

    def runCTAnalysis(self, sitkCtImage):
        (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        ) = af.resol_ct(
            sitkCtImage,
            self.startCoordCT,
            self.resolutionTrianglesCT,
            sitkCtImage.GetSpacing(),
        )
        finalSegmentation = self.getSitkImageFromNumpyArray(finalSegmentation, sitkCtImage)
        return (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        )

    def runT1Analysis(self, sitkMrImage):
        (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        ) = af.resol_mr_t1(
            sitkMrImage,
            self.startCoordT1,
            self.resolutionTrianglesT1,
            sitkMrImage.GetSpacing(),
        )
        finalSegmentation = self.getSitkImageFromNumpyArray(finalSegmentation, sitkMrImage)
        return (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        )

    def runT2Analysis(self, sitkMrImage):
        (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        ) = af.resol_mr_t2(
            sitkMrImage,
            self.startCoordT2,
            self.resolutionTrianglesT2,
        )
        finalSegmentation = self.getSitkImageFromNumpyArray(finalSegmentation, sitkMrImage)
        return (
            contrast_total,
            holes_total,
            rc_holes_total,
            section_names,
            normalized_profiles,
            startingCoords,
            newEndCoords,
            plotImageValuesList,
            plotThreshHolesValuesList,
            finalSegmentation,
            resultValuesForPlt2,
        )

    def run(self):
        resultCT = resultPET = resultT1 = resultT2 = None
        self.increaseStatusUpdateFunction(10, "Resolution")
        if self.petVolume and self.startCoordPET is not None:
            sitkPetImage = sitkUtils.PullVolumeFromSlicer(self.petVolume)
            resultPET = self.runPETAnalysis(sitkPetImage)
        self.increaseStatusUpdateFunction(40, "Resolution")
        if self.ctVolume and self.startCoordCT is not None:
            sitkCtImage = sitkUtils.PullVolumeFromSlicer(self.ctVolume)
            resultCT = self.runCTAnalysis(sitkCtImage)
        self.increaseStatusUpdateFunction(60, "Resolution")
        if self.t1Volume and self.startCoordT1 is not None:
            sitkT1Image = sitkUtils.PullVolumeFromSlicer(self.t1Volume)
            resultT1 = self.runT1Analysis(sitkT1Image)
        self.increaseStatusUpdateFunction(80, "Resolution")
        if self.t2Volume and self.startCoordT2 is not None:
            sitkT2Image = sitkUtils.PullVolumeFromSlicer(self.t2Volume)
            resultT2 = self.runT2Analysis(sitkT2Image)
        self.increaseStatusUpdateFunction(90, "Resolution")
        self.callbackFunction(resultCT, resultPET, resultT1, resultT2)
        self.increaseStatusUpdateFunction(100, "Resolution")


class RadiomicsRoutine(SegmentationRoutine):
    def __init__(
        self,
        ctVolume,
        petVolume,
        t1Volume,
        t2Volume,
        ctSegment,
        petSegment,
        t1Segment,
        t2Segment,
        callBackFunction=None,
        statusUpdateFunction=None,
        resourcePath=None,
    ):
        super().__init__(ctVolume, petVolume, t1Volume, t2Volume, callBackFunction, statusUpdateFunction)
        self.ctSegment = ctSegment
        self.petSegment = petSegment
        self.t1Segment = t1Segment
        self.t2Segment = t2Segment
        self.resourcePath = resourcePath

    def getBinWidht(self, sitkImage):
        image_array = sitk.GetArrayFromImage(sitkImage)
        unique_values = np.unique(image_array)
        data_range = np.max(unique_values) - np.min(unique_values)
        num_bins = 128
        bin_width = int(data_range / num_bins)
        return bin_width

    def runRadiomicsAnalysis(self, slicerImage, slicerSegmentation, modality):
        sitkImage = sitkUtils.PullVolumeFromSlicer(slicerImage)
        all_radiomics = {}
        segmentationsForRadiomics = {}
        if slicerSegmentation is not None:
            nameContourMapping = self.convertSlicerSegmentationToSitkSegmentations(slicerSegmentation, slicerImage)
            for contourName in nameContourMapping.keys():
                sitkContour = nameContourMapping[contourName]
                segm_radiomics = af.get_radiomics(
                    sitkImage,
                    sitkContour,
                    modality,
                    self.getBinWidht(sitkImage),
                    os.path.join(self.resourcePath, "Params.yaml"),
                )
                all_radiomics[contourName] = segm_radiomics
        else:
            referenceImageFileName = os.path.join(self.resourcePath, modality, f"{modality}.nrrd")
            im_ref = sitk.ReadImage(referenceImageFileName, sitk.sitkFloat32)  # reference image
            im_size = sitkImage.GetSize()[2]  # axial length of the analysis image
            ref_im_size = im_ref.GetSize()[2]

            if ref_im_size > im_size:
                # Reduce the reference image to the length of the analysis image
                im_ref = sitk.Extract(im_ref, [im_ref.GetSize()[0], im_ref.GetSize()[1], im_size], [0, 0, 0])
            registration_transform = af.register(sitk.Cast(sitkImage, sitk.sitkFloat32), im_ref)

            seg_ref_paths = [
                os.path.join(self.resourcePath, modality, d)
                for d in os.listdir(os.path.join(self.resourcePath, modality))
                if d != f"{modality}.nrrd" and d[-5:] == ".nrrd"
            ]

            for seg_ref in seg_ref_paths:
                segm_ref = sitk.ReadImage(seg_ref, sitk.sitkFloat32)  # reference segmentation
                if ref_im_size > im_size:
                    segm_ref = sitk.Extract(
                        segm_ref, [segm_ref.GetSize()[0], segm_ref.GetSize()[1], im_size], [0, 0, 0]
                    )  # Reduce the reference segmentation to the length of the analysis image
                name_segm_ref = os.path.splitext(os.path.basename(seg_ref))[0]

                moving_seg_resampled = sitk.Resample(
                    segm_ref, sitkImage, registration_transform, sitk.sitkLinear, 0.0, segm_ref.GetPixelID()
                )
                moving_seg_resampled_int = sitk.Cast(moving_seg_resampled, sitk.sitkUInt8)

                segmentationsForRadiomics[name_segm_ref] = moving_seg_resampled_int

                segm_radiomics = af.get_radiomics(
                    sitkImage,
                    moving_seg_resampled_int,
                    modality,
                    self.getBinWidht(sitkImage),
                    os.path.join(self.resourcePath, "Params.yaml"),
                )

                all_radiomics[name_segm_ref] = segm_radiomics

        return all_radiomics, segmentationsForRadiomics

    def convertSlicerSegmentationToSitkSegmentations(self, segmentationNode, referenceVolume):
        sitkSegmentsWithNames = {}
        if segmentationNode is not None:
            numberOfSegments = segmentationNode.GetSegmentation().GetNumberOfSegments()
            for i in range(0, numberOfSegments):
                currentSegment = segmentationNode.GetSegmentation().GetNthSegment(i)
                segmentName = currentSegment.GetName()
                segmentId = segmentationNode.GetSegmentation().GetSegmentIdBySegment(currentSegment)
                segmentIds = vtk.vtkStringArray()
                segmentIds.InsertNextValue(segmentId)
                labelmapVolumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode")
                slicer.vtkSlicerSegmentationsModuleLogic.ExportSegmentsToLabelmapNode(
                    segmentationNode, segmentIds, labelmapVolumeNode, referenceVolume
                )
                sitkSegmentation = sitkUtils.PullVolumeFromSlicer(labelmapVolumeNode)
                sitkSegmentsWithNames[segmentName] = sitkSegmentation
                slicer.mrmlScene.RemoveNode(labelmapVolumeNode)

        return sitkSegmentsWithNames

    def run(self):
        radiomicValues = {}
        radiomicSegmentations = {}
        self.increaseStatusUpdateFunction(10, "Radiomics")
        if self.petVolume is not None:
            radiomicValues["PET"], radiomicSegmentations["PET"] = self.runRadiomicsAnalysis(
                self.petVolume, self.petSegment, "PET"
            )
            self.increaseStatusUpdateFunction(30, "Radiomics")
        if self.ctVolume is not None:
            radiomicValues["CT"], radiomicSegmentations["CT"] = self.runRadiomicsAnalysis(
                self.ctVolume, self.ctSegment, "CT"
            )
            self.increaseStatusUpdateFunction(50, "Radiomics")
        if self.t1Volume is not None:
            radiomicValues["T1"], radiomicSegmentations["T1"] = self.runRadiomicsAnalysis(
                self.t1Volume, self.t1Segment, "T1"
            )
            self.increaseStatusUpdateFunction(70, "Radiomics")
        if self.t2Volume is not None:
            radiomicValues["T2"], radiomicSegmentations["T2"] = self.runRadiomicsAnalysis(
                self.t2Volume, self.t2Segment, "T2"
            )
            self.increaseStatusUpdateFunction(90, "Radiomics")
        self.callbackFunction(radiomicValues, radiomicSegmentations)
        self.increaseStatusUpdateFunction(100, "Radiomics")


class QuantificationRoutine(SegmentationRoutine):
    def __init__(
        self,
        ctVolume,
        petVolume,
        t1Volume,
        t2Volume,
        activityDecay,
        insertType,
        quantificationStartSlicePET,
        quantificationStartSliceCT,
        quantificationStartSliceT1,
        quantificationStartSliceT2,
        callBackFunction=None,
        statusUpdateFunction=None,
    ):
        super().__init__(ctVolume, petVolume, t1Volume, t2Volume, callBackFunction, statusUpdateFunction)
        self.activityDecay = activityDecay
        self.insertType = insertType
        self.quantificationStartSlicePET = quantificationStartSlicePET
        self.quantificationStartSliceCT = quantificationStartSliceCT
        self.quantificationStartSliceT1 = quantificationStartSliceT1
        self.quantificationStartSliceT2 = quantificationStartSliceT2

    @staticmethod
    def _getCenterSliceForVolume(imgVolume, innerMethod, postProcess):
        sitkPetVolume = sitkUtils.PullVolumeFromSlicer(imgVolume)
        sitkPetVolume = sitk.SmoothingRecursiveGaussian(sitkPetVolume, 2.0)
        size = list(sitkPetVolume.GetSize())
        size[2] = 0
        minMaxFilter = sitk.MinimumMaximumImageFilter()

        extractor = sitk.ExtractImageFilter()
        extractor.SetSize(size)
        diffTo7ConnComponents = []
        avgSliceRoundness = []
        sliceMaxima = []
        for i in range(sitkPetVolume.GetSize()[2]):
            index = [0, 0, i]
            extractor.SetIndex(index)
            sitkPetSlice = extractor.Execute(sitkPetVolume)
            minMaxFilter.Execute(sitkPetSlice)
            sliceMaximum = minMaxFilter.GetMaximum()
            sliceMaxima.append(sliceMaximum)

            diffTo7ConnComponent, roundNess = innerMethod(sitkPetSlice, sliceMaximum)
            diffTo7ConnComponents.append(diffTo7ConnComponent)

            avgSliceRoundness.append(roundNess)

        diffTo7ConnComponents = np.asarray(diffTo7ConnComponents)
        sliceMaxima = np.asarray(sliceMaxima)
        avgSliceRoundness = np.asarray(avgSliceRoundness)

        return postProcess(diffTo7ConnComponents, sliceMaxima, avgSliceRoundness)

    ##TODO: move to logic class
    @staticmethod
    def getCenterSlicerForCtInsert(imgVolume):
        def thresholdAndParamExtraction(sitkSlice, upperThreshold):
            lowerThreshold = -615
            segmentedSphereCandidates = sitk.BinaryThreshold(sitkSlice, lowerThreshold, upperThreshold)
            connCompImage = sitk.ConnectedComponent(segmentedSphereCandidates)

            label_shape_filter = sitk.LabelShapeStatisticsImageFilter()
            label_shape_filter.Execute(connCompImage)
            numberOfLabels = 0
            roundNess = 0.0
            for label in range(1, label_shape_filter.GetNumberOfLabels() + 1):
                labelSize = label_shape_filter.GetPhysicalSize(label)
                if labelSize > 500 and labelSize < 1500:
                    if np.abs(1.0 - label_shape_filter.GetRoundness(label)) < 0.5:
                        numberOfLabels = numberOfLabels + 1
                        roundNess = roundNess + np.abs(1.0 - label_shape_filter.GetRoundness(label))
            if numberOfLabels > 6:
                differnceTo7QuantInserts = numberOfLabels - 7
            else:
                differnceTo7QuantInserts = 200
            roundNess = roundNess / np.max([1.0, numberOfLabels])

            return (differnceTo7QuantInserts, roundNess)

        def postProcessing(diffTo7ConnComponents, sliceMaxima, avgSliceRoundness):
            tmpMinSliceRoundness = np.min(avgSliceRoundness[diffTo7ConnComponents == diffTo7ConnComponents.min()])

            centerSlice = np.argwhere(
                (diffTo7ConnComponents == diffTo7ConnComponents.min()) & (avgSliceRoundness == tmpMinSliceRoundness)
            )[0][0]
            minDiff = diffTo7ConnComponents[centerSlice]
            sliceToTake = centerSlice
            while np.abs(diffTo7ConnComponents[centerSlice] - minDiff) < 4:
                sliceToTake = centerSlice
                centerSlice = centerSlice - 1
            return sliceToTake

        return QuantificationRoutine._getCenterSliceForVolume(imgVolume, thresholdAndParamExtraction, postProcessing)

    ##TODO: mode to logic class
    @staticmethod
    def getCenterSlicerForPetInsert(imgVolume):
        def thresholdAndParamExtraction(sitkSlice, upperThreshold):
            lowerThreshold = upperThreshold * 0.5
            segmentedSphereCandidates = sitk.BinaryThreshold(sitkSlice, lowerThreshold, upperThreshold)
            connCompImage = sitk.ConnectedComponent(segmentedSphereCandidates)

            label_shape_filter = sitk.LabelShapeStatisticsImageFilter()
            label_shape_filter.Execute(connCompImage)
            differnceTo7QuantInserts = np.abs(7 - label_shape_filter.GetNumberOfLabels())
            roundNess = 1.0
            if label_shape_filter.GetNumberOfLabels() > 0:
                roundNess = 0.0
                for label in range(1, label_shape_filter.GetNumberOfLabels() + 1):
                    roundNess = roundNess + np.abs(1.0 - label_shape_filter.GetRoundness(label))
                roundNess = roundNess / np.max([1.0, label_shape_filter.GetNumberOfLabels()])

            return (differnceTo7QuantInserts, roundNess)

        def postProcessing(diffTo7ConnComponents, sliceMaxima, avgSliceRoundness):
            diffTo7ConnComponents[sliceMaxima < np.max(sliceMaxima) * 0.4] = 7
            tmpMinSliceRoundness = np.min(avgSliceRoundness[diffTo7ConnComponents == diffTo7ConnComponents.min()])

            centerSlice = np.argwhere(
                (diffTo7ConnComponents == diffTo7ConnComponents.min()) & (avgSliceRoundness == tmpMinSliceRoundness)
            )[0][0]
            minDiff = diffTo7ConnComponents[centerSlice]
            sliceToTake = centerSlice
            while diffTo7ConnComponents[centerSlice] == minDiff:
                sliceToTake = centerSlice
                centerSlice = centerSlice - 1
            return sliceToTake

        return QuantificationRoutine._getCenterSliceForVolume(imgVolume, thresholdAndParamExtraction, postProcessing)

    def runCtAnalysis(self, sitkCTImage, initalCoord):
        image = sitk.GetArrayFromImage(sitkCTImage)

        spacing = sitkCTImage.GetSpacing()

        _, final_coord = af.segmentation_region(sitkCTImage, initalCoord, float(spacing[2]), "CT", self.insertType)

        names_ct, _, inserts_ct = af.inserts_ct(sitkCTImage, initalCoord, final_coord, spacing, self.insertType)

        intensity_ct = []
        err_intensity_ct = []
        # Go through each insert
        for i in range(len(inserts_ct)):
            # Create a reduced analysis region inside the insert for HU quantification
            segm_den, _ = af.reduced_segmentation_cyl(inserts_ct[i], float(sitkCTImage.GetSpacing()[0]))

            intensity, err = af.quantification(image, sitk.GetArrayFromImage(segm_den))
            intensity_ct.append(intensity)
            err_intensity_ct.append(err)

        return names_ct, intensity_ct, err_intensity_ct

    def runPetAnalysis(self, sitkPetImage, initalCoord):
        spacing = sitkPetImage.GetSpacing()
        image = sitk.GetArrayFromImage(sitkPetImage)
        area_pet, final_coord_pet = af.segmentation_region(sitkPetImage, initalCoord, float(spacing[2]), "PET", "RC")
        names_pet, _, inserts_pet = af.inserts_pet(sitkPetImage, area_pet, initalCoord, final_coord_pet, spacing)

        diameter_pet = []
        intensity_pet = []
        err_intensity_pet = []
        rc_ac = []
        for i in range(len(inserts_pet)):
            segm_activity, diameter_insert = af.reduced_segmentation_cyl(inserts_pet[i], float(spacing[0]))

            if not np.isnan(diameter_insert):
                diameter_pet.append(diameter_insert)

                in_pet, err_pet = af.quantification(image, sitk.GetArrayFromImage(segm_activity))
                intensity_pet.append(in_pet)
                err_intensity_pet.append(err_pet)

                rc = in_pet / self.activityDecay
                rc_ac.append(rc)

        return names_pet, diameter_pet, intensity_pet, err_intensity_pet, rc_ac

    def runMRAnalysis(self, sitkImage, initialCoord):
        _, final_coord = af.segmentation_region(
            sitkImage, initialCoord, float(sitkImage.GetSpacing()[2]), "MR", self.insertType
        )

        names_mr, _, inserts_mr = af.inserts_mr(sitkImage, initialCoord, final_coord, sitkImage.GetSpacing())

        intensity_mr = []
        err_intensity_mr = []
        # Go through each insert
        imageA = sitk.GetArrayFromImage(sitkImage)
        for i in range(len(inserts_mr)):
            # Create a reduced analysis region inside the insert for intensity quantification
            segm_mr, _ = af.reduced_segmentation_cyl(inserts_mr[i], float(sitkImage.GetSpacing()[0]))

            intensity, err = af.quantification(imageA, sitk.GetArrayFromImage(segm_mr))
            intensity_mr.append(intensity)
            err_intensity_mr.append(err)

        return names_mr, intensity_mr, err_intensity_mr

    def run(self):
        names_pet = diameter_pet = intensity_pet = err_intensity_pet = rc_ac = None
        self.increaseStatusUpdateFunction(10, "Quantification")
        if self.petVolume is not None and self.quantificationStartSlicePET > -1:
            sitkPetImage = sitkUtils.PullVolumeFromSlicer(self.petVolume)
            (names_pet, diameter_pet, intensity_pet, err_intensity_pet, rc_ac) = self.runPetAnalysis(
                sitkPetImage, self.quantificationStartSlicePET
            )

        self.increaseStatusUpdateFunction(20, "Quantification")

        names_ct = intensity_ct = err_intensity_ct = None
        if self.ctVolume is not None and self.quantificationStartSliceCT > -1:
            sitkCtImage = sitkUtils.PullVolumeFromSlicer(self.ctVolume)
            (names_ct, intensity_ct, err_intensity_ct) = self.runCtAnalysis(
                sitkCtImage, self.quantificationStartSliceCT
            )

        self.increaseStatusUpdateFunction(40, "Quantification")
        names_t1 = intensity_t1 = err_intensity_t1 = None
        if (
            self.t1Volume is not None
            and self.quantificationStartSliceT1 > -1
            and self.insertType == "Materials Characterization"
        ):
            sitkT1Image = sitkUtils.PullVolumeFromSlicer(self.t1Volume)
            (names_t1, intensity_t1, err_intensity_t1) = self.runMRAnalysis(
                sitkT1Image, self.quantificationStartSliceT1
            )

        self.increaseStatusUpdateFunction(60, "Quantification")
        names_t2 = intensity_t2 = err_intensity_t2 = None
        if (
            self.t2Volume is not None
            and self.quantificationStartSliceT2 > -1
            and self.insertType == "Materials Characterization"
        ):
            sitkT2Image = sitkUtils.PullVolumeFromSlicer(self.t2Volume)
            (names_t2, intensity_t2, err_intensity_t2) = self.runMRAnalysis(
                sitkT2Image, self.quantificationStartSliceT2
            )
        self.increaseStatusUpdateFunction(80, "Quantification")
        self.callbackFunction(
            names_ct,
            intensity_ct,
            err_intensity_ct,
            names_pet,
            diameter_pet,
            intensity_pet,
            err_intensity_pet,
            rc_ac,
            names_t1,
            intensity_t1,
            err_intensity_t1,
            names_t2,
            intensity_t2,
            err_intensity_t2,
        )
        self.increaseStatusUpdateFunction(100, "Quantification")


#
# QAHybridTest
#


class QAHybridTest(ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.
    Uses ScriptedLoadableModuleTest base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """Do whatever is needed to reset the state - typically a scene clear will be enough."""
        slicer.mrmlScene.Clear()

    def runTest(self):
        """Run as few or as many tests as needed here."""
        self.setUp()
        self.test_QAHybrid1()

    def test_QAHybrid1(self):
        """Ideally you should have several levels of tests.  At the lowest level
        tests should exercise the functionality of the logic with different inputs
        (both valid and invalid).  At higher levels your tests should emulate the
        way the user would interact with your code and confirm that it still works
        the way you intended.
        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of your
        module.  For example, if a developer removes a feature that you depend on,
        your test should break so they know that the feature is needed.
        """

        self.delayDisplay("Starting the test")

        # Get/create input data

        self.delayDisplay("Test passed")
