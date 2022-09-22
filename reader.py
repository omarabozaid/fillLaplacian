from vtkmodules.vtkIOGeometry import vtkSTLReader,vtkSTLWriter
import vtk

class reader:
    def __init__(self,file:str)-> None:
        self.filename=file
        if self.filename.endswith('.stl'):
            self.readSTL()
        elif self.filename.endswith('.igs'):
            self.readIges()
        elif self.filename.endswith('.vtu'):
            self.readVTU()
        elif self.filename.endswith('.vtp'):
            self.readVTP()
        elif self.filename.endswith('.vtk'):
            self.readVTK()

    def readSTL(self)-> None:
        self.reader=vtkSTLReader()
        self.reader.SetFileName(self.filename)
        self.reader.Update()
        self.polydata=self.reader.GetOutput()

    def readVTU(self)-> None:
        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.reader.SetFileName(self.filename)
        self.reader.Update()  
        self.polydata = self.reader.GetOutput()
    
    def readVTK(self)-> None:
        self.reader = vtk.vtkPolyDataReader()
        self.reader.SetFileName(self.filename)
        self.reader.Update()  
        self.polydata = self.reader.GetOutput()
    
    def readVTP(self)-> None:
        self.reader = vtk.vtkXMLPolyDataReader()
        self.reader.SetFileName(self.filename)
        self.reader.Update()  
        self.polydata = self.reader.GetOutput()

    def polyData(self):
        return self.polydata
    
    def GetOutput(self):
        return self.reader.GetOutput()