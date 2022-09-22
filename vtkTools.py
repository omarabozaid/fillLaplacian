""" 
This class is an API for some vtk functions
"""

import numpy as np
import vtk

from vtkmodules.vtkFiltersCore import (
    vtkFeatureEdges,
    vtkCleanPolyData,
    vtkFeatureEdges,
    vtkIdFilter,
    vtkTriangleFilter,
    vtkPolyDataNormals
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonComputationalGeometry import vtkParametricRandomHills
from vtkmodules.vtkCommonCore import (
    VTK_DOUBLE,
    vtkIdList,
    vtkLookupTable,
    vtkVersion
)
from vtkmodules.vtkCommonTransforms import vtkTransform
from numpy_support import *

class vtkTools():
    def __init__(self)->None:
        self.Inf:float=float(1e6)
        self.NegInf:float=float(-1e6)

    """ 
    Given polydata return sub-polydata 
    containing the feature edges
    """
    def extractPolyDataEdges(self,
                            polydata
                            ):
        featureEdges = vtkFeatureEdges()
        featureEdges.SetInputData(polydata)
        featureEdges.BoundaryEdgesOn()
        featureEdges.FeatureEdgesOff()
        featureEdges.ManifoldEdgesOff()
        featureEdges.NonManifoldEdgesOff()
        featureEdges.ColoringOn()
        featureEdges.Update()
        return featureEdges.GetOutput()

    """ 
    Given polydata/or vtk mesh
    return sub-polydata/ or vtk mesh
    clipped by vector and point
    """
    def vtkClip(self,
            polyData,
            axis:np.ndarray,
            origin:np.ndarray
            ):
        plane = vtk.vtkPlane()
        plane.SetNormal(axis[0],axis[1],axis[2])
        plane.SetOrigin(origin[0],origin[1],origin[2])
        clip = vtk.vtkClipPolyData()
        clip.SetClipFunction(plane)
        clip.SetInputData(polyData)
        clip.Update()
        return clip.GetOutput()

    """ 
    Given polydata return points as np array
    """
    def points(self,polyData)->np.ndarray:
        nodes = polyData.GetPoints().GetData()
        return vtk_to_numpy(nodes)


    """ 
    Given polydata of feature edges
    return two arrays
        1. segments: edge=tuple (vertex_index,vertex_index)
        2. list of points as np.array [p1,p2,p3....,pn]
    """
    def polyEdgesToNumpy(self,
                        polydata
                        ):
        nLines=polydata.GetNumberOfCells()
        lines = polydata.GetLines()
        linesData = lines.GetData()
        nCols = linesData.GetNumberOfValues()//nLines
        numpy_lines = vtk_to_numpy(linesData)
        numpy_lines = numpy_lines.reshape((-1,nCols))
        segments=np.full((nLines,2),[-1,-1],dtype="int")
        for i in range(nLines):
            segments[i]=np.asarray([numpy_lines[i][1],numpy_lines[i][2]])
        return np.asarray(segments)

    """ 
    Given polydata of stl
    return each triangle as list of 3 indices
    """
    def triangles(self,
                polydata
                ):
        cells = polydata.GetPolys()
        nCells = cells.GetNumberOfCells()
        array = cells.GetData()
        nCols = array.GetNumberOfValues()//nCells
        numpy_cells = vtk_to_numpy(array)
        numpy_cells = numpy_cells.reshape((-1,nCols))
        return numpy_cells[:,1:4]

    """ 
    Given 
        list of points
        list of triangles
        return list of centers
    """
    def centersTriangles(self,
                points:np.ndarray,
                triangles:np.ndarray
                ):
        centers=[None]*len(triangles)
        i:int=0
        for t in triangles:
            center=np.full(3,0.0,dtype="float")
            for j in range(3):
                center+=points[t[j]]
            centers[i]=center*(1.0/3.0)
            i+=1
        return centers

    """ 
    Given polydata of stl
    return each triangle as list of 3 indices
    """
    def cells(self,
            polydata
            ):
        cells = polydata.GetCells()
        nCells = polydata.GetNumberOfCells()
        array = cells.GetData()
        nCols = array.GetNumberOfValues()//nCells
        numpy_cells = vtk_to_numpy(array)
        numpy_cells = numpy_cells.reshape((-1,nCols))
        return numpy_cells[:,1:5]

    """ 
    Given polydata of stl
    return each triangle as list of 3 indices
    """
    def cellsP(self,
            polydata
            ):
        cells = polydata.GetPolys()
        nCells = polydata.GetNumberOfCells()
        array = cells.GetData()
        nCols = array.GetNumberOfValues()//nCells
        numpy_cells = vtk_to_numpy(array)
        numpy_cells = numpy_cells.reshape((-1,nCols))
        return numpy_cells[:,1:5]

    """ 
    Given polydata
    return generic cells
    """
    def cells2(self,
            polydata
            ):
        cells = polydata.GetPolys()
        array = cells.GetData()
        numpy_cells = vtk_to_numpy(array)
        cells=[]
        i:int=0
        while i<len(numpy_cells):
            nPts=numpy_cells[i]
            cell=[]
            for j in range(1,nPts+1):
                cell.append(numpy_cells[i+j])
            if(len(cell)>2):
                cells.append(np.asarray(cell))
            i+=(nPts+1)
        return cells

    """ 
    Given array of edges
    compute the angle w.r.t an axis
    extract only polylines (edges) 
    below some tolerance
    """
    def npExtractPlaneEdges(self,
                            axis:np.ndarray,
                            segments:np.ndarray,
                            points:np.ndarray,
                            tol:float=1e-1
                            )->np.ndarray:
        extractedSegmentsIndices=[]
        for i in range(len(segments)):
            segment=segments[i]
            pStartIndex=segment[0]
            pEndIndex=segment[1]
            vector=points[pStartIndex]-points[pEndIndex]
            vector/=np.linalg.norm(vector)
            projection=abs(np.dot(vector, axis))
            if(projection<tol):
                extractedSegmentsIndices.append(i)
        return np.asarray(extractedSegmentsIndices)

    """ 
    Compute center of each edge
    given np array of edges 
    """
    def edgesCenters(self,
                    segments:np.ndarray,
                    points:np.ndarray
                    )->np.ndarray:
        nEdges=len(segments)
        centers=np.full(
            (nEdges,3),[-self.Inf,-self.Inf,-self.Inf],dtype="float"
            )
        i=0
        for segment in segments:
            centers[i]=0.5*(points[segment[0]]+points[segment[0]])
            i+=1
        return centers    

    """
    Compute curvature
    return it as numpy array
    """
    def curvature(self,polydata,curvature_type:str="Mean"):
        curvature = vtk.vtkCurvatures()
        curvature.SetInputData(polydata)
        if curvature_type=="Mean":
            curvature.SetCurvatureTypeToMean()
        elif curvature_type=="Gauss":
            curvature.SetCurvatureTypeToGauss()
        curvature.Update()
        return curvature.GetOutput()
    
    """
    Return decimated polydata
    """
    def decimation(self,polydata):
        decimate = vtk.vtkDecimatePro()
        decimate.SetInputData(polydata)
        decimate.SetTargetReduction(.10)
        decimate.Update()
        return decimate.GetOutput()

    """
    Return np array as field
    """
    def npField(self,
        vtkArr
        )->np.ndarray:
        return vtk_to_numpy(vtkArr)
        
    """
    Convert graph edges into polydata
    """
    def npEdgesToPolyData(self,
        npEdges:np.ndarray,
        npPoints:np.ndarray
        ):

        nPoints=len(npPoints)
        nLines=len(npEdges)
        points = vtk.vtkPoints()
        for npPoint in npPoints:
            points.InsertNextPoint(npPoint)

        lines = vtk.vtkCellArray()
        for i in range(nLines):
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, npEdges[i][0])
            line.GetPointIds().SetId(1, npEdges[i][1])
            lines.InsertNextCell(line)
        linesPolyData = vtk.vtkPolyData()
        linesPolyData.SetPoints(points)
        linesPolyData.SetLines(lines)
        return linesPolyData

    """
    Convert points into polydata
    """
    def npPointsToPolyData(self,
        npPoints:np.ndarray
        ):
        points = vtk.vtkPoints()
        pd = vtk.vtkPolyData()
        vertices = vtk.vtkCellArray()
        for pt in npPoints:
            id=points.InsertNextPoint(pt)
            vertices.InsertNextCell(1)
            vertices.InsertCellPoint(id)
        pd.SetPoints(points)
        pd.SetVerts(vertices)

        return pd

    """
    convert numpy array to a vtk array
    """
    def npToVtk(self,
        arr:np.ndarray,
        name:str,
        grid,
        nodal:bool=False,
        cell:bool=True
        ):

        vtkField = numpy_to_vtk(
            num_array=arr,
            deep=True,
            array_type=VTK_DOUBLE
            )
        vtkField.SetName(name)
        if nodal:
            grid.GetPointData().AddArray(vtkField)
            grid.GetPointData().SetActiveScalars(name)
        elif cell:
            grid.GetCellData().AddArray(vtkField)
            grid.GetCellData().SetActiveScalars(name)
        return grid

    """
    convert numpy array to a vtk array
    """
    def npToVtkMultiBlock(self,
        arr:np.ndarray,
        name:str,
        grid,
        nodal:bool=False,
        cell:bool=True
        ):

        nBlocks=len(arr)
        for i in range(nBlocks):
            vtkField = numpy_to_vtk(
                num_array=arr[i],
                deep=True,
                array_type=VTK_DOUBLE
                )
            vtkField.SetName(name)
            if nodal:
                grid.GetBlock(i).GetPointData().AddArray(vtkField)
                grid.GetBlock(i).GetPointData().SetActiveScalars(name)
            if cell:
                grid.GetBlock(i).GetCellData().AddArray(vtkField)
                grid.GetBlock(i).GetCellData().SetActiveScalars(name)

    """
    return unit normal vector field
    of polydata
    """
    def npUnitNormals(self,
        polydata
        ):
        normalsFilter=vtkPolyDataNormals()
        normalsFilter.SetInputData(polydata)
        normalsFilter.ComputePointNormalsOff()
        normalsFilter.ComputeCellNormalsOn()
        normalsFilter.Update()
        normalsVtkArray = normalsFilter.GetOutput().GetCellData().GetVectors('Normals')
        return self.npField(normalsVtkArray)
