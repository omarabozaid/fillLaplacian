import geometry as gm
import vtk
import copy
import numpy as np
import math as math
import gc
import vtkTools as vtkTools

class surfaceMesh:
    
    def __init__(self,
        nodesList,
        facesList,
        buildTopology:bool=True,
        name:str="patch",
        connectivityLists=[None]
    ) -> None:
        self.nodesList=np.asarray(nodesList)
        self.facesList=facesList
        self.name=name
        self.leftCorner=np.full(3,1e6,dtype="float")
        self.rightCorner=np.full(3,-1e6,dtype="float")
        self.centerOfMass=np.full(3,0.0,dtype="float")
        for pt in self.nodesList:
            self.leftCorner[0]=min(self.leftCorner[0],pt[0])
            self.leftCorner[1]=min(self.leftCorner[1],pt[1])
            self.leftCorner[2]=min(self.leftCorner[2],pt[2])
            self.rightCorner[0]=max(self.rightCorner[0],pt[0])
            self.rightCorner[1]=max(self.rightCorner[1],pt[1])
            self.rightCorner[2]=max(self.rightCorner[2],pt[2])
            self.centerOfMass+=pt
        self.centerOfMass/=len(self.nodesList)
        if buildTopology:
            self.__buildTopology()  
        else:
            self.nodesEdgesIndices = connectivityLists[0]
            self.nodesNodesList = connectivityLists[1]
            self.edgesList = connectivityLists[2]
            self.edgesFaces = connectivityLists[3]
            self.isBoundaryEdge = connectivityLists[4]
            self.isBoundaryNode = connectivityLists[5]
            self.facesEdgesIndices = connectivityLists[6]
            self.elementTypes = connectivityLists[7]
            self.centersList = connectivityLists[8]
    
    def topologyInformation(self):
        return (
                [
                    self.nodesEdgesIndices, 
                    self.nodesNodesList, 
                    self.edgesList ,
                    self.edgesFaces ,
                    self.isBoundaryEdge,
                    self.isBoundaryNode,
                    self.facesEdgesIndices,
                    self.elementTypes,
                    self.centersList
                ]
        )
        
    def __buildTopology(self)->None:

        nodesEdgesList=[None]*len(self.nodesList)
        nodesEdgesIndices=[None]*len(self.nodesList)
        nodesNodesList=[None]*len(self.nodesList)
        edgesList=[]
        edgesFaces=[]
        isBoundaryEdge=[]
        isBoundaryNode=[False]*len(self.nodesList)
        facesEdgesIndices=[None]*len(self.facesList)
        elementTypes=[None]*len(self.facesList)
        centersList=[None]*len(self.facesList)
        for i in range(len(self.facesList)):
            if len(self.facesList[i])==3:
                elementTypes[i]="triangle"
            elif len(self.facesList[i])==4:
                elementTypes[i]="quad"
            else:
                elementTypes[i]="polygon"
            c=np.full(3,0.0,dtype="float")
            for n in self.facesList[i]:
                c+=self.nodesList[n]
            c/=len(self.facesList[i])
            centersList[i]=c
            faceEdges=[]
            for j in range(len(self.facesList[i])-1):
                v1=self.facesList[i][j]
                v2=self.facesList[i][j+1]
                if v2<v1:
                    v1,v2=v2,v1
                faceEdges.append([v1,v2])
            v1=self.facesList[i][-1]
            v2=self.facesList[i][0]
            if v2<v1:
                v1,v2=v2,v1
            faceEdges.append([v1,v2])
            for e in faceEdges:
                startIdx:int=e[0]
                endIdx:int=e[1]
                if nodesEdgesList[startIdx] == None:
                    edgeIndex=len(edgesList)
                    nodesEdgesList[startIdx]=[e]
                    nodesEdgesIndices[startIdx]=[edgeIndex]
                    if nodesEdgesList[endIdx] == None:
                        nodesEdgesList[endIdx]=[e]
                        nodesEdgesIndices[endIdx]=[edgeIndex]
                    else:
                        nodesEdgesList[endIdx].append(e)
                        nodesEdgesIndices[endIdx].append(edgeIndex)
                    edgesList.append(e)
                    edgesFaces.append([i])
                    if facesEdgesIndices[i]== None:
                        facesEdgesIndices[i]=[edgeIndex]
                    else:
                        facesEdgesIndices[i].append(edgeIndex)
                else:
                    isRepeated:bool=False
                    for j in range(len(nodesEdgesList[startIdx])):
                        e2=nodesEdgesList[startIdx][j]
                        if e[1]==e2[1]:
                            isRepeated=True
                            origEdgeIndex=nodesEdgesIndices[startIdx][j]
                            edgesFaces[origEdgeIndex].append(i)
                            if facesEdgesIndices[i]== None:
                                facesEdgesIndices[i]=[origEdgeIndex]
                            else:
                                facesEdgesIndices[i].append(origEdgeIndex)
                    if not isRepeated:
                        edgeIndex=len(edgesList)
                        nodesEdgesList[startIdx].append(e)
                        nodesEdgesIndices[startIdx].append(
                            edgeIndex
                        )
                        if nodesEdgesList[endIdx] == None:
                            nodesEdgesList[endIdx]=[e]
                            nodesEdgesIndices[endIdx]=[edgeIndex]
                        else:
                            nodesEdgesList[endIdx].append(e)
                            nodesEdgesIndices[endIdx].append(edgeIndex)
                        edgesList.append(e)
                        edgesFaces.append([i])
                        if facesEdgesIndices[i]== None:
                                facesEdgesIndices[i]=[edgeIndex]
                        else:
                            facesEdgesIndices[i].append(edgeIndex)
        i:int=0                
        for edge in edgesList:
            n0,n1=edge[0],edge[1]
            if nodesNodesList[n0]== None:
                nodesNodesList[n0]=[n1]
            else:
                nodesNodesList[n0].append(n1)
            if nodesNodesList[n1]== None:
                nodesNodesList[n1]=[n0]
            else:
                nodesNodesList[n1].append(n0)
            if len(edgesFaces[i])==1:
                isBoundaryEdge.append(True)
                isBoundaryNode[n0]=True
                isBoundaryNode[n1]=True
            else:
                isBoundaryEdge.append(False)
            i+=1

        self.edgesList=edgesList
        self.nodesEgesList=nodesEdgesList
        self.nodesEdgesIndices=nodesEdgesIndices
        self.edgesFaces=edgesFaces
        self.facesEdgesIndices=facesEdgesIndices
        self.nodesNodesList=nodesNodesList
        self.isBoundaryEdge=isBoundaryEdge
        self.isBoundaryNode=isBoundaryNode
        self.elementTypes=elementTypes
        self.centersList=centersList

    def __isAdjacent(self,firstEdgeIndex:int,secondEdgeIndex:int)->bool:
        e1=self.edgesList[firstEdgeIndex]
        e2=self.edgesList[firstEdgeIndex]
        return isAdjacent(e1=e1, e2=e2)

    def surfaceMeshToVTK(self,
                    arrays=None,
                    arrays_names=None,
                    is_nodal=None
                ):

        block = vtk.vtkUnstructuredGrid()
        vtkPoints = vtk.vtkPoints()
        i:int=0
        for node in self.nodesList:
            vtkPoints.InsertPoint(i,node[0], node[1], node[2])
            i+=1
        for face in self.facesList:
            aFace = vtk.vtkPolygon()
            aFace.GetPointIds().SetNumberOfIds(len(face))
            for q in range(len(face)):
                nodeIndex:int=face[q]
                aFace.GetPointIds().SetId(q, nodeIndex)
            block.InsertNextCell(aFace.GetCellType(), aFace.GetPointIds())
        block.SetPoints(vtkPoints)

        if arrays_names!= None:
            for i in range(len(arrays_names)):
                block=vtkTools.vtkTools().npToVtk(
                    arr=arrays[i],
                    name=arrays_names[i],
                    grid=block,
                    nodal=is_nodal[i],
                    cell=(not is_nodal[i])
                    )
        return block

    def writeVTK(self,
                    name:str="surfaceMesh.vtu",
                    arrays=None,
                    arrays_names=None,
                    is_nodal=None
                )->None:

        block = self.surfaceMeshToVTK(arrays,arrays_names,is_nodal)
        writer=vtk.vtkXMLDataSetWriter()
        writer.SetFileName("build/levels/"+name)
        writer.SetInputData(block)
        writer.Write()

    def triangulate(self):
        triangles=[]
        nPts=len(self.nodesList)
        i:int=0
        for f in self.facesList:
            for edgeIndex in self.facesEdgesIndices[i]:
                triangles.append(
                    [
                        self.edgesList[edgeIndex][0],
                        self.edgesList[edgeIndex][1],
                        nPts+i
                        ]
                )
            i+=1
        trianglesNodes=list(copy.deepcopy(self.nodesList))
        for c in self.centersList:
            trianglesNodes.append(c)
        return surfaceMesh(nodesList=trianglesNodes, facesList=triangles)

    def meanEdgeLength(self):
        l=0.0
        for e in self.edgesList:
            l+=np.linalg.norm(self.nodesList[e[0]]-self.nodesList[e[1]])
        return (l/len(self.edgesList))
    
    def __isCollinear(self,twoEdgesIndices,tol=1e-5):
        e1=self.edgesList[twoEdgesIndices[0]]
        e2=self.edgesList[twoEdgesIndices[1]]
        v1=self.nodesList[e1[0]]-self.nodesList[e1[1]]
        v2=self.nodesList[e2[0]]-self.nodesList[e2[1]]
        if (np.linalg.norm(v1)*np.linalg.norm(v2))>tol:
            check=abs(np.dot(v1,v2))/(np.linalg.norm(v1)*np.linalg.norm(v2))
        else:
            check=0.0
        return True if (check > 1-tol) else False

    def edgeLength(self,edgeIndex:int)->float:
        edge=self.edgesList[edgeIndex]
        n1=edge[0]
        n2=edge[1]
        p1=np.asarray(self.nodesList[n1])
        p2=np.asarray(self.nodesList[n2])
        l=np.linalg.norm(p1-p2)
        return l

    def __edgeLocalInFace(self,faceIdx,edgeIndex):
        n1=self.edgesList[edgeIndex][0]
        n2=self.edgesList[edgeIndex][1] 
        
        iLocal1=-1
        iLocal2=-1

        i:int=0
        for n in self.facesList[faceIdx]:
            if n == n1:
                iLocal1=i
                break
            i+=1
        
        i:int=0
        for n in self.facesList[faceIdx]:
            if n == n2:
                iLocal2=i
                break
            i+=1
        
        return np.sort([iLocal1,iLocal2])

    def __averageStencil(self,nodeIdx:int):
        stencil=self.nodesNodesList[nodeIdx]
        p=np.full(shape=3, fill_value=0.0, dtype="float")
        w=1/len(stencil)
        for neighbor in stencil:
            p+=w*(np.asarray(self.nodesList[neighbor]))
        return p

    def fixOrientation(self,refCenter):
        refCenter=refCenter
        for faceIdx in range(len(self.facesList)):
            orientedFace=gm.reOrderPoints(
                ptsIndices=self.facesList[faceIdx],
                points=[self.nodesList[i] for i in self.facesList[faceIdx]],
                refCtr=refCenter
            )
            self.facesList[faceIdx]=orientedFace

    def __triangulateFace(self,faceIdx:int):

        if len(self.facesList[faceIdx])<3:
            raise ValueError("Invalid face rank")
            
        else:
            triangles=[]
            for edgeIndex in self.facesEdgesIndices[faceIdx]:
                    triangles.append(
                        [
                            self.edgesList[edgeIndex][0],
                            self.edgesList[edgeIndex][1],
                            len(self.nodesList)+faceIdx
                            ]
                    )
            return triangles
    
    def boundaryNodes(self):
        nodes=[]
        i:int=0
        for isBoundary in self.isBoundaryNode:
            if isBoundary:
                nodes.append(self.nodesList[i])
            i+=1
        return nodes

    def isTriangularMesh(self):
        for face in self.facesList:
            if len(face) != 3:
                return False
        return True

    def nodesToFacesInterpolation(self,array):
        self.barycentricTirangleCenters()
        interpolatedField=[None]*len(self.facesList)
        i=0
        for face in self.facesList:
            value=0.0
            j=0
            for node in face:
                value+=self.barycentricWeights[i][j]*array[node]
                j+=1
            interpolatedField[i]=value
            i+=1
        return interpolatedField

    def barycentricTirangleCenters(self):
        nTriangles=len(self.facesList)
        self.barycentricWeights=[None]*nTriangles
        for i in range(nTriangles):
            triangle=self.facesList[i]
            center=self.centersList[i]
            p0=self.nodesList[triangle[0]]
            p1=self.nodesList[triangle[1]]
            p2=self.nodesList[triangle[2]]
            areaTotal=triangleArea(p0, p1, p2)
            area0=triangleArea(center, p1, p2)
            area1=triangleArea(center, p0, p2)
            area2=triangleArea(center, p0, p1)
            self.barycentricWeights[i]=[area0/areaTotal,area1/areaTotal,area2/areaTotal]

def triangulateFace(faceVerticesList):
    if len(faceVerticesList)<3:
        raise ValueError("Invalid face rank")
        
    else:
        center=np.full(3,0.0,dtype="float")
        for pt in faceVerticesList:
            center+=pt
        center/=len(faceVerticesList)

        triangles=[]
        i:int=0
        for pt in faceVerticesList:
                triangles.append(
                    [
                        pt,
                        faceVerticesList[i+1//len(faceVerticesList)],
                        center
                        ]
                )
        return triangles

def isAdjacent(e1,e2)->bool:
    if e1[0] in e2:
        return True
    if e1[1] in e2:
        return True
    return False

def triangleArea(p0,p1,p2):
    ab=p1-p0
    ac=p2-p0
    area=0.5*np.linalg.norm(np.cross(ab, ac))
    return area+1e-7