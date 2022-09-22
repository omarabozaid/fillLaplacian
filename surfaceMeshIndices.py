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
        refNodesEdgesList=None,
        refNodesEdgesIndices=None,
        refNodesNodesList=None,
        refEdgesList=None,
        refEdgesFaces=None,
        refIsBoundaryEdge=None,
        refIsBoundaryNode=None,
        refFacesEdgesIndices=None,
        refElementTypes=None,
        refCentersList=None
    ) -> None:
        self.nodesList=np.asarray(nodesList)
        self.facesList=np.asarray(facesList)
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
        self.__inflateBB(1.1)
        if buildTopology:
            self.__buildTopology()  
        else:
            self.nodesEdgesList = refNodesEdgesList
            self.nodesEdgesIndices = refNodesEdgesIndices
            self.nodesNodesList = refNodesNodesList
            self.edgesList = refEdgesList
            self.edgesFaces = refEdgesFaces
            self.isBoundaryEdge = refIsBoundaryEdge
            self.isBoundaryNode = refIsBoundaryNode
            self.facesEdgesIndices = refFacesEdgesIndices
            self.elementTypes = refElementTypes
            self.centersList = refCentersList
    
    def topologyInformation(self):
        return (
                self.nodesEdgesList, 
                self.nodesEdgesIndices, 
                self.nodesNodesList, 
                self.edgesList ,
                self.edgesFaces ,
                self.isBoundaryEdge,
                self.isBoundaryNode,
                self.facesEdgesIndices,
                self.elementTypes,
                self.centersList
        )
        
    def __inflateBB(self,lamda:float=1.1):
        center=0.5*(self.leftCorner+self.rightCorner)
        self.leftCorner=center-lamda*(center-self.leftCorner)
        self.rightCorner=center+lamda*(self.rightCorner-center)

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

    def smoothMesh(self,
    surfaceMeshSmoothingDict:dict,
    isFree
    )->None:
        iter:int=0
        smoother=surfaceMeshSmoothingDict.get('smoother')
        maxSmoothingIter=surfaceMeshSmoothingDict.get('maxSmoothingIter')
        smoothBoundaries=surfaceMeshSmoothingDict.get('smoothBoundaries')

        if smoother == "power":
            rf=surfaceMeshSmoothingDict.get('rf')
            exponent=surfaceMeshSmoothingDict.get('exponent')
            while(iter<maxSmoothingIter):
                #print("SMOOTHING ITER ",iter)
                self.__powerSmoothMesh(
                    rf=rf,exponent=exponent,smoothBoundaries=smoothBoundaries,
                    isFree=isFree
                )
                iter+=1
        elif smoother == "weighted":
            rf=surfaceMeshSmoothingDict.get('rf')
            exponent=surfaceMeshSmoothingDict.get('exponent')
            while(iter<maxSmoothingIter):
                #print("SMOOTHING ITER ",iter)
                self.__weightedLaplacian(
                        rf=rf,
                        isFree=isFree
                )
                iter+=1
        
        elif smoother == "laplacian":
            rf=surfaceMeshSmoothingDict.get('rf')
            while(iter<maxSmoothingIter):
                #print("SMOOTHING ITER ",iter)
                self.__laplacian(rf=rf,isFree=isFree)
                iter+=1
        
        elif smoother == "HC":
            alpha=surfaceMeshSmoothingDict.get('alpha')
            beta=surfaceMeshSmoothingDict.get('beta')
            origPts=copy.deepcopy(self.nodesList)
            while(iter<maxSmoothingIter):
                #print("SMOOTHING ITER ",iter)
                self.__HCSmoothMesh(o=origPts,alpha=alpha,beta=beta,isFree=isFree)
                gc.collect()
                iter+=1
            gc.collect()

    def __laplacian(self,rf:float=0.25,isFree=[None]):
        nPts=len(self.nodesList)
        p=[None]*nPts
        i:int=0
        for stencil in self.nodesNodesList:
            p[i]=np.full(shape=3, fill_value=0.0, dtype="float")
            w=1/len(stencil)
            for neighbor in stencil:
                p[i]+=w*(np.asarray(self.nodesList[neighbor]))
            i+=1
        for i in range(nPts):
            self.nodesList[i]=np.asarray(rf*p[i]+(1-rf)*np.asarray(self.nodesList[i]))
            if isFree[i]==False:
                self.nodesList[i][1]=0.0

    def __powerSmoothMesh(self,rf:float=0.25,exponent:"float"=1,smoothBoundaries:bool=True,isFree=[None])->None:
        newNodesList=[None]*len(self.nodesList)
        i:int=0
        for stencil in self.nodesNodesList:
            den=0
            newPos=np.full(3,0.0,dtype="float")
            vertex=self.nodesList[i]
            for neighbor in stencil:
                d=np.power(np.linalg.norm(vertex-self.nodesList[neighbor]),exponent)
                den+=d
                newPos+=d*self.nodesList[neighbor]
            if den>1e-5:
                newPos/=den
            else:
                newPos=vertex
            if smoothBoundaries:
                if isFree[i]:
                    newNodesList[i]=newPos*rf+vertex*(1-rf)
            else:
                if self.isBoundaryNode[i]:
                    newNodesList[i]=vertex
                elif isFree[i]:
                    newNodesList[i]=newPos*rf+vertex*(1-rf)
            i+=1
        self.nodesList=newNodesList

    def __HCSmoothMesh(self,o,alpha,beta,isFree=[None]):
        nPts=len(self.nodesList)

        p=[None]*nPts

        alpha=alpha
        beta =beta

        b=[None]*nPts
        d=[None]*nPts
        
        i:int=0
        for stencil in self.nodesNodesList:
            p[i]=0.0*self.nodesList[i]
            w=1/len(stencil)
            for neighbor in stencil:
                p[i]+=w*(self.nodesList[neighbor])
            b[i]=p[i]-(alpha*o[i]+(1-alpha)*self.nodesList[i])
            i+=1

        i:int=0
        for stencil in self.nodesNodesList:
            d[i]=beta*b[i]
            for neighbor in stencil:
                d[i]+=(((1-beta)/len(stencil))*b[neighbor])
            d[i]*=-1
            if isFree[i]:
                self.nodesList[i]=p[i]+d[i]
            i+=1
    
    def repairSelfIntersectingFaces(self,maxIter:int=50)->None:
        iter:int=0
        while iter<=maxIter:
            selfIntersectingFacesList=self.__selfIntersectingFaces()
            if len(selfIntersectingFacesList)==0:
                break
            else:
                for faceIndex in selfIntersectingFacesList:
                    self.__fixSelfIntersectingFace(faceIndex=faceIndex)
                self.__buildTopology()
            iter+=1

    def naiveRepairSelfIntersectingFaces(self,
        surfaceMeshControlsDict:dict,
        ctrs
    )->None:
        iter:int=0
        while iter<surfaceMeshControlsDict.get("maxIterFaceCorrection"):
            print("SELF-INTERSECTION ITERATION ",iter)
            for faceIndex in range(len(self.facesList)):
                self.__fixSelfIntersectingFace(faceIndex=faceIndex,refCtr=ctrs[faceIndex])
            self.__buildTopology()
            iter+=1

    def __fixSelfIntersectingFace(self,faceIndex,refCtr)->None:
        wrongFace=self.facesList[faceIndex]
        pts=[]
        for i in range(len(wrongFace)):
            pts.append(self.nodesList[wrongFace[i]])
        correctedFace=gm.reOrderPoints(
            ptsIndices=wrongFace,
            points=pts,
            refCtr=refCtr
        )
        self.facesList[faceIndex]=correctedFace

    def __selfIntersectingFaces(self)->None:
        facesList=self.facesList
        nFaces=len(facesList)
        selfIntersectingFaces=[]
        for i in range(nFaces):
            if self.__isSelfIntersecting(faceIndex=i):
                selfIntersectingFaces.append(i)
        return selfIntersectingFaces
    
    def __isSelfIntersecting(self,faceIndex)->bool:
        edgesIndices=self.facesEdgesIndices[faceIndex]
        edges=[]
        for edgeIndex in edgesIndices:
            edges.append(self.edgesList[edgeIndex])
        for i in range(len(edges)):
            node1Index=edges[i][0]
            node2Index=edges[i][1]
            line1=[
                self.nodesList[node1Index],self.nodesList[node2Index]
                ]
            edge1=edges[i]
            for j in range(len(edges)):
                edge2=edges[j]
                if not self.__isAdjacent(e1=edge1,e2=edge2):
                    node3Index=edges[j][0]
                    node4Index=edges[j][1]
                    line2=[
                        self.nodesList[node3Index],self.nodesList[node4Index]
                        ]
                    if(
                        gm.isIntersecting(
                        line1=line1,
                        line2=line2
                            )
                        ):
                        return True
        return False

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

    def serialize(self,dx):
        self.dx=dx
        self.n=ptIJK(pt=self.rightCorner,left=self.leftCorner,dx=self.dx)
        nObjects=(self.n[0]+1)*(self.n[1]+1)*(self.n[2]+1)
        serialized=[None]*nObjects
        j:int=0
        for pt in self.nodesList:
            ijk=ptIJK(pt=pt,left=self.leftCorner,dx=self.dx)
            ptIndex=ijkToIndex(ijk=ijk,n=self.n)
            if serialized[ptIndex]==None:
                serialized[ptIndex]=[j]
            else:
                serialized[ptIndex].append(j)
            j+=1
        self.serialized=serialized
        nObjects=(self.n[0]+1)*(self.n[1]+1)*(self.n[2]+1)
        serialized=[None]*nObjects
        j:int=0
        for pt in self.centersList:
            ijk=ptIJK(pt=pt,left=self.leftCorner,dx=self.dx)
            ptIndex=ijkToIndex(ijk=ijk,n=self.n)
            if serialized[ptIndex]==None:
                serialized[ptIndex]=[j]
            else:
                serialized[ptIndex].append(j)
            j+=1
        self.serializedTriangles=serialized

    def meanEdgeLength(self):
        l=0.0
        for e in self.edgesList:
            l+=np.linalg.norm(self.nodesList[e[0]]-self.nodesList[e[1]])
        return (l/len(self.edgesList))

    def isEnclosed(self,testPoint):
        
        checkRights:bool=[False]*3
        checkLefts:bool=[False]*3

        testIJK=ptIJK(
                pt=testPoint,
                left=self.leftCorner,
                dx=self.dx
            )
        isInBox:bool=True
        k:int=0
        for idx in testIJK:
            if idx>=self.n[k] or idx<0:
                isInBox=False
                break
            k+=1

        if isInBox:
            for axis in range(3):
                testIJK=ptIJK(
                            pt=testPoint,
                            left=self.leftCorner,
                            dx=self.dx
                        )
                testIndex=ijkToIndex(ijk=testIJK,n=self.n)

                if self.serialized[testIndex] != None:
                    indices=self.serialized[testIndex]
                    pt=[]
                    for idx in indices:
                        pt.append(self.nodesList[idx])
                    nearestPt=gm.closestPoint(testPoint,pt)
                    if nearestPt[axis]>=testPoint[axis]:
                        checkRights[axis]=True
                if checkRights[axis] == False:
                    i=testIJK[axis]+1
                    while i<self.n[axis]:
                        testIJK[axis]=i
                        testIndex=ijkToIndex(ijk=testIJK,n=self.n)
                        if self.serialized[testIndex] != None:
                            checkRights[axis]=True
                            break
                        i+=1
            
                testIJK=ptIJK(
                            pt=testPoint,
                            left=self.leftCorner,
                            dx=self.dx
                        )
                testIndex=ijkToIndex(ijk=testIJK,n=self.n)

                if self.serialized[testIndex] != None:
                    indices=self.serialized[testIndex]
                    pt=[]
                    for idx in indices:
                        pt.append(self.nodesList[idx])
                    nearestPt=gm.closestPoint(testPoint,pt)
                    if nearestPt[axis]<=testPoint[axis]:
                        checkLefts[axis]=True
                if checkLefts[axis] == False: 
                    i=testIJK[axis]-1
                    while i>=0:
                        testIJK[axis]=i
                        testIndex=ijkToIndex(ijk=testIJK,n=self.n)
                        if self.serialized[testIndex] != None:
                            checkLefts[axis]=True
                            break
                        i-=1
        Probs=0
        for check in checkRights:
            if check == True:
                Probs+=1
        for check in checkLefts:
            if check == True:
                Probs+=1
        return True if Probs>3 else False

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
    return area