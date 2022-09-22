from geometry import isIntersecting
from loader import *

from shapely.geometry import Point,Polygon

def makeVolMeshTwoSurfs(upperMesh, lowerMesh, nLayers):
    nLines = len(upperMesh.nodesList)
    pathes= []
    
    for i in range(nLines):
        path = []
        start_pt = upperMesh.nodesList[i]
        end_pt = lowerMesh.nodesList[i]
        line = end_pt-start_pt
        dL = line/nLayers
        for j in range(nLayers+1):
            cut_pt = start_pt+j*dL
            path.append(cut_pt)
        pathes.append(path) 
    
    levels = [[None]*nLines for _ in range(nLayers+1)]
    for j in range(nLines):
        for i in range(nLayers+1):
            levels[i][j]=pathes[j][i]
    
    points=[]
    for level in levels:
        for node in level:
            points.append(np.asarray(node))

    cells = []
    for face in upperMesh.facesList:
        for layer in range(nLayers):
            upper_ij = layer*nLines
            lower_ij = (layer+1)*nLines
            upper_face = []
            lower_face = []
            for node in face:
                upper_face.append(node+upper_ij)
                lower_face.append(node+lower_ij)
            cell = upper_face
            for node in lower_face:
                cell.append(node)
            cellPoints = [points[idx] for idx in cell]
            cell=fixCellOrientation(cell, cellPoints)
            cells.append(cell)
    
    boundaryFaces=[f for f in upperMesh.facesList]
    boundaryFacesRegions = [1 for i in range(len(upperMesh.facesList))]
    for face in upperMesh.facesList:
        lower_ij = nLayers*nLines
        lower_face = []
        for node in face:
            lower_face.append(node+lower_ij)
        boundaryFaces.append(lower_face)
        boundaryFacesRegions.append(2)
    
    nNodes = int(len(pathes[0])*len(pathes))
    nCells = int(len(cells))
    nBoundaryFaces = int(len(boundaryFaces))
    nElements = nBoundaryFaces+nCells

    with open('build/levels/boundary_volume.msh', 'w') as data:
        fw = data.write
        fw("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fw("$PhysicalNames\n3\n2 1 \"cut_faces_laplacian_side\"\n2 2 \"boat\"\n3 3 \"fluid\"\n$EndPhysicalNames\n")
        fw("$Nodes\n%i\n"%nNodes)
        counter = 1
        for level in levels:
            for node in level:
                fw("%i %f %f %f \n"%(counter,node[0],node[1],node[2]))
                counter+=1
        fw("$EndNodes\n")
        fw("$Elements\n%i\n"%(nElements))
        counter = 1
        for face in boundaryFaces:
            elemRegion = boundaryFacesRegions[counter-1]
            elemType = -1
            if(len(face)==3):
                elemType=2
                fw("%i %i 1 %i %i %i %i \n"%(counter,elemType,elemRegion,face[0]+1,face[1]+1,face[2]+1))
            elif(len(face)==4):
                elemType=3
                fw("%i %i 1 %i %i %i %i %i \n"%(counter,elemType,elemRegion,face[0]+1,face[1]+1,face[2]+1,face[3]+1))
            else:
                print("wrong face")
            counter+=1
        for cell in cells:
            elemRegion = 3
            elemType = -1
            if(len(cell)==6):
                elemType=6
                fw("%i %i 1 %i %i %i %i %i %i %i \n"%(counter,elemType,elemRegion,cell[0]+1,cell[1]+1,cell[2]+1,cell[3]+1,cell[4]+1,cell[5]+1))
            elif(len(cell)==8):
                elemType=5
                fw("%i %i 1 %i %i %i %i %i %i %i %i %i \n"%(counter,elemType,elemRegion,cell[0]+1,cell[1]+1,cell[2]+1,cell[3]+1,cell[4]+1,cell[5]+1,cell[6]+1,cell[7]+1))
            else:
                print("wrong cell")
            counter+=1
        fw("$EndElements\n")

def fixCellOrientation(nodesIndices,cellPoints):
    if(len(nodesIndices)==8):
        u=cellPoints[1]-cellPoints[0]
        u/=np.linalg.norm(u)
        v=cellPoints[2]-cellPoints[1]
        v/=np.linalg.norm(v)
        xC=np.full(3,0,dtype='float')
        xF=np.full(3,0,dtype='float')
        
        for pt in cellPoints:
            xC+=pt

        for counter in range(4):
            xF+=cellPoints[counter]

        xC/=8
        xF/=4
        w=xC-xF
        w/=np.linalg.norm(w)
        n=np.cross(u, v)
        check = np.dot(n, w)
        if(check<0):
            return [nodesIndices[3],nodesIndices[2],nodesIndices[1],nodesIndices[0],nodesIndices[7],nodesIndices[6],nodesIndices[5],nodesIndices[4]]
        else:
            return nodesIndices
    elif(len(nodesIndices)==6):
        u=cellPoints[1]-cellPoints[0]
        u/=np.linalg.norm(u)
        v=cellPoints[2]-cellPoints[0]
        v/=np.linalg.norm(v)
        w=cellPoints[3]-cellPoints[0]
        w/=np.linalg.norm(w)
        n=np.cross(u, v)
        check = np.dot(n, w)
        if(check<0):
            return [nodesIndices[2],nodesIndices[1],nodesIndices[0],nodesIndices[5],nodesIndices[4],nodesIndices[3]]
        else:
            return nodesIndices

def fixSelfIntersections(refMesh, fixMesh):
    edges = refMesh.edgesList
    refPoints = refMesh.nodesList
    fixPoints = fixMesh.nodesList
    nFixedEdges = 0
    for edge in edges:
        n0=edge[0]
        n1=edge[1]
        a=copy.deepcopy(refPoints[n0])
        b=copy.deepcopy(refPoints[n1])
        c=copy.deepcopy(fixPoints[n0])
        d=copy.deepcopy(fixPoints[n1])
        line1=[a,c]
        line2=[b,d]
        if gm.isIntersecting(line1,line2):
            print("self intersection")
            nFixedEdges+=1
            fixPoints[n0]=d
            fixPoints[n1]=c


    return fixMesh, nFixedEdges


def main():
    start()
    tools = vtkTools.vtkTools()
    nLayers = 5
    nMaxSwappingIter = 1

    cutMeshFile = "data/dtmb/ref_mesh.vtk"
    advectedMeshFile = "data/dtmb/advected_mesh.vtk"
    cutMeshPolyData = (reader.reader(cutMeshFile)).polyData()
    advectedMeshPolyData = (reader.reader(advectedMeshFile)).polyData()
    cutVertices = tools.points(cutMeshPolyData)
    cutElements = tools.cells2(cutMeshPolyData)
    advectedVertices = tools.points(advectedMeshPolyData)
    advectedElements = tools.cells2(advectedMeshPolyData)
    for i in range(len(advectedVertices)):
        if cutVertices[i][1] == 0 :
            advectedVertices[i][1]=0
    cutMesh = sM.surfaceMesh(cutVertices,cutElements)
    advectedMesh = sM.surfaceMesh(advectedVertices,advectedElements)
    
    iter = 0
    while(iter<nMaxSwappingIter):
        advectedMesh, nFixedEdges = fixSelfIntersections(cutMesh, advectedMesh)
        if nFixedEdges>0:
            print("number of swapped edges is %i in the %i-th iteration"%(nFixedEdges,iter))
        if nFixedEdges == 0:
            break
        iter+=1
    
    if iter> 0:
        advectedMesh.writeVTK("swapped.vtu")
    makeVolMeshTwoSurfs(cutMesh, advectedMesh, nLayers)


if __name__ == '__main__':
    main()