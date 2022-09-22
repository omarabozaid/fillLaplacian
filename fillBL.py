from geometry import isIntersecting
from loader import *

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
            if not isValidCell(cell, cellPoints):
                cell = fixCellOrientation(cell)
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

def isValidCell(nodesIndices,cellPoints):
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
            return False
        else:
            return True
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
            return False
        else:
            return True

def fixCellOrientation(nodesIndices):
    if(len(nodesIndices)==8):
        return [nodesIndices[3],nodesIndices[2],nodesIndices[1],nodesIndices[0],nodesIndices[7],nodesIndices[6],nodesIndices[5],nodesIndices[4]]
    elif(len(nodesIndices)==6):
        return [nodesIndices[2],nodesIndices[1],nodesIndices[0],nodesIndices[5],nodesIndices[4],nodesIndices[3]]
    
def fixFlippedEdges(refMesh, fixMesh):
    for _ in range(2):
        edges = fixMesh.edgesList
        fixPoints = fixMesh.nodesList
        refPoints = refMesh.nodesList
        for edge in edges:
            u = refPoints[edge[0]]-refPoints[edge[1]]
            v = fixPoints[edge[0]]-fixPoints[edge[1]]
            if (np.dot(u,v)<0):
                fixPoints[edge[0]],fixPoints[edge[1]]=fixPoints[edge[1]],fixPoints[edge[0]]
        fixMesh.nodesList = fixPoints
    return fixMesh

def fixBadProjection(refMesh, fixMesh):
    fixPoints = fixMesh.nodesList
    refPoints = refMesh.nodesList
    faces = fixMesh.facesList
    for face in faces:
        if (len(face)==3):
            n0, n1, n2 = face
            faceEdges = [[n0,n1], [n1,n2], [n2,n0]]
            corners = [n2,n0,n1]
            for counter in range(3):
                a=copy.deepcopy(refPoints[faceEdges[counter][0]])
                b=copy.deepcopy(refPoints[faceEdges[counter][1]])
                c=copy.deepcopy(fixPoints[faceEdges[counter][0]])
                d=copy.deepcopy(fixPoints[faceEdges[counter][1]])
                refC = refPoints[corners[counter]]
                fixC = fixPoints[corners[counter]]
                e=gm.projectPointOntoLine(refC, [a,b])
                f=gm.projectPointOntoLine(fixC, [c,d])
                if (gm.isOnLine([a,b], e) and gm.isOnLine([c,d], f)):
                    v1=refC-e
                    v1/=np.linalg.norm(v1)
                    v2=fixC-f
                    v2/=np.linalg.norm(v2)
                    if np.dot(v1,v2) < 0 :
                        direction = refC-e
                        direction /= np.linalg.norm(direction)
                        mag = np.linalg.norm(fixC-f)
                        fixPoints[corners[counter]] = f+mag*direction
    fixMesh.nodesList=fixPoints
    return fixMesh

def repairMesh(refMesh, fixMesh, nIterations):
    for _ in range(nIterations):
        fixMesh = fixFlippedEdges(refMesh,fixMesh)
        fixMesh = fixBadProjection(refMesh,fixMesh)
    return fixMesh

def main():
    start()
    f = open('parameters.json')
    data = json.load(f)
    f.close()
    tools = vtkTools.vtkTools()
    nLayers = int(data['nLayers'])
    nMaxRepairIter = int(data['nMaxRepairIter'])
    cutMeshFile = data['ref_mesh']
    advectedMeshFile = data['advected_mesh']
    cutMeshPolyData = (reader.reader(cutMeshFile)).polyData()
    advectedMeshPolyData = (reader.reader(advectedMeshFile)).polyData()
    cutVertices = tools.points(cutMeshPolyData)
    cutElements = tools.cells2(cutMeshPolyData)
    advectedVertices = tools.points(advectedMeshPolyData)
    advectedElements = tools.cells2(advectedMeshPolyData)
    if(data['sym']):
        axis = -1
        if data["symAxis"]=="x":
            axis = 0
        elif data["symAxis"]=="y":
            axis = 1
        else:
            axis=2
        for i in range(len(advectedVertices)):
            if cutVertices[i][axis] == 0 :
                advectedVertices[i][axis]=0
    cutMesh = sM.surfaceMesh(cutVertices,cutElements)
    advectedMesh = sM.surfaceMesh(advectedVertices,advectedElements)
    advectedMesh = repairMesh(cutMesh, advectedMesh, nMaxRepairIter)
    advectedMesh.writeVTK("fixed.vtu")
    makeVolMeshTwoSurfs(cutMesh, advectedMesh, nLayers)


if __name__ == '__main__':
    main()
