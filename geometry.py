import numpy as np
from numpy.core.umath_tests import inner1d

"""
Hard wired stability const
"""
epsilon=1e-10


def projectPointOntoLine(
        point:np.ndarray,
        line:np.ndarray
    )->np.ndarray:
    
    A:np.ndarray=line[0]
    B:np.ndarray=line[1]
    AB:np.ndarray=B-A
    P:np.ndarray=point
    d:np.ndarray = (AB) / np.linalg.norm(AB)
    v:np.ndarray = P - B
    t:float = np.dot(v,d)
    return B + t * d

def projectionSegmentLine(
        segment:np.ndarray,
        line:np.ndarray
    )-> float:

    p1=projectPointOntoLine(point=segment[0],line=line)
    p2=projectPointOntoLine(point=segment[1],line=line)
    return np.linalg.norm(p1-p2)

def nearestLine(lines:np.ndarray,point:np.ndarray)->int:
    dist=np.full(len(lines),0.0,dtype="float")
    i:int=0
    for i in range(len(lines)):
        projectedPoint=projectPointOntoLine(line=lines[i],point=point)
        if isOnLine(line=lines[i],point=projectedPoint) is False:
            dist[i]=1e20
        else:
            dist[i]=np.linalg.norm(projectedPoint-point)
    return np.argmin(dist)

def isOnLine(line:np.ndarray,point:np.ndarray):
    x=[line[0][0],line[1][0]]
    y=[line[0][1],line[1][1]]
    z=[line[0][2],line[1][2]]
    if point[0]>= np.min(x) and point[0]<= np.max(x):
        if point[1]>= np.min(y) and point[1]<= np.max(y):
            if point[2]>= np.min(z) and point[2]<= np.max(z):
                return True
    else:
        return False

def isInternalRegion(
        cellCenter,
        triangleCenter,
        triangleNormal
        )->bool:
        sign=np.dot(
            triangleNormal, 
            cellCenter-triangleCenter
        )
        return True if sign <=0 else False

def closestPoint(point,points):
    dist=[]
    for pt in points:
        dist.append(np.linalg.norm(np.asarray(point)-np.asarray(pt)))
    idxMin=np.argmin(dist)
    return points[idxMin]

def closestPointIdx(point,points):
    dist=[]
    for pt in points:
        dist.append(np.linalg.norm(np.asarray(point)-np.asarray(pt)))
    return np.argmin(dist)

def projectPointTriangles(point,triangles,points):
    projectedPoints=[]
    for triangle in triangles:
        pt=projectPointTriangle(
            point=point,
            triangle=triangle,
            points=points
        )
        projectedPoints.append(pt)
    return projectedPoints

def projectPointTriangle(point,triangle,points):
        
    
    v0=points[triangle[0]]
    v1=points[triangle[1]]
    v2=points[triangle[2]]

    diff = v0-point
    edge0=v1-v0
    edge1=v2-v0

    a00 = np.dot(edge0, edge0)
    a01 = np.dot(edge0, edge1)
    a11 = np.dot(edge1, edge1)
    b0 = np.dot(diff, edge0)
    b1 = np.dot(diff, edge1)
    det =max(a00 * a11 - a01 * a01, 0.0)
    s = a01 * b1 - a11 * b0
    t = a01 * b0 - a00 * b1

    if (s + t <= det):
        if (s < 0.0):
            if (t < 0.0):
                if (b0 < 0.0):
                    t = 0.0
                    if (-b0 >= a00) :
                        s = 1.0
                    else:
                        s = -b0 / a00
                else:
                    s = 0.0
                    if (b1 >= 0.0):
                        t = 0.0
                    elif (-b1 >= a11):
                        t = 1.0
                    else:
                        t = -b1 / a11
                    
            else:
                s = 0.0
                if (b1 >= 0.0) :
                    t = 0.0
                elif (-b1 >= a11) :
                    t = 1.0
                else:
                    t = -b1 / a11

        elif (t < 0.0) :
            t = 0.0
            if (b0 >= 0.0):
                s = 0.0
            elif (-b0 >= a00):
                s = 1.0
            else:
                s = -b0 / a00
        else:
            s /= (det+1e-9)
            t /= (det+1e-9)
    else:
        if (s < 0.0):
            tmp0 = a01 + b0
            tmp1 = a11 + b1
            if (tmp1 > tmp0):
                numer = tmp1 - tmp0
                denom = a00 - 2.0 * a01 + a11
                if (numer >= denom):
                    s = 1.0
                    t = 0.0
                else:
                    s = numer / denom
                    t = 1.0 - s
            else:
                s = 0.0
                if (tmp1 <= 0.0):
                    t = 1.0
                elif (b1 >= 0.0):
                    t = 0.0
                else:
                    t = -b1 / a11

        elif (t < 0.0):
            tmp0 = a01 + b1
            tmp1 = a00 + b0
            if (tmp1 > tmp0):
                numer = tmp1 - tmp0
                denom = a00 - 2.0 * a01 + a11
                if (numer >= denom):
                    t = 1.0
                    s = 0.0
                else:
                    t = numer / denom
                    s = 1.0 - t
            else:
                t = 0.0
                if (tmp1 <= 0.0):
                    s = 1.0
                elif (b0 >= 0.0):
                    s = 0.0
                else:
                    s = -b0 / a00
        else:
            numer = a11 + b1 - a01 - b0
            if (numer <= 0.0):
                s = 0.0
                t = 1.0
            else:
                denom = a00 - 2.0 * a01 + a11
                if (numer >= denom):
                    s = 1.0
                    t = 0.0
                else:
                    s = numer / denom
                    t = 1.0 - s
    return (v0 + s * edge0 + t * edge1)

def projectPointPlane(point,triangle,points):
    vertices=points[triangle]
    planeCenter=(1/3)*(vertices[0]+vertices[1]+vertices[2])
    planeNormal=np.cross(
        vertices[0]-vertices[1],vertices[0]-vertices[2]
    )
    planeNormal/=(np.linalg.norm(planeNormal))
    return point - np.dot(point - planeCenter, planeNormal) * planeNormal

def isPointInTriangle(point,triangle,points,tol=1e-9):
    vertices=points[triangle]
    area=0.5*np.linalg.norm(np.cross(
        vertices[0]-vertices[1],vertices[0]-vertices[2]
            )
        )
    alpha=np.linalg.norm(np.cross(
        point-vertices[1],point-vertices[2]
            )
        )/(2*area)
    if alpha<=(1+tol) and alpha>=(-tol):
        beta=np.linalg.norm(np.cross(
        point-vertices[2],point-vertices[0]
            )
        )/(2*area)
        if beta<=(1+tol) and beta>=-tol:
            gamma=1-alpha-beta
            if gamma<=(1+tol) and gamma>=-tol:
                return True
    return False
    
def isIntersecting(line1,line2):
    p1,p2=line1[0],line1[1]
    p3,p4=line2[0],line2[1]
    mua = 0
    mub = 0
    p13 = p1 - p3
    p43 = p4 - p3
    p21 = p2 - p1
    d1343 = np.dot(p13,p43)
    d4321 = np.dot(p21,p43)
    d1321 = np.dot(p13,p21)
    d4343 = np.dot(p43,p43)
    d2121 = np.dot(p21,p21)
    denom = d2121 * d4343 - d4321 * d4321
    if (abs(denom) < epsilon):
      return False
    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343
    pa = p1 + p21 * mua
    pb = p3 + p43 * mub
    if isOnLine(line1,pa) and isOnLine(line2,pb):
        if np.linalg.norm(pa-pb) < epsilon:
            return True
    return False

def reOrderPoints(
    ptsIndices,
    points,
    refCtr
):
    center,normal=averagedNormal(points=points)
    refDir=(refCtr-center)/np.linalg.norm(refCtr-center)
    if np.dot(refDir,normal)<0.0:
        normal*=-1.0

    thetas=projectPointsOnPlane(
        points=points,
        planeCenter=center,
        planeNormal=normal
    )
    for i in range(len(thetas)):
        if thetas[i]<0.0:
            thetas[i]+=2*np.pi
    sortedArr=np.argsort(thetas)
    sortedIndices=[]
    for i in sortedArr:
        sortedIndices.append(ptsIndices[i])
    return sortedIndices

def reOrderPoints2(
    ptsIndices,
    points,
    normal,
    center
):

    thetas=projectPointsOnPlane(
        points=points,
        planeCenter=center,
        planeNormal=normal
    )
    for i in range(len(thetas)):
        if thetas[i]<0.0:
            thetas[i]+=2*np.pi
    sortedArr=np.argsort(thetas)
    sortedIndices=[]
    for i in sortedArr:
        sortedIndices.append(ptsIndices[i])
    return sortedIndices,thetas

def projectPointsOnPlane(
    points,
    planeCenter,
    planeNormal
):
    pts=[None]*len(points)
    radii=[None]*len(points)
    thetas=[None]*len(points)

    for i in range(len(points)):
        v=points[i]-planeCenter
        dist=np.dot(v,planeNormal)
        pts[i]=points[i]-dist*planeNormal
        radii[i]=np.linalg.norm(
            pts[i]-planeCenter
        )
        du=pts[i][0]-planeCenter[0]
        dv=pts[i][2]-planeCenter[2]
        thetas[i]=np.arctan2(dv,du)
    return thetas

def projectPointOnPlane(
    point,
    planeCenter,
    planeNormal
):
    return (point-(np.dot(point-planeCenter,planeNormal))*planeNormal)

def centerOfMass(points):
    center=np.full(3,0,dtype="float")
    for pt in points:
        center+=pt
    center/=len(points)
    return center

def averagedNormal(points):
    """
    https://en.wikipedia.org/wiki/Singular_value_decomposition
    https://math.stackexchange.com/questions/3869/what-is-the-intuitive-relationship-between-svd-and-pca/3871#3871
    """
    ctr=centerOfMass(points=points)
    X=np.full((3,len(points)),0.0,dtype="float")
    i:int=0
    for pt in points:
        X[0][i]=pt[0]-ctr[0]
        X[1][i]=pt[1]-ctr[1]
        X[2][i]=pt[2]-ctr[2]
        i+=1

    normal=np.full(3,0,dtype="float")
    M = np.dot(X, X.T)
    normal=np.linalg.svd(M)[0][:,-1]
    normal/=np.linalg.norm(normal)
    refDir=ctr/np.linalg.norm(ctr)
    if np.dot(refDir,normal)>0:
        normal*=-1
    return ctr,normal

def isRHS(axis,ptsCycle)->bool:
    p0=ptsCycle[0]
    p1=ptsCycle[1]
    pMid=0.5*(p0+p1)
    tangential=p1-p0
    
    ctr=centerOfMass(ptsCycle)
    inPlaneNormal=(pMid-ctr)/np.linalg.norm(pMid-ctr)

    normal=np.cross(inPlaneNormal, tangential)
    normal/=(np.linalg.norm(normal))

    axis/=np.linalg.norm(axis)
    
    check=np.dot(axis,normal)
    if check>=-1e-5:
        return True
    return False
