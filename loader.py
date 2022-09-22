import sys
sys.setrecursionlimit(100000000)

import numpy as np 
import pandas as pd
import time 
import scipy
import copy
import vtk as vtk
import json

import vtkTools as vtkTools
import surfaceMeshIndices as sM
import reader as reader
import geometry as gm


def start()->None:
    print("----------------------------------")
    print("-------FILL BOUNDARY LAYER--------")
    print("----------------------------------")