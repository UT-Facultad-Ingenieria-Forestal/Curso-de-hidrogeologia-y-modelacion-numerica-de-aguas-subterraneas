import math
import numpy as np
from scipy.interpolate import griddata

class vtkProcess:
    def __init__(self, name):
        self.name = name
        
    def arrayWaterTableObject(modDis,hdsObject):
        wtObject=np.zeros([modDis['cellRows'],modDis['cellCols']])
        for row in range(modDis['cellRows']):
            for col in range(modDis['cellCols']):
                heads=hdsObject[:,row,col]
                if heads[heads>-1.0000e+20].size > 0:
                    wtObject[row,col]=round(heads[heads>-1.0000e+20][0],2)
                else:
                    wtObject[row,col]=-1.0000e+20
        return wtObject   
    
    #function that return a dictionary of z values on the vertex
    def interpolateCelltoVertex(modDis,item):
        dictZVertex = {}
       
        #arrange to hace positive heads in all vertex of an active cell
        for lay in range(len(modDis[item].keys())):
            
            #print(lay)
            matrix = np.zeros([modDis['vertexRows'],modDis['vertexCols']])
            
            modDisItemMatrix = np.array(modDis[item]['lay'+str(lay)]).reshape([modDis['cellRows'],modDis['cellCols']])
            #print(matrix.shape)
            #print(modDisItemMatrix.shape)
            matrix[0,:-1] = modDisItemMatrix[0,:]
            matrix[:-1,0] = modDisItemMatrix[:,0]
            
            for row in range(1,modDis['cellRows']):
                for col in range(1,modDis['cellCols']):
                    headLay = modDisItemMatrix
                    neighcartesianlist = [headLay[row-1,col-1],headLay[row-1,col],headLay[row,col-1],headLay[row,col]]

                    headMean = sum(neighcartesianlist)/len(neighcartesianlist)
            
                    matrix[row,col]=headMean
            
            
            matrix[-1,:-1] = modDisItemMatrix[-1,:]
            matrix[:-1,-1] = modDisItemMatrix[:,-1]
            matrix[-1,-1] = modDisItemMatrix[-1,-1]

            dictZVertex['lay'+str(lay)]=matrix
            
        return dictZVertex
    
    def simplifiedVertexHead(modDis,headObject):
        matrix = np.zeros([modDis['vertexLays'],modDis['vertexRows'],modDis['vertexCols']])
        
        matrix[1:,0,:-1]  = headObject[:,0,:]  #first row
        matrix[1:,-1,:-1] = headObject[:,-1,:] #last row
        matrix[1:,:-1,0]  = headObject[:,:,0]  #first col
        matrix[1:,:-1,-1] = headObject[:,:,-1] #last col
        matrix[1:,-1,-1]  = headObject[:,-1,-1] #lower right value
        
        for lay in range(1,modDis['vertexLays']):
            for row in range(1,modDis['vertexRows']-1):
                for col in range(1,modDis['vertexCols']-1):
                    neighcartesianlist = [headObject[lay-1,row-1,col-1],
                                          headObject[lay-1,row-1,col],
                                          headObject[lay-1,row,col-1],
                                          headObject[lay-1,row,col]]
                    headList = []
                    for item in neighcartesianlist: #accumulator that only add positive values
                        if item > -1e+20:
                            headList.append(item)
                    if len(headList) > 0:
                        headMean = sum(headList)/len(headList)
                    else:
                        headMean = -1e+20
                        
                    matrix[lay,row,col]=headMean
        
        matrix[0,:,:] = matrix[1,:,:]
        
        return list(matrix.flatten())
    

    def listWaterTableCellFunction(modDis,wtObject):
        wtHeads = list(wtObject[wtObject>-1.0000e+20])
        return wtHeads
#    def listWaterTableCellFunction(modDis,headObject):
#        wtHeads = []
#        for row in range(modDis['cellRows']):
#            for col in range(modDis['cellCols']):
#                heads=headObject[:,row,col]
#                if heads[heads>-1.0000e+20].size > 0:
#                    wtHeads.append(round(heads[heads>-1.0000e+20][0],2))
                
        

    
    def listWaterTableVertexFunction(modDis,wtObject):
        wtVertex = []
        
        for row in range(modDis['vertexRows']):
                for col in range(modDis['vertexCols']):
                    #neighboringCells
                    nL=np.array([[row,col],[row-1,col],[row-1,col-1],[row,col-1]])
                    #define applicable cells - non negative, no greater than matrix dimensions
                    appCells = nL[(nL[:,0]>-1) & (nL[:,0]<modDis['cellRows']) & 
                                 (nL[:,1]>-1) & (nL[:,1]<modDis['cellCols'])]
                    
                    #empty numpy array to store data
                    cellValues = np.array([])
                    
                    #store all data for neigh cells
                    #print(appCells)
                    for cell in appCells:
                        
                        if wtObject[cell[0],cell[1]]>-1.0000e+20:
                            cellValues = np.append(cellValues,wtObject[cell[0],cell[1]])
                            
                    #check is neigh cells are empty
                    if cellValues.shape[0]>0:
                        wtVertex.append(cellValues.mean())
                    else: 
                        wtVertex.append(-1.0000e+20)
                    
        return wtVertex
    
    
##########################################    
# Lists for the VTK file of Model Geometry  
##########################################  
    
    ####Unvariability Check####
    def vertexXYZPointsFunction(modDis):
        #empty list to store all vertex XYZ
        vertexXYZPoints = []

        #definition of xyz points for all vertex
        for lay in range(modDis['vertexLays']):
            for row in range(modDis['vertexRows']):
                for col in range(modDis['vertexCols']):
                    xyz=[
                        modDis['vertexEastingArray1D'][col], 
                        modDis['vertexNorthingArray1D'][row],
                        modDis['vertexZGrid']['lay'+str(lay)][row, col]
                        ]
                    vertexXYZPoints.append(xyz)
                    
        return vertexXYZPoints
    
    
    ####Unvariability Check####
    def vertexWaterTableXYZPointsFunction(listWaterTableVertex,modDis):
        #empty list to store all vertex Water Table XYZ
        vertexWaterTableXYZPoints = []
        
        gridWaterTableVertex = np.array(listWaterTableVertex).reshape(modDis["vertexRows"],modDis["vertexCols"])

        #definition of xyz points for all vertex
        for row in range(modDis['vertexRows']):
            for col in range(modDis['vertexCols']):
                waterTable = gridWaterTableVertex[row, col]
                
                xyz=[modDis['vertexEastingArray1D'][col], 
                     modDis['vertexNorthingArray1D'][row],
                     waterTable]
                
                vertexWaterTableXYZPoints.append(xyz)    
                
        return vertexWaterTableXYZPoints
    
    
    ####Unvariability Check####
    def bcCellsListFunction(modFile,keyName,listHexaSequence,modDis,modBas):
        #definition of cells
        anyGrid = np.concatenate(modBas['active']).flatten().reshape([modDis["cellLays"],modDis["cellRows"],modDis["cellCols"]])
        
        for item in modFile[keyName]:
            if anyGrid[int(item[0]),int(item[1]),int(item[2])] > 0:
                anyGrid[int(item[0]),int(item[1]),int(item[2])] = 2
              
        listBcCellsIO = list(anyGrid[anyGrid>0])
        
        listBcCellsIODef = []
        listBcCellsSecuenceDef = []
        
        for item in range(len(listBcCellsIO)):
            if listBcCellsIO[item] == 2:
                listBcCellsIODef.append(listBcCellsIO[item])
                listBcCellsSecuenceDef.append(listHexaSequence[item])
        
        
        return [listBcCellsIODef,listBcCellsSecuenceDef]
    
##################################################   
# Hexahedrons and Quads sequences for the VTK File
##################################################  

    def listLayerQuadSequenceFunction(modDis,modBas,wtObject):
        anyQuadList = []
        
        #definition of hexahedrons cell coordinates
        for row in range(modDis['cellRows']):
            for col in range(modDis['cellCols']):
                if modBas['active'][0][row,col]==1 and wtObject[row,col]>-1.0000e+20:
                    pt0 = modDis['vertexCols']*(row+1)+col
                    pt1 = modDis['vertexCols']*(row+1)+col+1
                    pt2 = modDis['vertexCols']*(row)+col+1
                    pt3 = modDis['vertexCols']*(row)+col
                    anyList = [pt0,pt1,pt2,pt3]
                    anyQuadList.append(anyList)
                
        return anyQuadList
    
    def listHexaSequenceFunction(modDis,modBas):
        #empty list to store cell coordinates
        listHexaSequence = []

        #definition of hexahedrons cell coordinates
        for lay in range(modDis['cellLays']):
            for row in range(modDis['cellRows']):
                for col in range(modDis['cellCols']):
                    if modBas['active'][0][row,col]==1:
                        pt0 = modDis['vertexPerLay']*(lay+1)+modDis['vertexCols']*(row+1)+col
                        pt1 = modDis['vertexPerLay']*(lay+1)+modDis['vertexCols']*(row+1)+col+1
                        pt2 = modDis['vertexPerLay']*(lay+1)+modDis['vertexCols']*(row)+col+1
                        pt3 = modDis['vertexPerLay']*(lay+1)+modDis['vertexCols']*(row)+col
                        pt4 = modDis['vertexPerLay']*(lay)+modDis['vertexCols']*(row+1)+col
                        pt5 = modDis['vertexPerLay']*(lay)+modDis['vertexCols']*(row+1)+col+1
                        pt6 = modDis['vertexPerLay']*(lay)+modDis['vertexCols']*(row)+col+1
                        pt7 = modDis['vertexPerLay']*(lay)+modDis['vertexCols']*(row)+col
                        anyList = [pt0,pt1,pt2,pt3,pt4,pt5,pt6,pt7]
                        listHexaSequence.append(anyList)
                    
        return listHexaSequence