''' FUNCTIONS TO READ THE DEFORMATION FROM FILE '''

def beam_configuration(input_file):
    
    beam_coordinates = {'x':[],'y':[],'z':[]}
    for line in input_file:
        coordinate = line.split()
        if coordinate[0] == 'NodeNo' and coordinate[1]=='X' and coordinate[2]=='Y' and coordinate[3]=='Z':
            print 'I cannot convert string to float'     
        else :
        
            beam_coordinates['y'].append(float(coordinate[2]))
            beam_coordinates['x'].append(float(coordinate[1]))
            beam_coordinates['z'].append(float(coordinate[3]))
            beam = []
            for i in range(0,len(beam_coordinates['x'])):
                for j in range(0,len(beam_coordinates['y'])):
                    for k in range(0,len(beam_coordinates['z'])):   
                        if i==j==k:
                            beam.append((beam_coordinates['x'][i],
                                         beam_coordinates['y'][j],
                                         beam_coordinates['z'][k]))
    
    return beam


def deformation_applied(input_file_2):

    deformation_translation =  {'dX':[],'dY':[],'dZ':[]}
    deformation_rotation =  {'ThetaX':[],'ThetaY':[],'ThetaZ':[]}

    for line in input_file_2:
        displacement = line.split()
    
        if displacement[0]=='NodeNo' and displacement[1]=='dX' and displacement[2]=='dY' and displacement[3]=='dZ'  and displacement[7]=='Thx' and displacement[8]=='Thy' and displacement[9]=='Thz' :
            print 'I cannot convert string to float'
        else:
            deformation_translation['dX'].append(float(displacement[1]))
            deformation_translation['dY'].append(float(displacement[2]))
            deformation_translation['dZ'].append(float(displacement[3]))
            translation_applied = []
            for i in range(0,len(deformation_translation['dX'])):
                for j in range(0,len(deformation_translation['dY'])):
                    for k in range(0,len(deformation_translation['dZ'])):
                        if i==j==k:
                            translation_applied.append((deformation_translation['dX'][i],
                                                        deformation_translation['dY'][j],
                                                        deformation_translation['dZ'][k]))
                        
                    
        
            deformation_rotation ['ThetaX'].append(float(displacement[7]))
            deformation_rotation ['ThetaY'].append(float(displacement[8]))
            deformation_rotation ['ThetaZ'].append(float(displacement[9]))
        
            rotation_applied=[]
        
            for thetax in range(0,len(deformation_rotation['ThetaX'])):
                for thetay in range(0,len(deformation_rotation['ThetaY'])):
                    for thetaz in range(0,len(deformation_rotation['ThetaZ'])):
                    
                        if thetax==thetay==thetaz:
                            rotation_applied.append((deformation_rotation['ThetaX'][thetax],
                                                     deformation_rotation['ThetaY'][thetay],
                                                     deformation_rotation['ThetaZ'][thetaz]))
                        
    return [translation_applied,rotation_applied]

