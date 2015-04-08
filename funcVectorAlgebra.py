import math 
import os
import sys
''' FUNCTION TO RETURN THE DISTANCE BETWEEN TWO POINTS'''

def distance_points(tuple1,tuple2):
    distance=math.sqrt((tuple2[0]-tuple1[0])**2.0+
                       (tuple2[1]-tuple1[1])**2.0+
                       (tuple2[2]-tuple1[2])**2.0)
    return distance



'''function to return the legth of the segments'''
def squared_distance(my_list_of_points):
    squared_distance=np.zeros(len(my_list_of_points)-1)
    for i in range(0,len(squared_distance)):
        squared_distance[i]=math.sqrt((my_list_of_points[i+1][0]-my_list_of_points[i][0])**2.0 +
                                      (my_list_of_points[i+1][1]-my_list_of_points[i][1])**2.0 +
                                      (my_list_of_points[i+1][2]-my_list_of_points[i][2])**2.0)
        
    return squared_distance  


'''function to return a vector'''
def construct_vector(my_list_of_points):
    construct_vector={}
    counter=0
    for i in range(0,len(my_list_of_points)-1):
        construct_vector[i]=(my_list_of_points[i+1][0]-my_list_of_points[i][0],
                             my_list_of_points[i+1][1]-my_list_of_points[i][1],
                             my_list_of_points[i+1][2]-my_list_of_points[i][2])
    
    return construct_vector


'''function to construct a vector'''
def construct_vector_PtoS(my_list_1,my_list_2):
    '''loop over the main dictionaries: output two entries'''
    '''my_list_1 is the list of points in space and my_list_2 is the list_forming the segment'''
    construct_vector_PtoS={}
    counter = 0
    for i in range(0,len(my_list_1)):
        construct_vector_PtoS[counter]=[]
        for j in range(0,len(my_list_2)):
            #construct_vector_PtoS[counter_1]=[]
            construct_vector_PtoS[counter].append((my_list_1[i][0]-my_list_2[j][0],
                                                   my_list_1[i][1]-my_list_2[j][1],
                                                   my_list_1[i][2]-my_list_2[j][2]))
        counter+=1
    #print 'the dictionary is = ', construct_vector_PtoS
    
    return construct_vector_PtoS



'''define a function to calculate distances  from points to the segments'''
def dist_and_neighbours(my_list_1,my_list_2):
    dist_and_neighbours={}
    count=0
    for i in range(0,len(my_list_1)):
        dist_and_neighbours[count]=[]
        for j in range(0,len(my_list_2)):
           
            dist_and_neighbours[count].append(math.sqrt((my_list_1[i][0]-my_list_2[j][0])**2.0 +
                                                        (my_list_1[i][1]-my_list_2[j][1])**2.0 +
                                                        (my_list_1[i][2]-my_list_2[j][2])**2.0))
            
        count+=1
    #print dist     
    return dist_and_neighbours



'''define function to find the parameter t '''
def mapping(my_list_1,my_list_2):
    parameter_t={}
    count=0
    distance=squared_distance(my_list_2)
    
    
    for i in range(0,len(my_list_1)):
        parameter_t[count]=[]
        for j in range(0,len(my_list_2)-1):
            parameter_t[count].append((((my_list_1[i][0]-my_list_2[j][0])*(my_list_2[j+1][0]-my_list_2[j][0]))+
                                      ((my_list_1[i][1]-my_list_2[j][1])*(my_list_2[j+1][1]-my_list_2[j][1]))+
                                      ((my_list_1[i][2]-my_list_2[j][2])*(my_list_2[j+1][2]-my_list_2[j][2])))/
                                      (distance[j]**2.0))
        count+=1
    distances_and_neighbours=dist_and_neighbours(my_list_1,my_list_2)
        

    '''store the projections '''
    projections={}
    nearest_points={}
    projected_points={}
    
    count_1=0
    count_2=0
    count_3=0
    
    for key in parameter_t:
        projections[count_1]=[]
        nearest_points[count_2]=[]
        projected_points[count_3]=[]
        for t in range(0,len(parameter_t[key])):
            
            if parameter_t[key][t]<0.0:
                
                '''It falls BEFORE LEFT point of the segment. return the point '''
                '''(my_list_2[t]),my_list_1[key]'''
                projections[count_1].append(distances_and_neighbours[key][t])
                nearest_points[count_2].append(my_list_2[t])
                projected_points[count_3].append(my_list_1[key])
            elif parameter_t[key][t]>1.0:
                
                '''It falls AFTER RIGHT point of the segment. return the point '''
                '''(my_list_2[t+1]),my_list_1[key]'''
                projections[count_1].append(distances_and_neighbours[key][t+1])
                nearest_points[count_2].append(my_list_2[t+1])
                projected_points[count_3].append(my_list_1[key])
                    
            else:
                mypoint_q = ()
                mypoint_q= mypoint_q + (my_list_2[t][0]+parameter_t[key][t]*(my_list_2[t+1][0]-my_list_2[t][0]),
                                        my_list_2[t][1]+parameter_t[key][t]*(my_list_2[t+1][1]-my_list_2[t][1]),
                                        my_list_2[t][2]+parameter_t[key][t]*(my_list_2[t+1][2]-my_list_2[t][2])) 
                #print 'succesfull projected points     ', mypoint_q
                '''IT FALLS ON THE SEGMENT '''
                '''(mypoint_q),my_list_1[key]'''
                projections[count_1].append(distance_points(mypoint_q,my_list_1[key]))
                nearest_points[count_2].append(mypoint_q)
                projected_points[count_3].append(my_list_1[key])
                
        count_1+=1
        count_2+=1
        count_3+=1
                
    return projections,nearest_points,projected_points



'''function to return a segment between two points'''

def segments(segment_points):
    
    segments=[]
    
    for i in range(0,len(segment_points)-1):
        
        segments.append([segment_points[i],segment_points[i+1]])
    return segments




def project(point,segment):
    #global parameter
    paramater = ( (point[0]-segment[0][0])*(segment[1][0]-segment[0][0]) +
                  (point[1]-segment[0][1])*(segment[1][1]-segment[0][1]) +
                  (point[2]-segment[0][2])*(segment[1][2]-segment[0][2]) ) / ((distance_points(segment[0],segment[1]))**2.0)

    #print 'I am the calculated parameter     ' , parameter

    if paramater<0.0:
        #print '`I am am before the segment',paramater
        point_projection = segment[0]
        return point_projection

    elif paramater>1.0:

        #print 'I am after the segment',paramater
        point_projection = segment[1]
        return point_projection

    else:
        #print 'between the end points of the segment',paramater
        x_p = segment[0][0] + paramater * (segment[1][0]-segment[0][0])
        y_p = segment[0][1] + paramater * (segment[1][1]-segment[0][1])
        z_p = segment[0][2] + paramater * (segment[1][2]-segment[0][2])
        return (x_p,y_p,z_p)     
    
'''function to calculate the coordinates of the moved point in respect with the distance from the interpolated_point to the  projected point'''
def x_interpolation(interpolated_point,surf_point,projected_point):
    
    x_coordinate = surf_point[0] - projected_point[0] + interpolated_point[0]

    return x_coordinate

def y_interpolation(interpolated_point,surf_point,projected_point):

    y_coordinate = surf_point[1] - projected_point[1] + interpolated_point[1]

    return  y_coordinate

def z_interpolation(interpolated_point,surf_point,projected_point):
    
    z_coordinate = surf_point[2] - projected_point[2] + interpolated_point[2]
    
    return z_coordinate

'''function to rotate a point in 3D'''

def pointRotate3D(p1, p2, p0, theta):
    from math import cos, sin, sqrt
     
    # Translate so axis is at origin  
    p=[0.0,0.0,0.0]
    p[0]= p0[0] - p1[0]
    p[1]= p0[1] - p1[1]
    p[2]= p0[2] - p1[2]
    
    p = (p[0] , p[1] , p[2])
    # Initialize point q
    q = [0.0,0.0,0.0]
    
    N = [0.0,0.0,0.0]
    
    N[0] = p2[0]-p1[0]
    N[1] = p2[1]-p1[1]
    N[2] = p2[2]-p1[2]
    
    N = (N[0] , N[1] , N[2])
    
    # calculate the distance in order to normalize
    Nm = sqrt(N[0]**2 + N[1]**2 + N[2]**2)
    
    # Rotation axis unit vector
    n=[0.0,0.0,0.0]
    n[0] = N[0]/Nm 
    n[1] = N[1]/Nm 
    n[2] = N[2]/Nm
    
    n = (n[0],n[1],n[2])
        
    # Matrix common factors     
    c = cos(theta)
    t = (1 - cos(theta))
    s = sin(theta)
    X = n[0]
    Y = n[1]
    Z = n[2]

    # Matrix 'M'
    d11 = t*X**2 + c
    d12 = t*X*Y - s*Z
    d13 = t*X*Z + s*Y
    d21 = t*X*Y + s*Z
    d22 = t*Y**2 + c
    d23 = t*Y*Z - s*X
    d31 = t*X*Z - s*Y
    d32 = t*Y*Z + s*X
    d33 = t*Z**2 + c

    #            |p.x|
    # Matrix 'M'*|p.y|
    #            |p.z|
    q[0] = d11*p[0] + d12*p[1] + d13*p[2]
    q[1] = d21*p[0] + d22*p[1] + d23*p[2]
    q[2] = d31*p[0] + d32*p[1] + d33*p[2]
    new_point=[0.0,0.0,0.0]
    new_point[0] = q[0] + p1[0]
    new_point[1] = q[1] + p1[1]
    new_point[2] = q[2] + p1[2]
    new_point = (new_point[0],new_point[1],new_point[2])
    # Translate axis and rotated point back to original location    
    return new_point

''' function to covert from degrees to radians'''

def degToRad(angle_degrees):
    
    angle_radians = angle_degrees * 2.0*math.pi/360.0
    return angle_radians


'''function to add halos to the beam sgment in order to cover all the suface points'''


def add_halo(beam_segments):
    
    left_factor = 3.0
    right_factor = 3.0

    first_segment = beam_segments[0]
    last_segment =  beam_segments[-1]
    #print first_segment
    #print last_segment
    #xq = xp + sqrt(5)(xp-xa)
    
    halo_left_x = first_segment[0][0] - left_factor * (first_segment[1][0]-first_segment[0][0])
    halo_left_y = first_segment[0][1] - left_factor * (first_segment[1][1]-first_segment[0][1])
    halo_left_z = first_segment[0][2] - left_factor * (first_segment[1][2]-first_segment[0][2])
    halo_left = (halo_left_x,halo_left_y,halo_left_z)
    #print 'lest halo',halo_left
    halo_right_x = last_segment[1][0] + right_factor * (last_segment[1][0]-last_segment[0][0])
    halo_right_y = last_segment[1][1] + right_factor * (last_segment[1][1]-last_segment[0][1])
    halo_right_z = last_segment[1][2] + right_factor * (last_segment[1][2]-last_segment[0][2])
    
    halo_right = (halo_right_x,halo_right_y,halo_right_z)
    #print 'right halo', halo_right
    
    
    beam_segments.insert(0,[halo_left,first_segment[0]])
    beam_segments.insert(len(beam_segments),[last_segment[1],halo_right])
    #print beam_geometry
    return beam_segments


def create_surface_points(beam):
    surface_points = []

    for i in range (1,len(beam)-1):

        surface_points.append((beam[i][0]+0.01 ,beam[i][1],beam[i][2]+0.5,))


    return surface_points


def theta_nodes(theta_list_of_angles):
    theta_segments = []

    for i in range(0,len(theta_list_of_angles)-1):

        theta_segments.append([theta_list_of_angles[i],theta_list_of_angles[i+1]])
    return theta_segments


def theta_interpolation(theta_segment,point,segment):

    theta_paramater_x = (point[0]-segment[0][0])/(segment[1][0]-segment[0][0])
    #print 'They are all equal'
    if theta_paramater_x<0.0 :
        #print 'The parameter is select the left node' ,theta_paramater_x
        theta_interp = theta_segment[0]
        return theta_interp

    elif theta_paramater_x>1.0 :  
            #print 'The parameter is , selenct the right node' ,theta_paramater_x
            theta_interp = theta_segment[1]
            return theta_interp

    else:

        #print 'between the end points of the segment',paramater
            #print 'The parameter is ' ,theta_paramater_x
        theta_interp = (1.0 - theta_paramater_x) * theta_segment[0] + theta_paramater_x * theta_segment[1]
        
        return theta_interp      
        




