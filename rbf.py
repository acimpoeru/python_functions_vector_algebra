
# coding: utf-8

# In[1]:

get_ipython().magic(u'pylab inline')
import pylab as pl
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D





# In[2]:

def plot_function1(my_list_of_points):
    for point in range(0,len(my_list_of_points)):
        #print list_points[point]
        xs = np.array([my_list_of_points[point][0]])
        ys = np.array([my_list_of_points[point][1]])
        zs = np.array([my_list_of_points[point][2]])
        plt.plot(xs,ys,marker='o', markersize=10,color='b',label='surface points')
def plot_function2(my_list_of_points):
    for point in range(0,len(my_list_of_points)):
        #print list_points[point]
        xs = np.array([my_list_of_points[point][0]])
        ys = np.array([my_list_of_points[point][1]])
        zs = np.array([my_list_of_points[point][2]])
        plt.plot(xs,ys,marker='o', markersize=10,color='r',label='line points')


# In[148]:

'''function to return the dot product, function to return the square distance, function to return the projection'''
'''Add function to return the distance between two points. Assumming point are distributed as tuples'''
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


'''define a translation'''

def translation(my_list_2):
    translation=[]
    for i in range(0,len(my_list_2)):
        translation.append((my_list_2[i][0],my_list_2[i][1]+i,my_list_2[i][2]))
        
    return translation

'''function to return a segment between two points'''

def segments(segment_points):
    
    segments=[]
    
    for i in range(0,len(segment_points)-1):
        
        segments.append([segment_points[i],segment_points[i+1]])
    return segments

def project(point,segment):
    paramater = ( (point[0]-segment[0][0])*(segment[1][0]-segment[0][0]) +
                 (point[1]-segment[0][1])*(segment[1][1]-segment[0][1]) +
                 (point[2]-segment[0][2])*(segment[1][2]-segment[0][2]) ) / (distance_points(segment[0],segment[1]))
    if paramater<0.0:
        point_projection = segment[0]
        return point_projection
    elif paramater>1.0:
        point_projection = segment[1]
        return point_projection
    else:
        x_p = segment[0][0] + paramater * (segment[1][0]-segment[0][0])
        y_p = segment[0][1] + paramater * (segment[1][1]-segment[0][1])
        z_p = segment[0][2] + paramater * (segment[1][2]-segment[0][2])
        return (x_p,y_p,z_p)     
    
'''function to calculate the solutions of the quadratics'''
def x_quadratic(interpolated_point,surf_point,projected_point):
    
    a = 1.0
    b = -2.0 * surf_point[0]

    c = (surf_point[0]**2.0 - projected_point[0]**2.0 + 2.0*projected_point[0]*interpolated_point[0] - 
         interpolated_point[0]**2.0)

    d = b**2-4*a*c 

    if d < 0.0:
        print 'No solutions'
    elif d == 0.0:
        x = -b / (2*a)
        print 'The sole solution is',x
    else: # if d > 0
        x1 = (-b + math.sqrt(d)) / (2*a)
        x2 = (-b - math.sqrt(d)) / (2*a)
        if x1 > x2 :
            x=x1
        else :
            x=x2
    return x
def y_quadratic(interpolated_point,surf_point,projected_point):
    
    a = 1.0
    b = -2.0 * surf_point[1]

    c = (surf_point[1]**2.0 - projected_point[1]**2.0 + 2.0*projected_point[1]*interpolated_point[1] - 
         interpolated_point[1]**2.0)

    d = b**2-4*a*c 

    if d < 0.0:
        print 'No solutions'
    elif d == 0.0:
        y = -b / (2*a)
        print 'The sole solution is',y
    else: # if d > 0
        y1 = (-b + math.sqrt(d)) / (2*a)
        y2 = (-b - math.sqrt(d)) / (2*a)
        if y1 > y2 :
            y=y1
        else :
            y=y2
    return y
def z_quadratic(interpolated_point,surf_point,projected_point):
    
    a = 1.0
    b = -2.0 * surf_point[2]

    c = (surf_point[2]**2.0 - projected_point[2]**2.0 + 2.0*projected_point[2]*interpolated_point[2] - 
         interpolated_point[2]**2.0)

    d = b**2-4*a*c 

    if d < 0.0:
        print 'No solutions'
    elif d == 0.0:
        z = -b / (2*a)
        print 'The sole solution is',z
    else: # if d > 0
        z1 = (-b + math.sqrt(d)) / (2*a)
        z2 = (-b - math.sqrt(d)) / (2*a)
        if z1 > z2 :
            z=z1
        else :
            z=z2
    return z


# In[35]:




# In[156]:


surface_points = [(-0.5,4.0,0.0),(1.5, 4.0, 0.0),(2.5,4.0,0.0),(3.5,4.0,0.0),(7.5,4.0,0.0)]

segment_points=[(0.0,0.0,0.0),(1.0,0.0,0.0),(2.0,0.0,0.0),(3.0,0.0,0.0),(4.0,0.0,0.0),(5.0,0.0,0.0)]
line_segments=segments(segment_points)
# apply the traslation 
demo_translation = translation(segment_points)
# transform the points into line segments
translation_demo_translation=segments(demo_translation)
#print demo_translation
print '                           '
#print translation_demo_translation
print '                           '
for surf_pt in surface_points:
    min_distance = 1e06
    for segment in line_segments:
        
        projected_point = project(surf_pt,segment)
        
        distance = distance_points(surf_pt,projected_point)
        if distance < min_distance:
            min_distance=distance
            
            index1 = line_segments.index(segment)
            
            min_seg = line_segments[index1]
            # search for the corresponding segment into the translated line segment
            translation_min_seg = translation_demo_translation[index1]
            
            index2 = surface_points.index(surf_pt)
            
            q = projected_point
            
    #print 'The dist is' ,min_distance,'from point C',surf_pt,'to point D',q,'for segment A'B'',min_seg
    
    #print 'The corresponding translated  line segment is AB =  ', translation_min_seg
    
    ''' Intersection of two segment is a system of simultaneous equations of the form Ax=B . x will contain the
    solution for the paramaters t and s . Numpy can solve that'''
    
    if surface_points.index(surf_pt)==0:
        print 'Deformation not applied'
        print 'deform by    ',0.0
        deformation = 0.0
        ''' work out the coordinates for the new point'''
        deform_point_coord = surface_points[surface_points.index(surf_pt)]
        print deform_point_coord
        
    elif surface_points.index(surf_pt)==len(surface_points)-1:
        print 'Found the end of the beam'
        print 'Apply the minimum distance to the nearest point '
        
        deformation = distance_points(segment_points[surface_points.index(surf_pt)+1],
                                      demo_translation[surface_points.index(surf_pt)+1])
        print 'Deform by translation ',deformation
        
        ''' work out the coordinates for the new point'''
        p = surf_pt
        q_point = segment_points[surface_points.index(surf_pt)+1]
        p_second = demo_translation[surface_points.index(surf_pt)+1]
        
        deform_x_coordinate_b = x_quadratic(p_second,p,q_point)
        deform_y_coordinate_b = y_quadratic(p_second,p,q_point)
        deform_z_coordinate_b = z_quadratic(p_second,p,q_point)
        
        new_deformed_point_b = (deform_x_coordinate_b,deform_y_coordinate_b,deform_z_coordinate_b)
        
        print 'The coord of the new point beam are ', new_deformed_point_b
        #deform_point_x_coord = 
        #deform_point_y_coord = 
        #deform_point_z_coord = 
        
        
    else:
        
        AB = translation_min_seg
        CD = [surf_pt,q]
        a_matrix = np.array([[(CD[1][0]-CD[0][0]),-(AB[1][0]-AB[0][0])],[(CD[1][1]-CD[0][1]),-(AB[1][1]-AB[0][1])]]) 
        b_matrix = np.array([[(AB[0][0]-CD[0][0])],[(AB[0][1]-CD[0][1])]])
        x = np.linalg.solve(a_matrix,b_matrix)
        s_paramater = x[0][0]
        t_parameter = x[1][0]
        print x
        ''' find the coordinates of the interpolated point '''
        x_interp = (AB[1][0]-AB[0][0]) * t_parameter +  AB[0][0]
        y_interp = (AB[1][1]-AB[0][1]) * t_parameter +  AB[0][1]
        z_interp = (AB[1][2]-AB[0][2]) * t_parameter +  AB[0][2]
    
        interp_point = (x_interp,y_interp,z_interp)
        print 'the interpoled point is  ', interp_point
        deformation = distance_points(interp_point,q)
        print 'deform by translation    ', deformation
        
        ''' work out the coordinates for the new point'''
        p_second = interp_point
        q_point = q
        p = surf_pt
        deform_x_coordinate = x_quadratic(p_second,p,q_point)
        deform_y_coordinate = y_quadratic(p_second,p,q_point)
        deform_z_coordinate = z_quadratic(p_second,p,q_point)
        
        new_deformed_point = (deform_x_coordinte,deform_y_coordinte,deform_z_coordinte)
        
        print   'The coord of the new point are  ',new_deformed_point
        
        

surface_points = [(-0.5,4.0,0.0),(1.5, 4.0, 0.0),(2.5,4.0,0.0),(3.5,4.0,0.0),(7.5,4.0,0.0)]
#points_in_space_tr = [(-0.5,4.0,0.0),(1.5, 4.5, 0.0),(2.5,5.0,0.0),(3.5,5.5,0.0),(8.0,6.0,0.0)]

segment_points=[(0.0,0.0,0.0),(1.0,0.0,0.0),(2.0,0.0,0.0),(3.0,0.0,0.0),(4.0,0.0,0.0),(5.0,0.0,0.0)]

deformed_points=[(1.5, 5.5, 0.0),(2.5, 6.5, 0.0),(3.5, 7.5, 0.0),(7.5, 9.0, 0.0)]
fig = plt.figure(figsize=(30, 15),dpi=150, facecolor='w', edgecolor='k')
axes = plt.gca()
axes.set_xlim([-1.0,10.0])
axes.set_ylim([-1.0,15.0])
#ax = fig.gca(projection='3d')
plot_function2(segment_points)
plot_function1(surface_points)
#plot_function1(points_in_space_tr)
my_translation=translation(segment_points)
plot_function2(my_translation)
plot_function1(deformed_points)
plt.legend()
plt.show()
#mapping(points_in_space,list_points)
