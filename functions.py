#!/usr/bin/python3
import numpy as np
import math
import matplotlib.pyplot as plt


#Create the flat/upper slab coordinates, split into segments of gridwidth
def slab_coordinates(thickness, length, gridwidth):
    num_segments = int(length / gridwidth)  # Calculate the number of segments
    print(num_segments)
    segment_length = length / num_segments  # Calculate the length of each segment

    coords = []  # List to store the coordinates

    for i in range(num_segments):
        x1 = i * segment_length         # Calculate the x-coordinate of the bottom left corner
        y1 = 0                          # y-coordinate of the bottom left corner
        x2 = (i + 1) * segment_length   # x-coordinate of the bottom right corner
        y2 = 0                          # y-coordinate of the bottom right corner
        x3 = (i + 1) * segment_length   # x-coordinate of the top right corner
        y3 = -thickness                 # y-coordinate of the top right corner
        x4 = i * segment_length         # x-coordinate of the top left corner
        y4 = -thickness                 # y-coordinate of the top left corner

        rectangle_coords = [(x1, y1), (x2, y2), (x3, y3), (x4, y4),(x1, y1)]  # Coordinates of the current rectangle
        coords.append(rectangle_coords)                              # Add the coordinates to the list

    return coords


#The below function is currently super low resolution, needs to be fixed for circle plotting and arclength calculation
#Create the curve coordinates, split into segments of gridangle and stopping at a given dip relative to the surface
def semicircle_coords(thickness, length, radius, gridangle, dip):
    #Y Value for radius of curvature center
    radiuscurvey = -radius
    #X Value for radius of curvature center
    radiuscurvex = length
    #Circumference segment in radians


    circlecoords = []
    connectpoint = []
    num_segments_curve = int(dip/gridangle)
    
    for i in range(num_segments_curve):
        x4 = radiuscurvex + radius * np.sin(math.radians(i * gridangle))  # same order: 1,2,3,4:bottomleft,bottomright,topright,topleft
        y4 = -(radius - radius * np.cos(math.radians(i * gridangle)))  # 
        x1 = radiuscurvex + (radius - thickness) * np.sin(math.radians(i * gridangle))
        y1 = -(radius - (radius - thickness) * np.cos(math.radians(i * gridangle)))
        x2 = radiuscurvex + (radius - thickness) * np.sin(math.radians((i + 1) * gridangle))
        y2 = -(radius - (radius - thickness) * np.cos(math.radians((i + 1) * gridangle)))
        x3 = radiuscurvex + radius * np.sin(math.radians((i + 1) * gridangle))
        y3 = -(radius - radius * np.cos(math.radians((i + 1) * gridangle)))
         


        connectpoint.append((x2,y2))
        circle_rec = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]
        circlecoords.append(circle_rec)


    return  connectpoint, circlecoords
    
#Function to rotate and shift the dipping slab coordinates
def rotateandshift_coordinates(connectpoint, coords, slab_dip):
        rotated_coordinates = []
        for x, y in coords:
            # Convert angle to radians
            angle_rad = math.radians(slab_dip)
            
            # Apply rotation formula
            new_x = x * math.cos(angle_rad) + y * math.sin(angle_rad) +  connectpoint[len(connectpoint)-1][0]
            new_y = -x * math.sin(angle_rad) + y * math.cos(angle_rad) +  connectpoint[len(connectpoint)-1][1]
            
            rotated_coordinates.append((new_x, new_y))
        
        return rotated_coordinates
    


#Create the dipping slab coordinates
def dipping_slab_coordinates(connectpoint, thickness, slab_depth,slab_dip,gridwidth):
        slab2length = ((slab_depth +  connectpoint[len(connectpoint)-1][1]) / np.sin(math.radians(slab_dip)))

        
        num_segments_slab = int(slab2length / gridwidth)  # Calculate the number of segments
    
        segment_length_slab = slab2length / num_segments_slab  # Calculate the length of each segment

        coords2 = []  # List to store the coordinates

        for i in range(num_segments_slab):
          x1 = i * segment_length_slab        # Calculate the x-coordinate of the bottom left corner
          y1 = 0                          # y-coordinate of the bottom left corner
          x2 = (i+1) *  segment_length_slab   # x-coordinate of the bottom right corner
          y2 = 0                           # y-coordinate of the bottom right corner
          x3 = (i+1) *  segment_length_slab       # x-coordinate of the top right corner
          y3 = thickness                 # y-coordinate of the top right corner
          x4 = i * segment_length_slab         # x-coordinate of the top left corner
          y4 = thickness                 # y-coordinate of the top left corner

          rectangle_coords_slab = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]  # Coordinates of the current rectangle
          slab2_coords = rotateandshift_coordinates(connectpoint, rectangle_coords_slab, slab_dip)
          coords2.append(slab2_coords)                              # Add the coordinates to the list
          
        return coords2
       

        
#Plot the coordinates of the initial geometry
def plot_slab_coordinates(coords, coords2,  circlecoords):
    plt.figure(figsize=(10, 6))  # Set the figure size
    # plt.style.use('seaborn-darkgrid')  # Use seaborn-darkgrid style
    # plt.style.use('seaborn-darkgrid')  # Use seaborn-darkgrid style
    for rectangle_coords in coords:
        x_coords = [coord[0] for coord in rectangle_coords]
        y_coords = [coord[1] for coord in rectangle_coords]
        plt.plot(x_coords, y_coords, 'b-', linewidth=2)  # Increase line width

    for rotated_coords_slab in coords2:
        x_coords_slab = [coord[0] for coord in rotated_coords_slab]
        y_coords_slab = [coord[1] for coord in rotated_coords_slab]
        plt.plot(x_coords_slab, y_coords_slab, 'b-', linewidth=2)  # Increase line width

    for circle_rec in circlecoords:
        x_coords_curve = [coord[0] for coord in circle_rec]
        y_coords_curve = [coord[1] for coord in circle_rec]
        plt.plot(x_coords_curve, y_coords_curve, 'b-', linewidth=2)  # Increase line width

    plt.xlabel('Distance (km)', fontsize=14)  # Increase font size
    plt.ylabel('Depth (km)', fontsize=14)  # Increase font size
    plt.xlim(0, 400)
    plt.ylim(-600, 10)  # Modified to plot Y axis downwards
    plt.title('Slab Coordinates', fontsize=16)  # Increase font size
    plt.grid(True)  # Show grid
    plt.axis('equal')
    plt.tick_params(labelsize=12)  # Increase tick label size

    plt.show()


