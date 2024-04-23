import numpy as np
import math
import matplotlib.pyplot as plt

def generate_slab_geometry(thickness, length, gridwidth, radius,  dip, slab_depth):
    def slab_coordinates():
        num_segments = int(length / gridwidth)
        segment_length = length / num_segments

        coords = []

        for i in range(num_segments):
            x1 = i * segment_length
            y1 = 0
            x2 = (i + 1) * segment_length
            y2 = 0
            x3 = (i + 1) * segment_length
            y3 = -thickness
            x4 = i * segment_length
            y4 = -thickness

            rectangle_coords = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]
            coords.append(rectangle_coords)

        return coords

    def semicircle_coords():
        radiuscurvey = -radius
        radiuscurvex = length
        circlecoords = []
        connectpoint = []
        radians_gridangle = gridwidth/radius
        gridangle = math.degrees(radians_gridangle)
        num_segments_curve = int(dip / gridangle)

        for i in range(num_segments_curve):
            x4 = radiuscurvex + radius * np.sin(math.radians(i * gridangle))
            y4 = -(radius - radius * np.cos(math.radians(i * gridangle)))
            x1 = radiuscurvex + (radius - thickness) * np.sin(math.radians(i * gridangle))
            y1 = -(radius - (radius - thickness) * np.cos(math.radians(i * gridangle)))
            x2 = radiuscurvex + (radius - thickness) * np.sin(math.radians((i + 1) * gridangle))
            y2 = -(radius - (radius - thickness) * np.cos(math.radians((i + 1) * gridangle)))
            x3 = radiuscurvex + radius * np.sin(math.radians((i + 1) * gridangle))
            y3 = -(radius - radius * np.cos(math.radians((i + 1) * gridangle)))

            connectpoint.append((x2, y2))
            circle_rec = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]
            circlecoords.append(circle_rec)

        return connectpoint, circlecoords

    def rotateandshift_coordinates(connectpoint, coords, slab_dip):
        rotated_coordinates = []
        for x, y in coords:
            angle_rad = math.radians(slab_dip)
            new_x = x * math.cos(angle_rad) + y * math.sin(angle_rad) + connectpoint[-1][0]
            new_y = -x * math.sin(angle_rad) + y * math.cos(angle_rad) + connectpoint[-1][1]
            rotated_coordinates.append((new_x, new_y))

        return rotated_coordinates

    def dipping_slab_coordinates(connectpoint, thickness, gridwidth):
        slab2length = ((slab_depth + connectpoint[-1][1]) / np.sin(math.radians(dip)))
        num_segments_slab = int(slab2length / gridwidth)
        segment_length_slab = slab2length / num_segments_slab
        coords2 = []

        for i in range(num_segments_slab):
            x1 = i * segment_length_slab
            y1 = 0
            x2 = (i + 1) * segment_length_slab
            y2 = 0
            x3 = (i + 1) * segment_length_slab
            y3 = thickness
            x4 = i * segment_length_slab
            y4 = thickness

            rectangle_coords_slab = [(x1, y1), (x2, y2), (x3, y3), (x4, y4)]
            slab2_coords = rotateandshift_coordinates(connectpoint, rectangle_coords_slab, dip)
            coords2.append(slab2_coords)

        return coords2

    def plot_slab_coordinates(coords, coords2, circlecoords):
        for rectangle_coords in coords:
            x_coords = [coord[0] for coord in rectangle_coords]
            y_coords = [coord[1] for coord in rectangle_coords]
            plt.plot(x_coords, y_coords, 'b-')
        for rotated_coords_slab in coords2:
            x_coords_slab = [coord[0] for coord in rotated_coords_slab]
            y_coords_slab = [coord[1] for coord in rotated_coords_slab]
            plt.plot(x_coords_slab, y_coords_slab, 'b-')
        for circle_rec in circlecoords:
            x_coords_curve = [coord[0] for coord in circle_rec]
            y_coords_curve = [coord[1] for coord in circle_rec]
            plt.plot(x_coords_curve, y_coords_curve, 'b-')

        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (km)')
        plt.xlim(0, length)
        plt.ylim(-slab_depth, 10)  # Adjust the Y-axis limit accordingly
        plt.title('Slab Coordinates')
        plt.grid(False)
        plt.axis('equal')
        plt.show()

    # Main execution
    slab_coords = slab_coordinates()
    connectpoint, circlecoords = semicircle_coords()
    slab2_coords = dipping_slab_coordinates(connectpoint, thickness, gridwidth)

    # Plot the coordinates
    plot_slab_coordinates(slab_coords, slab2_coords, circlecoords)


# Test the combined function
generate_slab_geometry(thickness=50, length=800, gridwidth=20, radius=250, dip=70, slab_depth=500)
