'''map.map

This module contains the classes and enumerations for a map.
'''
import matplotlib.pyplot as plt
import numpy as np
import sys
from enum import Enum
from PIL import Image, ImageDraw
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d, delaunay_plot_2d # noqa
from map.util import centroid, generate_voronoi, plot_filtered_voronoi, plot_voronoi_delaunay # noqa
from map.region import Region


# Map Type Enumeration
class MapType(Enum):
    ISLAND = 1
    ISLANDS = 2

class Map():
    '''
    This class creates a map.

    Attributes:
        :param vor: (scipy.spactial.Voronoi) The Voronoi object containing
                    the map
        :param SEED: (int) The seed to generate the map

    '''
    # START FUNCTIONS
    def __init__(self, seed, relax_times=7, show_plots=False):
        '''
        This method initializes the Map class then procedurally generates a
        map using the inputted seed. Optionally, the number of relaxation
        iterations can be passed in as well as the option to show plots for
        debugging.

        Arguments:
            :param self: (Map) This object
            :param seed: (int) Seed for the random number generator
            :param relax_times: (int, default=7) Number of times to run Lloyd
                                algorithm
            :param show_plots: (bool, default=False) Show plots
        '''
        self.seed = seed
        self.relaxation_times = relax_times
        self.x_min = 0.0
        self.x_max = 1.0
        self.delta_x = self.x_max - self.x_min
        self.y_min = 0.0
        self.y_max = 1.0
        self.delta_y = self.y_max - self.y_min
        self.eps = sys.float_info.epsilon
        self.show_plots = show_plots

        np.random.seed(self.seed)
        self.num_points = np.random.randint(2000, 4000)
        self.data = np.random.rand(self.num_points, 2)
        self.bounds = [self.x_min, self.x_max, self.y_min, self.y_max]

        self.vor = generate_voronoi(self.data, self.bounds, self.eps)

        if(self.show_plots):
            voronoi_plot_2d(self.vor)
            plt.show()
            plot_filtered_voronoi(
                self.vor, self.bounds, self.delta_x, self.delta_y
            )

        # Now lets do some relaxation
        for i in range(self.relaxation_times):
            # Set up some variables
            centroids = np.empty_like(self.data)
            index = 0

            # Find the centroids of the region
            for region in self.vor.filtered_regions:
                vertices = self.vor.vertices[region + [region[0]], :]
                centroids[index] = centroid(vertices)
                index += 1
            
            # We only need to make the mappings for the last relaxation
            create_mapping = True if i+1 is self.relaxation_times else False

            # Create the Voronoi
            self.vor = generate_voronoi(centroids, self.bounds, self.eps, create_mapping)

            # Display plots if we need to
            if(self.show_plots):
                plot_filtered_voronoi(
                    self.vor, self.bounds, self.delta_x, self.delta_y
                )

        # After we generate the relaxed Voronoi plot, make the Delaunay
        # Triangulation
        self.delaunay = Delaunay(self.vor.filtered_points)

        if(self.show_plots):
            plot_voronoi_delaunay(
                self.vor, self.delaunay, self.bounds,
                self.delta_x, self.delta_y
            )

        # Now lets get the map type
        self.map_type = MapType(np.random.randint(1, len(MapType)))

        # Now generate the regions
        self.generate_regions()

    def generate_regions(self):
        '''
        This method takes the regions from the Voronoi object and converts
        them to Region objects and assigns a biome to each.
        '''
        self.regions = []
        region_index = 0
        for region in self.vor.filtered_regions:
            # Determine the biome of the region
            biome = np.random.randint(0, 10)

            # Next generate the ID of the region
            id = int('{}{}'.format(self.seed, len(self.regions)))

            # Now generate the region
            self.regions.append(
                Region(
                    id, self.vor.vertices[region, :], 
                    self.vor.filtered_region_point[region_index], biome
                )
            )

            # Increment index
            region_index += 0

    def draw_image(self, height, width, draw_outline=True):
        '''
        This method draws a PIL image of the map in the passed in resolution,
        then returns a PIL object.

        Arguments:
            :param self: (Map) This object
            :param height: (int) Height of final image in pixels
            :param width: (int) Width of final image in pixels

        Returns:
            PIL.Image: A PIL image object
        '''
        # Initialize the return value
        image = Image.new('RGB', (width, height))
        im_draw = ImageDraw.Draw(image)

        # Loop through the regions and draw them
        for region in self.regions:
            region.draw(im_draw, height, width)

        # Draw the outline as needed
        if(draw_outline):
            im_draw.rectangle(
                [0, 0, width - 1, height - 1],
                outline='#000000', width=1
            )

        # Done!
        return image
