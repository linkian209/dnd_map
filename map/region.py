'''map.region

This module contains the region class.
'''
import numpy as np


class Region():
    '''
    This class holds a region of the map. It contains biome information as well
    as the points making up the region. Each region is drawn differently based
    on its biome.

    Attributes:
        :param biome: (??) Biome type
        :param points: (list) Points bounding the region
    '''
    def __init__(self, id, points, center, biome):
        '''
        This method initializes the Region object with the inputted biome
        type and points.

        Arguments:
            :param self: (Region) This region
            :param id: (int) The id of this region, used for seeding
            :param points: (list) Points bounding this region
            :param center: (list) Centroid of the region
            :param biome: (???) The biome type
        '''
        self.id = id
        self.points = points
        self.tuple_points = [(x[0], x[1]) for x in points]
        self.biome = biome

    def draw(self, im_draw, height, width):
        '''
        This method draws the Region based on its biome. It will scale it based
        on the inputted height and width

        Arguments:
            :param self: (Region) This region object
            :param im_draw: (PIL.ImageDraw) The drawer object
            :param height: (int) Height of the final image
            :param width: (int) Width of the final image
        '''
        points = [(x[0] * width, x[1] * height) for x in self.tuple_points]
        color = 'rgb({},{},{})'.format(
            np.random.randint(0, 256), np.random.randint(0, 256),
            np.random.randint(0, 256)
        )
        im_draw.polygon(points, fill=color)
