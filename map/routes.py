'''map.routes

This module contains the blueprint for maps
'''
import numpy as np
from flask import Blueprint, send_file
from io import BytesIO
from map.map import Map


map_bp = Blueprint('map', __name__, url_prefix='/map')

@map_bp.route('/get/<int:id>', defaults={'height': 1000, 'width': 1000})
@map_bp.route('/get/<int:id>/<int:height>/<int:width>')
def get(id, height, width):
    '''
    This route gets a map based on ID. The client can also pass in the
    width and height of the map in the GET parameters.

    Arguments:
        :param id: (int) ID of the Map (used as Seed)
        :param height: (int, default=1000) Height of image
        :param width: (int, default=1000) Width of image

    Returns:
        response: The map
    '''
    img_io = BytesIO()
    map = Map(id)
    image = map.draw_image(height, width)
    image.save(img_io, 'PNG', quality=100)
    img_io.seek(0)
    return send_file(img_io, mimetype='image/png', attachment_filename='{}.png'.format(id))


@map_bp.route('/get/random/', defaults={'height': 1000, 'width': 1000})
@map_bp.route('/get/random/<int:height>/<int:width>')
def random_map(height, width):
    '''
    This route returns a random map. The client can also pass in the
    width and height of the map in the GET parameters.

    Arguments:
        :param height: (int, default=1000) Height of image
        :param width: (int, default=1000) Width of image

    Returns:
        response: The map
    '''
    # Make Map ID
    id = np.random.randint(np.iinfo(np.int).max)
    print(id)

    # Now make map
    img_io = BytesIO()
    map = Map(id)
    image = map.draw_image(height, width)
    image.save(img_io, 'PNG', quality=100)
    img_io.seek(0)
    return send_file(img_io, mimetype='image/png', attachment_filename='{}.png'.format(id))