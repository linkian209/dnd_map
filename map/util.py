'''map.util

This module contains functions to help generate maps. It also contains
functions to display them.
'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d # noqa


# in_bounds and generate_voronoi create a bounded set of voronoi regions
# https://stackoverflow.com/a/33602171/6095609
# It works by mirroring the inital points across all 4 sides so the edges of
# each region create the bounding box as a border
def in_bounds(points, bounds):
    '''
    This helper function checks if points are within the inputted bounds.

    Arguments:
        :param points: (list) List of points to check
        :param bounds: (list) List of bounds to check against

    Returns:
        bool: Whether the points are in bounds or not.
    '''
    return np.logical_and(
        np.logical_and(
            bounds[0] <= points[:, 0],
            points[:, 0] <= bounds[1]
        ),
        np.logical_and(
            bounds[2] <= points[:, 1],
            points[:, 1] <= bounds[3]
        )
    )


def generate_voronoi(points, bounds, eps, create_mapping=False):
    '''
    This function generates a voronoi plot that is bounded by the
    inputted bounds
    '''
    # Get whether data is in bounds or not
    in_or_not = in_bounds(points, bounds)

    # Now we will mirror the initial points across all 4 sides so we will have
    # sets of points:
    # - Center: Initial set of points
    # - Left: Initial points reflected across the line x = 0.0
    # - Right: Initial points reflected across the line x = 1.0
    # - Upper: Initial points reflected across the line y = 1.0
    # - Lower: Initial points reflected across the line y = 0.0
    points_center = points[in_or_not, :]

    points_left = np.copy(points_center)
    points_left[:, 0] = bounds[0] - (points_left[:, 0] - bounds[0])

    points_right = np.copy(points_center)
    points_right[:, 0] = bounds[1] + (bounds[1] - points_right[:, 0])

    points_upper = np.copy(points_center)
    points_upper[:, 1] = bounds[3] + (bounds[3] - points_upper[:, 1])

    points_lower = np.copy(points_center)
    points_lower[:, 1] = bounds[2] - (points_lower[:, 1] - bounds[2])

    input_data = np.append(
        points_center,
        np.append(
            np.append(points_left, points_right, axis=0),
            np.append(points_upper, points_lower, axis=0),
            axis=0
        ),
        axis=0
    )

    # Compute Voronoi and then filter the regions to the bounded ones
    vor = Voronoi(input_data)
    regions = []
    hashed_regions = np.array([hash_list(x) for x in vor.regions])
    filtered_point_region = np.zeros(len(points_center))
    filtered_region_point = np.empty(len(vor.regions))
    for region in vor.regions:
        is_in_bounds = True
        for index in region:
            # This point is at infinity
            if index is -1:
                is_in_bounds = False
                break
            # Else continue on
            else:
                x = vor.vertices[index, 0]
                y = vor.vertices[index, 1]

                # Check if it is in bounds
                if(not((bounds[0] - eps <= x and x <= bounds[1] + eps)
                   and (bounds[2] - eps <= y and y <= bounds[3] + eps))):
                    # This point is not in bounds
                    is_in_bounds = False
                    break

        if region != [] and is_in_bounds:
            if create_mapping:
                # First add this to the look up lists
                # Find the region index
                region_index = np.where(hashed_regions == hash_list(region))
                region_index = region_index[0][0]

                # Now find the point index
                point_index = np.where(vor.point_region == region_index)
                point_index = point_index[0][0]

                # Now find this point in the filtered point list
                filtered_point_index = np.where(points_center == vor.points[point_index])
                filtered_point_index = filtered_point_index[0][0]

                # Now add the points
                filtered_point_region[filtered_point_index] = len(regions)
                filtered_region_point[len(regions)] = filtered_point_index

            regions.append(region)

    # Add the filtered data to the vor object
    vor.filtered_points = points_center
    vor.filtered_regions = regions
    vor.filtered_point_region = filtered_point_region
    vor.filtered_region_point = filtered_region_point[:len(regions)]

    # Done!
    return vor


def plot_filtered_voronoi(vor, bounds, delta_x, delta_y):
    # Check first if we have the filtered data
    if(not hasattr(vor, 'filtered_points')
       and not hasattr(vor, 'filtered_regions')):
        return

    axes = plt.gca()

    # Plot initial points
    axes.plot(vor.filtered_points[:, 0], vor.filtered_points[:, 1], 'b.')

    # Plot ridges and ridge points and centroids
    for region in vor.filtered_regions:
        # Ridge Points
        vertices = vor.vertices[region, :]
        axes.plot(vertices[:, 0], vertices[:, 1], 'go')

        # Ridges
        vertices = vor.vertices[region + [region[0]], :]
        axes.plot(vertices[:, 0], vertices[:, 1], 'k-')

        # Centroids
        centroid_region = centroid(vertices)
        axes.plot(centroid_region[0], centroid_region[1], 'r.')

    axes.set_xlim([bounds[0] - (delta_x/10), bounds[1] + (delta_x/10)])
    axes.set_ylim([bounds[2] - (delta_y/10), bounds[3] + (delta_y/10)])
    plt.show()


def plot_voronoi_delaunay(vor, delaunay, bounds, delta_x, delta_y):
    # Check first if we have the filtered data
    if(not hasattr(vor, 'filtered_points')
       and not hasattr(vor, 'filtered_regions')):
        return

    axes = plt.gca()

    # Plot initial points
    axes.plot(vor.filtered_points[:, 0], vor.filtered_points[:, 1], 'b.')

    # Plot ridges and ridge points and centroids
    for region in vor.filtered_regions:
        # Ridge Points
        vertices = vor.vertices[region, :]
        axes.plot(vertices[:, 0], vertices[:, 1], 'go')

        # Ridges
        vertices = vor.vertices[region + [region[0]], :]
        axes.plot(vertices[:, 0], vertices[:, 1], 'k-')

        # Centroids
        centroid_region = centroid(vertices)
        axes.plot(centroid_region[0], centroid_region[1], 'r.')

    # Now plot Delaunay triangulations
    plt.triplot(
        delaunay.points[:, 0], delaunay.points[:, 1], delaunay.simplices
    )

    axes.set_xlim([bounds[0] - (delta_x/10), bounds[1] + (delta_x/10)])
    axes.set_ylim([bounds[2] - (delta_y/10), bounds[3] + (delta_y/10)])
    plt.show()


def centroid(points, alt_return=False):
    # Find Signed Area
    A = 0

    for i in range(len(points)):
        i1 = (i+1) % len(points)
        A += (points[i][0] * points[i1][1]) - (points[i1][0] * points[i][1])
    A /= 2
    if(not A):
        return alt_return
    # Now find centroid points
    Cx = 0
    Cy = 0
    for i in range(len(points)):
        i1 = (i+1) % len(points)
        Cx += (
            (points[i][0] + points[i1][0]) *
            ((points[i][0] * points[i1][1]) - (points[i1][0] * points[i][1]))
        )
        Cy += (
            (points[i][1] + points[i1][1]) *
            ((points[i][0] * points[i1][1]) - (points[i1][0] * points[i][1]))
        )
    Cx /= (6 * A)
    Cy /= (6 * A)
    return np.array([Cx, Cy])


def get_slope_and_intercept(a, b):
    slope = (b[1]-a[1])/(b[0]-[a[0]])
    inter = a[1] - slope * a[0]
    return (slope, inter)


def hash_list(a):
    return hash(str(a))