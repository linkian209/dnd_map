# Ian Implementation of Voronoi
import itertools
import heapq
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, hypot

class Point():
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __eq__(self, other):
        retval = self.x is other.x
        return retval and (self.y is other.y)
    def __str__(self):
        return '({}, {})'.format(self.x, self.y)

    __repr__ = __str__


class Segment():
    def __init__(self, point):
        self.start = point
        self.end = None
        self.done = False

    def finish(self, point):
        self.end = point
        self.done = True

    def __str__(self):
        return '{} to {}'.format(self.start, self.end)

    __repr__ = __str__

class Event():
    # Event Types
    SITE_EVENT = 1
    CIRCLE_EVENT = 2

    def __init__(self, event_type, point=None, arc=None):
        self.event_type = event_type
        self.point = point
        self.arc = arc
        self.valid = True

    def __str__(self):
        return '{} {} Event with point {} and arc {}'.format(
            'Valid' if self.valid else 'Invalid',
            'Site' if self.event_type is Event.SITE_EVENT else 'Circle',
            self.point, self.arc
        )

    __repr__ = __str__

class Arc():
    def __init__(self, point, key):
        self.point = point
        self.key = key
        self.right = None
        self.left = None
        self.seg_left = None
        self.seg_right = None
        self.event = None

    def __str__(self):
        return 'Arc: key=[{}] Point={}'.format(self.key, self.point)

    __repr__ = __str__
        
            
class TreeNode:
    def __init__(self, key, val):
        self.key = key
        self.val = val
        self.left = None
        self.right = None 
        
    def __iter__(self):
        """ return the iterator that iterates through the elements in the BST 
        rooted at this node in an inorder sequence """
        
        if self.left:
            # The following iterates through all the nodes in the left subtree. 
            
            # The first thing that python does when the for loop is encountered
            # is to obtain an iterator object for the left subtree.  
            # This is done ("under the covers") by recursively calling 
            # the __iter__ method on the left child. 
            for elt in self.left:         
                yield elt
                
        # at this point we "visit" the current node
        yield (self.key, self.val)
        
        if self.right:
            # we now visit all the nodes in the right subtree 
            for elt in self.right:
                yield elt
            
    def put(self, key, val):
        """ add a new mapping between key and value in the tree """
        if self.key == key:
            self.val = val              # replace existing value
        elif self.key > key:            # key belongs in left subtree 
            if self.left:
                self.left.put(key,val)
            else:                       # left subtree is empty
                self.left = TreeNode(key,val)
        else:                           # key belongs in right subtree 
            if self.right:
                self.right.put(key,val)
            else:                       # right subtree is empty
                self.right = TreeNode(key,val)
                
    def get(self, key):
        """ get the value associated with the key """
        if self.key == key:
            return self.val
        
        if self.key > key:              # key should be in the left subtree
            if self.left:
                return self.left.get(key)
            else:
                return None
        else:                           # key should be in the right subtree
            if self.right:
                return self.right.get(key)
            else:
                return None
            
    def delete(self, key):
        """ delete the node with the given key and return the 
        root node of the tree """
            
        if self.key == key:
            # found the node we need to delete
            
            if self.right and self.left: 
                
                # get the successor node and its parent 
                [psucc, succ] = self.right._findMin(self)
                
                # splice out the successor
                # (we need the parent to do this) 
                
                if psucc.left == succ:
                    psucc.left = succ.right
                else:
                    psucc.right = succ.right
                                
                # reset the left and right children of the successor
                
                succ.left = self.left
                succ.right = self.right
                
                return succ                
                
            else:
                # "easier" case
                if self.left:
                    return self.left    # promote the left subtree
                else:
                    return self.right   # promote the right subtree 
        else:
            if self.key > key:          # key should be in the left subtree
                if self.left:
                    self.left = self.left.delete(key)
                # else the key is not in the tree 
                    
            else:                       # key should be in the right subtree
                if self.right:
                    self.right = self.right.delete(key)
                    
        return self
    
    def _findMin(self, parent):
        """ return the minimum node in the current tree and its parent """

        # we use an ugly trick: the parent node is passed in as an argument
        # so that eventually when the leftmost child is reached, the 
        # call can return both the parent to the successor and the successor
        
        if self.left:
            return self.left._findMin(self)
        else:
            return [parent, self]

class PriorityQueue():
    def __init__(self):
        self.events = []
        self.entry_finder = {}
        self.counter = itertools.count()

    def push(self, event):
        if(event in self.entry_finder):
            return

        count = next(self.counter)
        entry = [event.point.y, count, event]
        self.entry_finder[event] = event
        heapq.heappush(self.events, entry)

    def remove_event(self, event):
        entry = self.entry_finder.pop(event)
        entry[-1] = "REMOVED"

    def pop(self):
        while self.events:
            priority, count, event = heapq.heappop(self.events)
            if(event is not "REMOVED"):
                del self.entry_finder[event]
                return event

    def empty(self):
        return not self.events

    def __str__(self):
        return '{}'.format(self.events)

class IanVoronoi():
    def __init__(self,xmin,xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.events = PriorityQueue()
        self.arcs = None
        self.segments = []
        self.vertices = []

    def process(self, points):
        # First insert points into the event queue as site events
        for point in points:
            # Points are np arrays
            new_point = Point(point[0], point[1])
            self.events.push(Event(Event.SITE_EVENT,point=new_point))
        
        # Now process events
        while not self.events.empty():
            cur_event = self.events.pop()
            if(cur_event.event_type is Event.SITE_EVENT):
                self.process_site_event(cur_event)
            else:
                self.process_circle_event(cur_event)

        # Clean up!
        self.clean_up()

    def circle(self, a, b, c):
        # Find center of circumcircle
        # Computational Geometry in C (2nd ed) pg 189

        # Look to see if AB -> BC is goes right
        if(((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) > 0):
            return None

        A = b.x - a.x
        B = b.y - a.y
        C = c.x - a.x
        D = c.y - a.y
        E = A * (a.x + b.x) + B * (a.y + b.y)
        F = C * (a.x + c.x) + D * (a.y + c.y)
        G = 2 * (A * (c.y - b.y) - B * (c.x - b.x))

        # Points are colinear
        if(G is 0):
            return None

        circle_x = (D * E - B * F) / G
        circle_y = (A * F - C * E) / G

        return Point(circle_x, circle_y)

    def check_for_circle(self, arc, point):
        if(arc.event is not None and arc.event.point.y is not self.ymax):
            arc.event.valid = False
        arc.event = None

        #print('Arc: {}, Left: {}, Right: {}'.format(arc, arc.left, arc.right))

        if(arc.left is None or arc.right is None):
            return None
        
        if(arc.left.point == arc.right.point):
            return None

        # Get center of circumcircle
        circle = self.circle(arc.left.point, arc.point, arc.right.point)

        print(circle)

        # Now, find the radius of this circle
        r = hypot(arc.point.x - circle.x, arc.point.y - circle.y)

        # If the max event is heigher than the top, no event
        if(circle.y + r > self.ymax):
            return None
        # Else we have an event
        else:
            arc.event = Event(Event.CIRCLE_EVENT, circle, arc)
            self.events.push(arc.event)

    def intersection(self, left_arc, right_arc, point):
        # Use quadratic equation to find the intersection
        # y = 1/4f x^2 - Vx/2f x + Vx/4f + Vy
        # f = (Fy - y0) / 2
        # Vx = Fx
        # Vy = Fy - f
        f1 = (left_arc.point.y - point.y) / 2
        f2 = (right_arc.point.y - point.y) / 2

        a1 = 1 / (4 * f1)
        a2 = 1 / (4 * f2)
        A = a2 - a1

        b1 = -left_arc.point.x / (2 * f1)
        b2 = -right_arc.point.x / (2 * f2)
        B = b2 - b1

        c1 = left_arc.point.x**2 / (4 * f1) + left_arc.point.y - f1
        c2 = right_arc.point.x**2 / (4 * f2) + right_arc.point.y - f2
        C = c2 - c1

        # We only want the following X value
        x = (-B - sqrt(B**2 - 4 * A * C)) / (2 * A)
        y = a1 * x**2 + b1 * x + c1

        return Point(x, y)

    def intersect(self, arc, point):
        if(arc is None or point is None):
            return None

        # If the y values are the same, they do not intersect
        if(arc.point.y is point.y):
            return None

        right = None
        left = None
        # See if we can intersect between the left and right
        if(arc.right is not None):
            right = self.intersection(arc, arc.right, point)

        if(arc.left is not None):
            left = self.intersection(arc.left, arc, point)

        # Check if we intersect. We intersect if both are true:
        # 1)
        #   a) The left arc is None
        #   or
        #   b) The intersection of the left and current is to the right
        # and
        # 2)
        #   a) the right arc is None
        #   b) The intersection of the current and right is to the left
        if((arc.left is None or left.x <= point.x) 
           and (arc.right is None or point.x <= right.x)):
            point_x = point.x
            point_y = ((1 / (2 * (arc.point.y - point.y))) * point_x**2 
                    - (arc.point.x / (arc.point.y - point.y)) * point_x
                    + (arc.point.x**2 / (2 * (arc.point.y - point.y)))
                    + (arc.point.y + (arc.point.y - point.y) / 2))
            return Point(point_x, point_y)

        return None

    def insert_arc_at(self, arc, point):
        # Create New Arcs
        left_arc = Arc(arc.point, point.x-.000000001)
        new_arc = Arc(point, point.x)
        right_arc = Arc(arc.point, point.x+.000000001)

        # Continue existing segments
        right_arc.seg_right = arc.seg_right

        # Create new line segments
        left_seg = Segment(point)
        self.segments.append(left_seg)
        left_arc.seg_right = new_arc.seg_left = left_seg

        right_seg = Segment(point)
        self.segments.append(right_seg)
        new_arc.seg_right = right_arc.seg_left = right_seg

        # Do pointer stuff
        if(arc.left is not None):
            arc.left.right = left_arc
            left_arc.left = arc.left

        left_arc.right = new_arc
        new_arc.left = left_arc
        new_arc.right = right_arc
        right_arc.left = new_arc

        if(arc.right is not None):
            right_arc.right = arc.right
            arc.right.left = right_arc

        print('-------------------------------------------------------------')
        print('Arc: {}, Left: {}, Right: {}'.format(arc, arc.left, arc.right))
        print('Arc: {}, Left: {}, Right: {}'.format(self.arcs, self.arcs.left, self.arcs.right))
        print('Arc: {}, Left: {}, Right: {}'.format(new_arc, new_arc.left, new_arc.right))
        print('Arc: {}, Left: {}, Right: {}'.format(left_arc, left_arc.left, left_arc.right))
        print('Arc: {}, Left: {}, Right: {}'.format(right_arc, right_arc.left, right_arc.right))

        # Check for circle events
        self.check_for_circle(left_arc, point)
        self.check_for_circle(arc, point)
        self.check_for_circle(right_arc, point)

        return True

    def site_event(self, arc, point):
        # We go left as far as we can first
        if(point.x < arc.key):
            # Is the left open?
            if(arc.left is None):
                # Check if we intersect
                inter_point = self.intersect(arc, point)

                # Do we intersect here?
                if(inter_point is not None):
                    return self.insert_arc_at(arc, point)
                else:
                    return False
            # Keep going left!
            else:
                return self.site_event(arc.left, point)
        # Now go as far right as we can
        else:
            # Is the right open?
            if(arc.right is None):
                # Check if we intersect
                inter_point = self.intersect(arc, point)

                # Do we intersect here?
                if(inter_point is not None):
                    return self.insert_arc_at(arc, point)
                else:
                    return False
            # Should we check for an intersection?    
            elif(arc.right.key >= point.x):
                inter_point = self.intersect(arc, point)

                # Do we intersect here?
                if(inter_point is not None):
                    # If the intersection is to the right of us, we intersect
                    # the current arc 
                    if(point.x < inter_point.x):
                        return self.insert_arc_at(arc, point)
                    else:
                        return self.insert_arc_at(arc.right, point)
                else:
                    return False
            # Keep going right!
            else:
                return self.site_event(arc.right, point)

    def process_site_event(self, event):
        # If there are no arcs, add this on top
        if(self.arcs is None):
            self.arcs = Arc(event.point, event.point.x)
        # Else, loop through and see where this goes
        else:
            done = self.site_event(self.arcs, event.point)

            # Check if we intersected anywhere
            if(not done):
                # Add it to the list where it should go.
                arc = self.arcs
                new_arc = None
                start_x = None
                while(arc is not None):
                    if(event.point < arc.key):
                        if(arc.left is None):
                            arc.left = new_arc = Arc(event.point, event.point.x)
                            new_arc.right = arc
                            start_x = (arc.point.x + event.point.x) / 2
                            arc = None
                        else:
                            arc = arc.left
                    else:
                        if(arc.right is None):
                            arc.right = new_arc = Arc(event.point, event.point.x)
                            new_arc.left = arc
                            start_x = (arc.point.x + event.point.x) / 2
                            arc = None
                        else:
                            arc = arc.right
                # Add a line segment starting at the edge of the area,
                # halfway between the two
                seg = Segment(Point(start_x, self.ymax))
                self.segments.append(seg)
                if(new_arc.left is None):
                    new_arc.left.seg_right = new_arc.seg_left = seg
                else:
                    new_arc.right.seg_left = new_arc.seg_right = seg

    def process_circle_event(self, event):
        # First check if this event is valid, if not return
        if(not event.valid):
            return

        arc = event.arc

        if(arc.left is None or arc.right is None):
            return

        # Check if this event's arc's neighbors have events.
        # If so, invalidate them
        if(arc.left.event is not None):
            arc.left.event = False
        
        if(arc.right.event is not None):
            arc.right.event = False

        # Get center of circumcircle
        circle = self.circle(arc.left.point, arc.point, arc.right.point)
        self.vertices.append(circle)

        # New segment at the center of this circle
        seg = Segment(circle)
        self.segments.append(seg)
        arc.seg_left.finish(circle)
        arc.seg_right.finish(circle)

        arc.left.seg_right = seg
        arc.left.right = arc.right
        arc.right.seg_left = seg
        arc.right.left = arc.left

        # Check for new events
        self.check_for_circle(arc.left, circle)
        self.check_for_circle(arc.right, circle)

    def clean_up(self):
        # Loop through segments and finish them all
        for seg in self.segments:
            seg.done = True