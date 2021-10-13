"""
Filter OSM file for Ways and Relations with given tags, then generate, simplify and export as ArduPilot Compatable Geo-Fence waypoint files
"""

import osmium as o
import numpy as np
import os
import time
import math

input_file = 'switzerland-padded.osm.pbf'
#input_file = 'wales-latest.osm.pbf'
#input_file = 'britain-and-ireland-latest.osm.pbf'
#input_file = 'europe-latest.osm.pbf'
#input_file = 'bigger_map.osm'
#input_file = 'Heidsee.osm'

# output directory
directory = 'Fences'

# search tags
tags = (('landuse', 'reservoir', None),
        ('natural', 'water', ('lake', 'reservoir', 'basin', 'lagoon', 'pond')))



# Only output fences with outer area larger than this
area_threshold = 1000 # m^2

# simplication area removal treshold, set None to disable
simp_area_threshold = 100 # m^2

# don't simplify to less than this number of nodes
simp_min_nodes = 50
simp_max_nodes = 250

# check if given tags are in list of desired
def check_tags(test_tags, desired_tags):
    for tag in desired_tags:
        if tag[0] in test_tags and test_tags[tag[0]] == tag[1]:
            return_tag = tag[1]
            if (tag[2] != None) and (tag[1] in test_tags):
                if test_tags[tag[1]] not in tag[2]:
                    break
                return_tag = test_tags[tag[1]]
            # found tag
            return return_tag
    return None

# find all ways or multipolygons with the given tags
class fence_search(o.SimpleHandler):

    def __init__(self):
        super(fence_search, self).__init__()

    def area(self, a):
        found_tag = check_tags(a.tags, tags)
        if found_tag == None:
            return
        name = 'unnamed'
        if 'name' in a.tags:
            name = a.tags['name']
        for outer in a.outer_rings():
            x = [None]
            y = [None]
            x[0], y[0], origin = get_polygon(outer, None)
            if x is None:
                # polygon not valid
                continue
            area = polygon_area(x[0], y[0])
            if area < area_threshold:
                # too small
                continue
            # add inner polygons
            for inner in a.inner_rings(outer):
                x_temp, y_temp, _ = get_polygon(inner, origin)
                if x_temp is None:
                    continue
                if polygon_polygon_intersection([x[0],x_temp], [y[0],y_temp]):
                    continue
                x.append(x_temp)
                y.append(y_temp)

            # simplify to the given thresholds, possibly turing polygons in to circles
            x, y, radius = simplify_poly(x, y)

            # convert back to lat lon and get center for file name, avoid odd wrap issues by doing in cartesian
            center = [0, 0]
            num_poly = len(x)
            num_nodes = 0
            lat = [None] * num_poly
            lon = [None] * num_poly
            for i in range(num_poly):
                center[0] += np.mean(x[i])
                center[1] += np.mean(y[i])
                lat[i], lon[i] = convert_from_cartesian(x[i], y[i], origin[0], origin[1])
                num_nodes += len(lat[i])
            center[0] /= num_poly
            center[1] /= num_poly
            center = convert_from_cartesian([center[0]], [center[1]], origin[0], origin[1])

            # create waypoint file
            file_name = save_to_file(name, found_tag, center, lat, lon, radius)

            # add to js index file, outer polygon only
            js_file.write('{\n')
            js_file.write('name: "%s",\n' % (name.replace('"', '\\"')))
            js_file.write('num_nodes: %i,\n' % (num_nodes))
            js_file.write('area: %f,\n' % (area))
            js_file.write('file_name: "%s",\n' % (file_name.replace('"', '\\"')))
            js_file.write('nodes: [')
            for i in range(len(lat[0])):
                js_file.write('[%f, %f],\n' % (lat[0][i], lon[0][i]))
            js_file.write(']},\n')

# get and check polygon from osmium structures
def get_polygon(nodes, origin):
    num_nodes = len(nodes)
    outer_lat = np.empty(num_nodes)
    outer_lon = np.empty(num_nodes)
    for i in range(num_nodes):
        if nodes[i].location.valid():
            outer_lat[i] = nodes[i].location.lat
            outer_lon[i] = nodes[i].location.lon
        else:
            outer_lat[i] = np.NaN
            outer_lon[i] = np.NaN
    # remove nan values
    outer_lat = outer_lat[~np.isnan(outer_lat)]
    outer_lon = outer_lon[~np.isnan(outer_lon)]
    # make sure end point are not the same
    if (outer_lat[0] == outer_lat[-1]) and (outer_lon[0] == outer_lon[-1]):
        outer_lat = outer_lat[0:-1]
        outer_lon = outer_lon[0:-1]

    if len(outer_lat) < 3:
        # not enough points
        return None, None, None
    if origin is None:
        origin = (outer_lat[0], outer_lon[0])
    x, y = convert_to_cartesian(outer_lat, outer_lon, origin[0], origin[1])
    if polygon_intersects_sweep(x, y):
        # self intersecting
        return None, None, None
    return x, y, origin

def wrap_180(diff):
    if diff > 180:
        wrap_180(diff - 360)
    elif diff < -180:
        wrap_180(diff + 360)
    return diff

def longitude_scale(lat):
    scale = math.cos(math.radians(lat))
    return max(scale, 0.01)

# convert lat lon to xy relative to origin point
LATLON_TO_M = 6378100 * (math.pi / 180)
def convert_to_cartesian(lat, lon, origin_lat, origin_lon):
    num_nodes = len(lat)
    x = np.empty(num_nodes)
    y = np.empty(num_nodes)
    for i in range(num_nodes):
        x[i] = (lat[i]-origin_lat) * LATLON_TO_M
        y[i] = wrap_180(lon[i] - origin_lon) * LATLON_TO_M * longitude_scale((lat[i]+origin_lat)*0.5)
    return x, y

# convert xy back to lat lon
def convert_from_cartesian(x, y, origin_lat, origin_lon):
    num_nodes = len(x)
    lat = np.empty(num_nodes)
    lon = np.empty(num_nodes)
    for i in range(num_nodes):
        dlat = x[i]/LATLON_TO_M
        lon[i] = wrap_180(origin_lon + ((y[i]/LATLON_TO_M) / longitude_scale(origin_lat+dlat/2)))
        lat[i] = origin_lat + dlat
    return lat, lon

# detect intersection between two points
def line_intersects(seg1_start, seg1_end, seg2_start, seg2_end):

    # do Y first, X will not trip during sweep line intersection in X axis
    min_y_1 = min(seg1_start[1], seg1_end[1])
    max_y_2 = max(seg2_start[1], seg2_end[1])

    if min_y_1 > max_y_2:
        return False

    max_y_1 = max(seg1_start[1], seg1_end[1])
    min_y_2 = min(seg2_start[1], seg2_end[1])

    if max_y_1 < min_y_2:
        return False

    min_x_1 = min(seg1_start[0], seg1_end[0])
    max_x_2 = max(seg2_start[0], seg2_end[0])

    if min_x_1 > max_x_2:
        return False

    max_x_1 = max(seg1_start[0], seg1_end[0])
    min_x_2 = min(seg2_start[0], seg2_end[0])

    if max_x_1 < min_x_2:
        return False

    # implementation borrowed from http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    r1 = (seg1_end[0] - seg1_start[0], seg1_end[1] - seg1_start[1])
    r2 = (seg2_end[0] - seg2_start[0], seg2_end[1] - seg2_start[1])
    r1xr2 = r1[0]*r2[1] - r1[1]*r2[0]
    if abs(r1xr2) < 1e-09:
        # either collinear or parallel and non-intersecting
        return False
    else:
        ss2_ss1 = (seg2_start[0] - seg1_start[0], seg2_start[1] - seg1_start[1])
        q_pxr = ss2_ss1[0]*r1[1] - ss2_ss1[1]*r1[0]
        # t = (q - p) * s / (r * s)
        # u = (q - p) * r / (r * s)
        t = (ss2_ss1[0]*r2[1] - ss2_ss1[1]*r2[0]) / r1xr2
        u = q_pxr / r1xr2
        if (u >= 0) and (u <= 1) and (t >= 0) and (t <= 1):
            # lines intersect
            # t can be any non-negative value because (p, p + r) is a ray
            # u must be between 0 and 1 because (q, q + s) is a line segment
            #intersection = seg1_start + (r1*t);
            return True
        else:
            # non-parallel and non-intersecting
            return False

# detect polygon self intersection
# https://github.com/rowanwins/sweepline-intersections
def polygon_intersects_sweep(x, y):
    # list of lines in polygon
    num_nodes = len(x)

    lines = {}
    event_que = 2*num_nodes*[None]
    for i in range(num_nodes):
        j = i+1
        if j >= num_nodes:
            j = 0
        lines[i] = ((x[i],y[i]), (x[j], y[j]))
        if x[i] <= x[j]:
            event_que[i*2] = (x[i],i,True)
            event_que[(i*2)+1] = (x[j],i,False)
        else:
            event_que[i*2] = (x[i],i,False)
            event_que[(i*2)+1] = (x[j],i,True)

    event_que = sorted(event_que, key=lambda event: event[0])
    active = {}
    for event in event_que:
        if event[2]:
            # adding new line, intersect with active items
            new_line = lines[event[1]]

            # don't compare adjacent lines in polygon
            next_line = event[1] + 1
            if next_line == num_nodes:
                next_line = 0
            prev_line = event[1] - 1
            if prev_line == -1:
                prev_line = num_nodes - 1

            for key, line in active.items():
                if (key == next_line) or (key == prev_line):
                    continue
                if line_intersects(new_line[0], new_line[1], line[0], line[1]):
                    return True
            active[event[1]] = new_line

        else:
            # remove line from active list
            active.pop(event[1])

    return False

# Check for intersections between multiple polygons useing sweep line
def polygon_polygon_intersection(x, y):
    num_poly = len(x)

    lines = {}
    event_que = []
    key = 0
    for k in range(num_poly):
        num_nodes = len(x[k])
        for i in range(num_nodes):
            j = i+1
            if j >= num_nodes:
                j = 0
            lines[key] = ((x[k][i], y[k][i]), (x[k][j], y[k][j]), k)
            right_to_left = x[k][i] <= x[k][j]
            event_que.append((x[k][i], key, right_to_left))
            event_que.append((x[k][j], key, not right_to_left))
            key += 1

    event_que = sorted(event_que, key=lambda event: event[0])
    active = {}
    for event in event_que:
        if event[2]:
            # adding new line, intersect with active items
            new_line = lines[event[1]]

            for line in active.values():
                if new_line[2] == line[2]:
                    # don't compare lines in same polygon
                    continue
                if line_intersects(new_line[0], new_line[1], line[0], line[1]):
                    return True
            active[event[1]] = new_line

        else:
            # remove line from active list
            active.pop(event[1])

    return False

# https://en.wikipedia.org/wiki/Shoelace_formula
def polygon_area(x, y):
    sum1 = x[-1] * y[0]
    sum2 = x[0] * y[-1]
    for i in range(len(x)-1):
        sum1 += x[i] * y[i+1]
        sum2 += y[i] * x[i+1]
    return abs(sum1 - sum2) * 0.5

# as above for thee points
def triangle_area(x, y):
    sum1 = x[2] * y[0] + x[0] * y[1] + x[1] * y[2]
    sum2 = x[0] * y[2] + y[0] * x[1] + y[1] * x[2]
    return abs(sum1 - sum2) * 0.5

# retrun true if point is outside of given polygon
def point_outside_polygon(point_x, point_y, poly_x, poly_y):
    # step through each edge pair-wise looking for crossings:
    num_nodes = len(poly_x)
    outside = True
    min_x = min(poly_x) - 1

    for i in range(num_nodes):
        j = i+1
        if j >= num_nodes:
            j = 0
        if (poly_y[i] > point_y) == (poly_y[j] > point_y):
            # both ends of line are on the same side of the point
            # no intersection possible
            continue
        if line_intersects((min_x, point_y), (point_x, point_y) , (poly_x[i], poly_y[i]), (poly_x[j], poly_y[j])):
            outside = not outside

    return outside

# simplify polygon using Visvalingam–Whyatt
# https://en.wikipedia.org/wiki/Visvalingam%E2%80%93Whyatt_algorithm
# will not create self intersecting polygon
def simplify_poly(x, y):
    poly_len = {}
    num_poly = len(x)
    minimum_polygon = [False] * num_poly
    circle_radius = [None] * num_poly
    for i in range(num_poly):
        poly_len[i] = len(x[i])
        if poly_len[i] <= 3:
            minimum_polygon[i] = True

    if all(minimum_polygon):
        # cant simplify any further
        return x, y, circle_radius

    if sum(poly_len.values()) <= simp_min_nodes:
        # already lower than min node threshold
        return x, y, circle_radius

    for i in range(num_poly):
        # try replacing polygon with circle
        center_x = np.mean(x[i])
        center_y = np.mean(y[i])
        radius = np.empty(poly_len[i])
        for j in range(poly_len[i]):
            radius[j] = math.sqrt(((x[i][j] - center_x) ** 2) + ((y[i][j] - center_y) ** 2))
        radius = np.mean(radius)
        circle_area = math.pi * (radius ** 2)
        poly_area = polygon_area(x[i], y[i])
        if abs(circle_area -  poly_area) < simp_area_threshold:
            x[i] = [center_x]
            y[i] = [center_y]
            circle_radius[i] = radius
            minimum_polygon[i] = True
            poly_len[i] = 1

    if all(minimum_polygon):
        # cant simplify any further
        return x, y, circle_radius

    area = [None] * num_poly
    for i in range(num_poly):
        area[i] = [None] * poly_len[i]
        for j in range(poly_len[i]):
            if minimum_polygon[i]:
                continue
            prev_point = j - 1
            if prev_point < 0:
                prev_point = poly_len[i] - 1

            next_point = j + 1
            if next_point >= poly_len[i]:
                next_point = 0

            area[i][j] = triangle_area([x[i][j], x[i][prev_point], x[i][next_point]], [y[i][j], y[i][prev_point], y[i][next_point]])

    while True:
        min_poly_val = math.inf
        min_poly_index = None
        index_min = None
        for i in range(num_poly):
            if minimum_polygon[i]:
                continue
            temp_index_min = min(range(len(area[i])), key=area[i].__getitem__)
            if area[i][temp_index_min] < min_poly_val:
                min_poly_index = i
                index_min = temp_index_min
                min_poly_val = area[i][index_min]

        if min_poly_val > area_threshold and ((simp_max_nodes == None) or (sum(poly_len.values()) <= simp_max_nodes)):
            # reached threshold, simplification complete
            break

        # test if removing this point will create a self intersections
        new_intersect = False
        prev_point = index_min - 1
        if prev_point < 0:
            prev_point = poly_len[min_poly_index] - 1

        prev_prev_point = prev_point - 1
        if prev_prev_point < 0:
            prev_prev_point = poly_len[min_poly_index] - 1

        next_point = index_min + 1
        if next_point >= poly_len[min_poly_index]:
            next_point = 0

        for i in range(poly_len[min_poly_index]):
            # compare all lines except the adjacent
            if i == prev_prev_point or i == prev_point or i == index_min or i == next_point:
                continue
            test_next_point = i + 1
            if test_next_point >= poly_len[min_poly_index]:
                test_next_point = 0
            if line_intersects((x[min_poly_index][i], y[min_poly_index][i]), (x[min_poly_index][test_next_point], y[min_poly_index][test_next_point]), (x[min_poly_index][prev_point], y[min_poly_index][prev_point]), (x[min_poly_index][next_point], y[min_poly_index][next_point])):
                new_intersect = True
                break

        if new_intersect:
            # cant remove this point without creating intersection
            area[min_poly_index][index_min] = math.inf
            continue

        x[min_poly_index] = np.concatenate((x[min_poly_index][:index_min], x[min_poly_index][index_min+1:]))
        y[min_poly_index] = np.concatenate((y[min_poly_index][:index_min], y[min_poly_index][index_min+1:]))
        del area[min_poly_index][index_min]

        poly_len[min_poly_index] -= 1
        if poly_len[min_poly_index] == 3:
            # cant simplify past 3 points
            minimum_polygon[min_poly_index] = True

        if all(minimum_polygon):
            # cant simplify any further
            break

        if sum(poly_len.values()) <= simp_min_nodes:
            # reached min node threshold
            break

        # test if inner polygons are outside outer due to simplification and remove
        # detect intersection between polygons and merge

        # recalculate area for adjacent points
        for j in [index_min-1, index_min]:

            if j >= poly_len[min_poly_index]:
                j = 0
            elif j < 0:
                j = poly_len[min_poly_index] - 1

            prev_point = j - 1
            if prev_point < 0:
                prev_point = poly_len[min_poly_index] - 1

            next_point = j + 1
            if next_point >= poly_len[min_poly_index]:
                next_point = 0

            area[min_poly_index][j] = triangle_area([x[min_poly_index][j], x[min_poly_index][prev_point], x[min_poly_index][next_point]], [y[min_poly_index][j], y[min_poly_index][prev_point], y[min_poly_index][next_point]])

        # recalculate any areas set to inf to avoid intersections
        for j in range(poly_len[min_poly_index]):
            if math.isinf(area[min_poly_index][j]):
                prev_point = j - 1
                if prev_point < 0:
                    prev_point = poly_len[min_poly_index] - 1
                next_point = j + 1
                if next_point >= poly_len[min_poly_index]:
                    next_point = 0
                area[min_poly_index][j] = triangle_area([x[min_poly_index][j], x[min_poly_index][prev_point], x[min_poly_index][next_point]], [y[min_poly_index][j], y[min_poly_index][prev_point], y[min_poly_index][next_point]])

    return x, y, circle_radius

def save_to_file(name, tag, center, lat, lon, radius):
    # sanitize name for use in file
    name = name.replace('/', '_')
    name = name.replace('\\', '_')

    file_name = '%s-%s:%f:%f.waypoints' % (name, tag, center[0], center[1])
    f = open(os.path.join(directory, file_name), "w")
    f.write('QGC WPL 110\n')
    total_points = 1
    for i in range(len(lat)):
        if i == 0:
            # first point is always inclusion
            poly_type = 5001
            circle_type = 5003
        else:
            # all others exclusion
            poly_type = 5002
            circle_type = 5004
        
        if radius[i] is None:
            # polygon points
            length = len(lat[i])
            for j in range(length):
                f.write('%i 0 3 %i %i 0 0 0 %f %f %i 1\n' % (total_points, poly_type, length, lat[i][j], lon[i][j], j))
                total_points += 1
        else:
            # Circle point
            f.write('%i 0 3 %i %f 0 0 0 %f %f 0 1\n' % (total_points, circle_type, radius[i], lat[i][0], lon[i][0]))
            total_points += 1
    f.close()
    return file_name


start_time = time.time()

# enable or disable profiling for speed
profiling = False
if profiling:
    import cProfile, pstats, io
    pr = cProfile.Profile()
    pr.enable()

# delete existing fences in directory
files = os.listdir(directory)
for _file in files:
    if _file.endswith(".waypoints"):
        os.remove(os.path.join(directory, _file))

# export a .js file containing outer points of all polygons
js_file = open(os.path.join(directory,'data.js'), "w") 
js_file.write('var fence_data = [\n')

# search input file
fence_search().apply_file(input_file)

# close js file
js_file.write(']\n')
js_file.close()

print("Took %0.2fs" % (time.time() - start_time))

if profiling:
    pr.disable()
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())
