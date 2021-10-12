"""
Filter all objects with a given tag, to speed up debugging of problem area in big OSM file
"""

import osmium as o
import os

class search(o.SimpleHandler):

    def __init__(self, tag_name, tag_value):
        super(search, self).__init__()
        self.ways = set()
        self.relations = set()
        self.tag_name = tag_name
        self.tag_value = tag_value

    def way(self, w):
        if (self.tag_name in w.tags) and (self.tag_value in w.tags[self.tag_name]):
            self.ways.add(w.id)

    def relation(self, r):
        if ('name' in r.tags) and (self.tag_value in r.tags[self.tag_name]):
            self.relations.add(r.id)
            for m in r.members:
                if m.type == 'w':
                   self.ways.add(m.ref)

class way_search(o.SimpleHandler):

    def __init__(self, ways):
        super(way_search, self).__init__()
        self.ways = ways
        self.nodes = set()

    def way(self, w):
        if w.id in self.ways:
            for n in w.nodes:
                self.nodes.add(n.ref)

class FeatureWriter(o.SimpleHandler):

    def __init__(self, writer, nodes, ways, relations):
        super(FeatureWriter, self).__init__()
        self.writer = writer
        self.nodes = nodes
        self.ways = ways
        self.relations = relations

    def node(self, n):
        if n.id in self.nodes:
            self.writer.add_node(n)

    def way(self, w):
        if w.id in self.ways:
            self.writer.add_way(w)

    def relation(self, r):
        if r.id in self.relations:
            self.writer.add_relation(r)


input_file = 'switzerland-padded.osm.pbf'

tag_name = 'name'
tag_value = 'Bodensee'

# Find Ways and Relations with correct name
features = search(tag_name, tag_value)
features.apply_file(input_file)

if (len(features.ways) + len(features.relations)) == 0:
    print('Could no find name %s in %s' % (name, input_file))
    exit()

output_name = input_file.split('.')[0] + ' - ' + tag_name + ' = ' + tag_value + '.osm'
try:
    os.remove(output_name)
except:
    pass

# get matching all nodes
way_searcher = way_search(features.ways)
way_searcher.apply_file(input_file)

# go through the file again and write out the marked data
writer = o.SimpleWriter(output_name)
FeatureWriter(writer, way_searcher.nodes, features.ways, features.relations).apply_file(input_file)
writer.close()

print('saved to: %s' % (output_name))