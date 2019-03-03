import bpy
import os
from mathutils import *
D = bpy.data
C = bpy.context
O = C.object

class Material:
    def __init__(self, color):
        self.color = color
    def emit(self, f):
        s = ("texture "
             "reflect(%f),diffuse(%f),specular(%f),color(%f,%f,%f),"
             "refract(%f),blur(%f)\n")
        (R,G,B) = self.color[:3]
        # TODO: fill in the missing inputs!
        out = s % (0,0,0,R,G,B,0,0)
        f.write(out)

class Face:
    def __init__(self, coords):
        assert len(coords) == 9 or len(coords) == 12, "tris or quads!"
        self.coords = coords
    def emit(self, f):
        f.write(' '.join([str(x) for x in self.coords]))
        f.write('\n')

class Mesh:
    def __init__(self, name):
        self.faces = []
        self.name = name
    def addFace(self, face):
        self.faces.append(face)
    def emit(self, f):
        f.write('file %s\n' % (self.name))
        with open(self.name+'.txt', 'w') as meshFile:
            for face in self.faces:
                face.emit(meshFile)

def getMesh(obj):
    faces = obj.data.polygons
    mesh = Mesh(obj.name)
    for face in faces:
        # these are indices into the obj vertices
        faceCoords = []
        for vertIdx in face.vertices:
            vertex = obj.data.vertices[vertIdx]
            vertexWorld = obj.matrix_world @ vertex.co
            (x, y, z) = vertexWorld
            faceCoords += [x, y, z]
        mesh.addFace(Face(faceCoords))

    return mesh

def getMaterialProperties(obj):
    # TODO: for when we support more materials!
    #for face in obj.polygons:
    #    matSlot = obj.material_slots[face.material_index]
    #    mat = matSlot.material
    #    color = mat.diffuse_color
    slot = obj.material_slots[0]
    mat = slot.material
    node = mat.node_tree.nodes['Principled BSDF']
    inputs = node.inputs
    # 4-tuple (RGBA normalized values)
    color = tuple(inputs['Base Color'].default_value)
    return Material(color=color)

def emitObject(f, obj):
    if obj.type == 'MESH':
        assert len(obj.material_slots) == 1, "Only 1 material slot for now!"
        mesh = getMesh(obj)
        mat  = getMaterialProperties(obj)
        mesh.emit(f)
        mat.emit(f)
    elif obj.type == 'LIGHT':
        pass
    elif obj.type == 'CAMERA':
        pass
    else:
        assert False, "unhandled type!"

# Just change to the blend file directory to dump the files
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
dname = os.path.dirname(dname)
os.chdir(dname)

allObjects = D.objects

with open('scene_data.txt', 'w') as f:
    for obj in allObjects:
        emitObject(f, obj)

print("Finished!")