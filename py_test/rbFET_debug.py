from netgen.occ import *
import sys
sys.path.append("../build_debug/")
from rigid_body_FEM import *
from rigid_body_FEM.bla import *
import pythreejs as p3
from IPython.display import display
import ipywidgets as widgets


def extract_vertices(obj: TopoDS_Shape):
    "extracts a p3js compatible vertex list from a netgen.occ TopoDS_Shape"
    
    data = obj._webgui_data()["Bezier_trig_points"]
    
    # for every face, each of the verts arrays holds one vertex
    verts1 = data[0]
    verts2 = data[1]
    verts3 = data[2]
    
    # corresponding normals
    normals1 = data[3]
    normals2 = data[4]
    normals3 = data[5]

    combined_vertices = []
    for i in range(0, len(verts1), 4):
        combined_vertices.append(verts1[i : i+3])
        combined_vertices.append(verts2[i : i+3])
        combined_vertices.append(verts3[i : i+3])

    combined_normals = []
    for i in range(0, len(normals1), 3):
        combined_normals.append(normals1[i : i+3])
        combined_normals.append(normals2[i : i+3])
        combined_normals.append(normals3[i : i+3])
            
    return combined_vertices, combined_normals


def body_from_solid(obj):
    "extracts the mass matrix of the TopAbs_ShapeEnum.SOLID obj, using the figures computed by netgen"

    # important: move the center of mass into the origin
    obj = obj.Move((-obj.center[0], -obj.center[1], -obj.center[2]))
    
    # copy the inertia matrix from netgen
    inertia_matrix = bla.Matrix(3,3)
    for i in range(3):
       for j in range(3):
           inertia_matrix[i, j] = obj.inertia[i, j]

    # inertia_matrix[0,0] = 1
    # inertia_matrix[0,1] = 0
    # inertia_matrix[0,2] = 0
    # inertia_matrix[1,0] = 0
    # inertia_matrix[1,1] = 1
    # inertia_matrix[1,2] = 0
    # inertia_matrix[2,0] = 0
    # inertia_matrix[2,1] = 0
    # inertia_matrix[2,2] = 1

    # copy the center of mass from netgen
    center_of_mass = bla.Vector(3)
    for i in range(3): center_of_mass[i] = obj.center[i]

    # rearrange it in C++ to make the mass matrix (the elegant way, using MatrixView)
    body = RigidBody_FEM()
    #for i in range(3) : body.center[i] = obj.center[i]
    #body.mass = obj.mass
    #body.inertia = inertia_matrix
    #body.recalcMassMatrix()

    body.vertices, body.normals = extract_vertices(obj)
    
    return body

def appendConnector(rbs, c,connectors):
    p = rbs.connectorPos(c)
    if(c.type == 0):
        color = 'green'
    else :
        color = 'black'
    connectors.append(
        p3.Mesh(p3.SphereBufferGeometry(0.2, 16, 16),
             p3.MeshStandardMaterial(color=color),
             position=(p[0], p[1], p[2])))

def initConnectors(rbs, l):
    connectors = []
    for s in l:
        cA = s.connectorA
        cB = s.connectorB
        appendConnector(rbs, cA,connectors);
        appendConnector(rbs, cB,connectors);
    return connectors
    

def updateConnectors(rbs, l, connectors):
    for i in range(len(l)):
        cA = l[i].connectorA
        cB = l[i].connectorB
        pA = rbs.connectorPos(cA);
        pB = rbs.connectorPos(cB);
        connectors[2*i].position = (pA[0], pA[1], pA[2])
        connectors[2*i+1].position = (pB[0], pB[1], pB[2])

def positionsOf(rbs, l):
    res = []
    for s in l:
        cA = s.connectorA
        cB = s.connectorB
        pA = rbs.connectorPos(cA);
        pB = rbs.connectorPos(cB);
        res.append ([ [pA[0], pA[1], pA[2]], [pB[0], pB[1], pB[2]] ])
    return res


#connectorsSprings = initConnectors(rbs.springs())
#springpos = positionsOf(rbs.springs())

#if springpos:
#    springgeo = p3.LineSegmentsGeometry(positions=springpos)
#    m2 = p3.LineMaterial(linewidth=3, color='cyan')
#    springs = p3.LineSegments2(springgeo, m2)

#connectorsBeams = initConnectors(rbs.beams())
#beampos = positionsOf(rbs.beams())

#if beampos:
#    beamgeo = p3.LineSegmentsGeometry(positions=beampos)
#    m2 = p3.LineMaterial(linewidth=4, color='blue')
#    beams = p3.LineSegments2(beamgeo, m2)

def display_system(rbs, step_size):
    connectorsBeams = initConnectors(rbs, rbs.beams())
    beampos = positionsOf(rbs, rbs.beams())
    if beampos:
      beamgeo = p3.LineSegmentsGeometry(positions=beampos)
      m2 = p3.LineMaterial(linewidth=4, color='blue')
      beams = p3.LineSegments2(beamgeo, m2)

    view_width = 1000
    view_height = 700
    buffergeos = []
    p3meshes = []

    # set up pythreejs 3d objects
    for body in rbs.bodies():
        buffergeom = p3.BufferGeometry(attributes = {"position" : p3.BufferAttribute(body.vertices), "normal" : p3.BufferAttribute(body.normals)})
        material = p3.MeshPhongMaterial(color='#ff3333', shininess=150, morphTargets=True, side="DoubleSide")
        p3mesh = p3.Mesh(buffergeom, material, position=(0,0,0))
        buffergeos.append(buffergeom)
        p3meshes.append(p3mesh)

    # extra scene contents
    camera = p3.PerspectiveCamera( position=[10, 6, 10], aspect=view_width/view_height)
    key_light = p3.DirectionalLight(position=[0, 10, 10])
    ambient_light = p3.AmbientLight()
    grid = p3.GridHelper(500, 500//5, "#2F4F4F", "#2F4F4F")
    axesHelper = p3.AxesHelper(5)

    # set up scene
    scene = p3.Scene(children=[camera, key_light, ambient_light, grid, axesHelper, *p3meshes, *connectorsBeams] + ([] if not beampos else [beams]))
    controller = p3.OrbitControls(controlling=camera)
    renderer = p3.Renderer(camera=camera, scene=scene, controls=[controller],
                        width=view_width, height=view_height, antialias=True) # if performance is bad, try removing antalias=True

    # play/reset button
    play = widgets.Play(
        value=0,
        min=0,
        max=100000,
        step=1,
        interval=1,
        description="Press play",
        disabled=False
    )

    for m in p3meshes:
        m.matrixAutoUpdate = False # make mesh movable
    
    eq = rbs.assemble(step_size)

    def refresh():
        "updates all pythreejs object transformations"
        for i in range(0, len(rbs.bodies())):
            p3meshes[i].matrix=rbs.bodies()[i].asTuple()
        #updateConnectors(rbs.springs(),connectorsSprings)
        #springpos = positionsOf(rbs.springs())
        #if springpos:
        #    springs.geometry = p3.LineSegmentsGeometry(positions=springpos)
        
        updateConnectors(rbs, rbs.beams(),connectorsBeams)
        beampos = positionsOf(rbs, rbs.beams())
        if beampos:
            beams.geometry = p3.LineSegmentsGeometry(positions=beampos)

    def update():
        "update function, gets called every timestep; quasi main loop"
        dt = 0.0005
        eq.step()
        #display_energy(rbs, dt)
        refresh()

    def observer(state):
        "event handler for clickable buttons"
        # if there is a change in time
        if state["name"] == "value":
            # it might be a reset
            if str(state["new"]) == "0":
                eq.reset()
                refresh()
            # or it might be a progress in time
            else:
                update()
        # repeat is used as an alias to reset
        elif state["name"] == "repeat":
            eq.reset()
            refresh()

    play.observe(observer)

    display(widgets.HBox([play, widgets.HTML("<b>click-and-drag to rotate, scroll to zoom, right-click-and-drag to move<b>")]))
    display(renderer)
    refresh()
