import gmsh
import sys

from numpy.f2py.capi_maps import lcb2_map
from numpy.polynomial.polynomial import polygrid3d

gmsh.initialize(sys.argv)

gmsh.model.add("mixed")

lc = 0.25

p1  = gmsh.model.geo.addPoint(0.15, 0.15, 0.0, lc)
p2  = gmsh.model.geo.addPoint(-0.15, 0.15, 0.0, lc)
p3 = gmsh.model.geo.addPoint(-0.15, -0.15, 0.0, lc)
p4 = gmsh.model.geo.addPoint(0.15, -0.15, 0.0, lc)

# outer Square
p5 = gmsh.model.geo.addPoint(0.5, 0.5, 0.0, lc)
p6 = gmsh.model.geo.addPoint(-0.5, 0.5, 0.0, lc)
p7 = gmsh.model.geo.addPoint(-0.5, -0.5, 0.0, lc)
p8 = gmsh.model.geo.addPoint(0.5, -0.5, 0.0, lc)


# inner Square
l1 = gmsh.model.geo.addLine(p4, p1)
l2 = gmsh.model.geo.addLine(p1, p2)
l3 = gmsh.model.geo.addLine(p2, p3)
l4 = gmsh.model.geo.addLine(p3, p4)

# outer square
l5 = gmsh.model.geo.addLine(p8, p5)
l6 = gmsh.model.geo.addLine(p5, p6)
l7 = gmsh.model.geo.addLine(p6, p7)
l8 = gmsh.model.geo.addLine(p7, p8)

# surface
cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])

s1  = gmsh.model.geo.addPlaneSurface([cl1,cl2])

gmsh.model.geo.mesh.setTransfiniteCurve(l1, 7)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, 7)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 7)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, 7)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 21)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 21)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 21)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, 21)


gmsh.model.geo.mesh.setRecombine(2, s1)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l5, l6, l7, l8], 11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4],       12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1],13 , name="fluid")

# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("grado2.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()