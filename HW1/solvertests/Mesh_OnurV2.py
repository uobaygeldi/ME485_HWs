import gmsh
import sys

from numpy.f2py.capi_maps import lcb2_map

gmsh.initialize(sys.argv)

gmsh.model.add("mixed")

lc = 0.25

# Center point
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)

# Inner boundary
p2 = gmsh.model.geo.addPoint(0.1, -0.1, 0, lc)
p3 = gmsh.model.geo.addPoint(0.1, 0.1, 0, lc)
p4 = gmsh.model.geo.addPoint(-0.1, 0.1, 0, lc)
p5 = gmsh.model.geo.addPoint(-0.1, -0.1, 0, lc)

# Outer boundary
p10 = gmsh.model.geo.addPoint(1, -1, 0, lc)
p11 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p12 = gmsh.model.geo.addPoint(-1, 1, 0, lc)
p13 = gmsh.model.geo.addPoint(-1, -1, 0, lc)

# square
p6 = gmsh.model.geo.addPoint(0.5, -0.5, 0, lc)
p7 = gmsh.model.geo.addPoint(0.5, 0.5, 0, lc)
p8 = gmsh.model.geo.addPoint(-0.5, 0.5, 0, lc)
p9 = gmsh.model.geo.addPoint(-0.5, -0.5, 0, lc)

# Outer circle
l1 = gmsh.model.geo.addCircleArc(p10, p1, p11)
l2 = gmsh.model.geo.addCircleArc(p11, p1, p12)
l3 = gmsh.model.geo.addCircleArc(p12, p1, p13)
l4 = gmsh.model.geo.addCircleArc(p13, p1, p10)

# Inner circle
l5 = gmsh.model.geo.addCircleArc(p2, p1, p3)
l6 = gmsh.model.geo.addCircleArc(p3, p1, p4)
l7 = gmsh.model.geo.addCircleArc(p4, p1, p5)
l8 = gmsh.model.geo.addCircleArc(p5, p1, p2)

#spuare
l9 = gmsh.model.geo.addLine(p6, p7)
l10 = gmsh.model.geo.addLine(p7, p8)
l11 = gmsh.model.geo.addLine(p8, p9)
l12 = gmsh.model.geo.addLine(p9, p6)

#square to outer circle
l13 = gmsh.model.geo.addLine(p6, p10)
l14 = gmsh.model.geo.addLine(p7, p11)
l15 = gmsh.model.geo.addLine(p8, p12)
l16 = gmsh.model.geo.addLine(p9, p13)

#surface definition
cl1 = gmsh.model.geo.addCurveLoop([l1, -l14, -l9, l13])
cl2 = gmsh.model.geo.addCurveLoop([l2, -l15, -l10, l14])
cl3 = gmsh.model.geo.addCurveLoop([l3, -l16, -l11, l15])
cl4 = gmsh.model.geo.addCurveLoop([l4, -l13, -l12, l16])
cl5 = gmsh.model.geo.addCurveLoop([l9, l10, l11, l12])
cl6 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])

s5 = gmsh.model.geo.addPlaneSurface([cl5,cl6]) #surface 5
s1 = gmsh.model.geo.addPlaneSurface([cl1]) #surface 1
s2 = gmsh.model.geo.addPlaneSurface([cl2]) #surface 2
s3 = gmsh.model.geo.addPlaneSurface([cl3]) #surface 3
s4 = gmsh.model.geo.addPlaneSurface([cl4]) #surface 4

gmsh.model.geo.mesh.setTransfiniteCurve(l1, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l9, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l11, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l12, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l13, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l15, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l16, 10)


gmsh.model.geo.mesh.setTransfiniteSurface(s1)
gmsh.model.geo.mesh.setTransfiniteSurface(s2)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setTransfiniteSurface(s4)

gmsh.model.geo.mesh.setRecombine(2, s1)
gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)


gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l9, l10, l11, l12], 11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1,l2,l3,l4],       12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5],13 , name="fluid")

# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("gradok2.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()


















