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

#Midpoints
p9 = gmsh.model.geo.addPoint(0.15, 0.5, 0.0, lc)
p10 = gmsh.model.geo.addPoint(-0.15, 0.5, 0.0, lc)
p11 = gmsh.model.geo.addPoint(-0.5, 0.15, 0.0, lc)
p12 = gmsh.model.geo.addPoint(-0.5, -0.15, 0.0, lc)

p13 = gmsh.model.geo.addPoint(-0.15, -0.5,0,lc)
p14 = gmsh.model.geo.addPoint(0.15, -0.5,0,lc)
p15 = gmsh.model.geo.addPoint(0.5, -0.15,0,lc)
p16 = gmsh.model.geo.addPoint(0.5, 0.15,0,lc)



# Region 1
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
# Region 2
l5 = gmsh.model.geo.addLine(p1, p9)
l6 = gmsh.model.geo.addLine(p9, p5)
l7 = gmsh.model.geo.addLine(p5, p16)
l8 = gmsh.model.geo.addLine(p16, p1)
# Region 3
l9 = gmsh.model.geo.addLine(p9, p10)
l10 = gmsh.model.geo.addLine(p10, p2)
# Region 4
l11 = gmsh.model.geo.addLine(p10, p6)
l12 = gmsh.model.geo.addLine(p6, p11)
l13 = gmsh.model.geo.addLine(p11, p2)
# Region 5
l14 = gmsh.model.geo.addLine(p11, p12)
l15 = gmsh.model.geo.addLine(p12, p3)
# Region 6
l16 = gmsh.model.geo.addLine(p12, p7)
l17 = gmsh.model.geo.addLine(p7, p13)
l18 = gmsh.model.geo.addLine(p13, p3)
# Region 7
l19 = gmsh.model.geo.addLine(p13, p14)
l20 = gmsh.model.geo.addLine(p14, p4)
# Region 8
l21 = gmsh.model.geo.addLine(p4, p15)
l22 = gmsh.model.geo.addLine(p15, p8)
l23 = gmsh.model.geo.addLine(p8, p14)
# Region 9
l24 = gmsh.model.geo.addLine(p15, p16)

# surface

cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
'''
s1  = gmsh.model.geo.addPlaneSurface([cl1])
'''
cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
s2 = gmsh.model.geo.addPlaneSurface([cl2])

cl3 = gmsh.model.geo.addCurveLoop([l5, l9, l10, -l1])
s3 = gmsh.model.geo.addPlaneSurface([cl3])

cl4 = gmsh.model.geo.addCurveLoop([-l10, l11, l12, l13])
s4 = gmsh.model.geo.addPlaneSurface([cl4])

cl5 = gmsh.model.geo.addCurveLoop([-l13, l14, l15, -l2])
s5 = gmsh.model.geo.addPlaneSurface([cl5])

cl6 = gmsh.model.geo.addCurveLoop([-l15, l16, l17, l18])
s6 = gmsh.model.geo.addPlaneSurface([cl6])

cl7 = gmsh.model.geo.addCurveLoop([-l18, l19, l20, -l3])
s7 = gmsh.model.geo.addPlaneSurface([cl7])

cl8 = gmsh.model.geo.addCurveLoop([l21, l22, l23, l20])
s8 = gmsh.model.geo.addPlaneSurface([cl8])

cl9 = gmsh.model.geo.addCurveLoop([l21, l24, l8, -l4])
s9 = gmsh.model.geo.addPlaneSurface([cl9])


gmsh.model.geo.mesh.setTransfiniteCurve(l1, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l8, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l9, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l11, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l12, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l13, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l15, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l16, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l17, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l18, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l19, 9)
gmsh.model.geo.mesh.setTransfiniteCurve(l20, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l21, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l22, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l23, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l24, 9)

gmsh.model.geo.mesh.setTransfiniteSurface(s2)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setTransfiniteSurface(s4)
gmsh.model.geo.mesh.setTransfiniteSurface(s5)
gmsh.model.geo.mesh.setTransfiniteSurface(s6)
gmsh.model.geo.mesh.setTransfiniteSurface(s7)
gmsh.model.geo.mesh.setTransfiniteSurface(s8)
gmsh.model.geo.mesh.setTransfiniteSurface(s9)


gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)
gmsh.model.geo.mesh.setRecombine(2, s5)
gmsh.model.geo.mesh.setRecombine(2, s6)
gmsh.model.geo.mesh.setRecombine(2, s7)
gmsh.model.geo.mesh.setRecombine(2, s8)
gmsh.model.geo.mesh.setRecombine(2, s9)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l6, l9, l11, l12, l14, l16, l17, l19 ,l23, l22, l24, l7], 11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1,l2,l3,l4],       12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s2, s3, s4, s5, s6, s7, s8, s9],13 , name="fluid")

# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("gradok2.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()