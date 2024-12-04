import gmsh
import sys

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize(sys.argv)

# Next add a new model named "cavity"
gmsh.model.add("mixed")
lc = 0.25

#Points
p1  = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc);

p2  = gmsh.model.geo.addPoint(-0.138,-0.19, 0,lc);
p3  = gmsh.model.geo.addPoint( 0.138,-0.19, 0,lc);
p4  = gmsh.model.geo.addPoint( 0.223, 0.073, 0, lc);
p5  = gmsh.model.geo.addPoint(-0.223, 0.073, 0,lc);
p6  = gmsh.model.geo.addPoint(0, 0.235, 0,lc);

p7  = gmsh.model.geo.addPoint(0,-0.5, 0,lc);
p8  = gmsh.model.geo.addPoint( 0.476,-0.155, 0,lc);
p9  = gmsh.model.geo.addPoint( 0.294, 0.405, 0,lc);
p10  = gmsh.model.geo.addPoint(-0.294, 0.405, 0,lc);
p11  = gmsh.model.geo.addPoint(-0.476, -0.155, 0,lc);

p12 = gmsh.model.geo.addPoint(0, -1, 0, lc);
p13 = gmsh.model.geo.addPoint(1, 0, 0, lc);
p14 = gmsh.model.geo.addPoint(0, 1, 0, lc);
p15 = gmsh.model.geo.addPoint(-1, 0, 0, lc);

p16 = gmsh.model.geo.addPoint(0.6, 0, 0, lc);
p17 = gmsh.model.geo.addPoint(0, 0.6, 0, lc);
p18 = gmsh.model.geo.addPoint(-0.6, 0, 0, lc);
p19 = gmsh.model.geo.addPoint(0, -0.6, 0, lc);

#As a general rule, elementary entity tags in Gmsh have to be unique per geometrical dimension.
# Star
l1   = gmsh.model.geo.addLine(p7  , p2)
l2   = gmsh.model.geo.addLine(p2  , p11)
l3   = gmsh.model.geo.addLine(p11 , p5)
l4   = gmsh.model.geo.addLine(p5  , p10)
l5   = gmsh.model.geo.addLine(p10 , p6)
l6   = gmsh.model.geo.addLine(p6  , p9)
l7   = gmsh.model.geo.addLine(p9  , p4)
l8   = gmsh.model.geo.addLine(p4  , p8)
l9   = gmsh.model.geo.addLine(p8  , p3)
l10  = gmsh.model.geo.addLine(p3  , p7)

'''# Circling the Star
l11  = gmsh.model.geo.addCircleArc(p7 , p1, p11)
l12  = gmsh.model.geo.addCircleArc(p11, p1, p10)
l13  = gmsh.model.geo.addCircleArc(p10, p1, p9)
l14  = gmsh.model.geo.addCircleArc(p9 , p1, p8)
l15  = gmsh.model.geo.addCircleArc(p8 , p1, p7)'''

# Outer boundaries
l16  = gmsh.model.geo.addCircleArc(p12, p1, p15)
l17  = gmsh.model.geo.addCircleArc(p15, p1, p14)
l18  = gmsh.model.geo.addCircleArc(p14, p1, p13)
l19  = gmsh.model.geo.addCircleArc(p13, p1, p12)

l20 = gmsh.model.geo.addCircleArc(p19 , p1, p18)
l21 = gmsh.model.geo.addCircleArc(p18, p1, p17)
l22 = gmsh.model.geo.addCircleArc(p17, p1, p16)
l23 = gmsh.model.geo.addCircleArc(p16, p1, p19)

l24 = gmsh.model.geo.addLine(p12  , p19)
l25 = gmsh.model.geo.addLine(p15  , p18)
l26 = gmsh.model.geo.addLine(p14  , p17)
l27 = gmsh.model.geo.addLine(p13  , p16)

#Surfaces
cl1 = gmsh.model.geo.addCurveLoop([l24, l20, -l25, -l16])
s1  = gmsh.model.geo.addPlaneSurface([cl1])

cl2 = gmsh.model.geo.addCurveLoop([l25, l21, -l26, -l17])
s2  = gmsh.model.geo.addPlaneSurface([cl2])

cl3 = gmsh.model.geo.addCurveLoop([l26, l22, -l27, -l18])
s3  = gmsh.model.geo.addPlaneSurface([cl3])

cl4 = gmsh.model.geo.addCurveLoop([l27, l23, -l24, -l19])
s4  = gmsh.model.geo.addPlaneSurface([cl4])

'''cl5 = gmsh.model.geo.addCurveLoop([l1, l2, -l11])
cl6 = gmsh.model.geo.addCurveLoop([l3, l4, -l12])
cl7 = gmsh.model.geo.addCurveLoop([l5, l6, -l13])
cl8 = gmsh.model.geo.addCurveLoop([l7, l8, -l14])
cl9 = gmsh.model.geo.addCurveLoop([l9, l10, -l15])
s5  = gmsh.model.geo.addPlaneSurface([cl5,cl6,cl7,cl8,cl9])'''

cl10 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7, l8, l9, l10])
cl11 = gmsh.model.geo.addCurveLoop([l20, l21, l22, l23])
s6   = gmsh.model.geo.addPlaneSurface([cl10, cl11])

# # The `setTransfiniteCurve()' meshing constraints explicitly specifies the
# # location of the nodes on the curve. For example, the following command forces
# # 10 uniformly placed nodes on curve 2 (including the nodes on the two end
# # points):
gmsh.model.geo.mesh.setTransfiniteCurve(l16, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l17, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l18, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l19, 10)

gmsh.model.geo.mesh.setTransfiniteCurve(l20, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l21, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l22, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l23, 10)

gmsh.model.geo.mesh.setTransfiniteCurve(l24, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l25, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l26, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l27, 10)


gmsh.model.geo.mesh.setTransfiniteCurve(l1 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l2 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l3 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l4 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l5 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l6 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l7 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l8 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l9 , 5)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 5)

gmsh.model.geo.mesh.setTransfiniteSurface(s1)
gmsh.model.geo.mesh.setTransfiniteSurface(s2)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setTransfiniteSurface(s4)

gmsh.model.geo.mesh.setRecombine(2, s1)
gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l16, l17, l18, l19],11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10],12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s6],13 , name="fluid")

# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("grad_star.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()