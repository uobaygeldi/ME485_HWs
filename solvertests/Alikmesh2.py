import gmsh
import sys
import math

#from share.doc.gmsh.examples.api.naca_boundary_layer_3d import l_2_12

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize(sys.argv)

# Next add a new model named "cavity"
gmsh.model.add("mixed")
lc = 0.25

# ALIK MESH
p1  = gmsh.model.geo.addPoint(0,0,0, lc)
p2  = gmsh.model.geo.addPoint(3.0, 0.0, 0.0, lc);
p3  = gmsh.model.geo.addPoint(3+0.3*0.6, -0.3*0.8, 0,lc);
p4  = gmsh.model.geo.addPoint(3+0.3*0.6,  0.3*0.8, 0,lc);
p5  = gmsh.model.geo.addPoint(-3.0, 0.0, 0.0, lc);
p6  = gmsh.model.geo.addPoint(-3-0.3*0.6, -0.3*0.8, 0,lc);
p7  = gmsh.model.geo.addPoint(-3-0.3*0.6,  0.3*0.8, 0,lc);
p8  = gmsh.model.geo.addPoint(4*math.sqrt(2)/2, 4*math.sqrt(2)/2, 0,lc);
p9  = gmsh.model.geo.addPoint(0,4,0);
p10 = gmsh.model.geo.addPoint(-4*math.sqrt(2)/2, 4*math.sqrt(2)/2, 0,lc);
p11 = gmsh.model.geo.addPoint(-4*math.sqrt(2)/2, -4*math.sqrt(2)/2, 0,lc);
p12 = gmsh.model.geo.addPoint(0,-4,0);
p13 = gmsh.model.geo.addPoint(4*math.sqrt(2)/2, -4*math.sqrt(2)/2, 0,lc);

p14 = gmsh.model.geo.addPoint(5,5,0)
p15 = gmsh.model.geo.addPoint(-5,5,0)
p16 = gmsh.model.geo.addPoint(-5,-5,0)
p17 = gmsh.model.geo.addPoint(5,-5,0)


# Curves
#Inner eye
l1 = gmsh.model.geo.addCircleArc(p3, p2, p4);
l2 = gmsh.model.geo.addCircleArc(p6, p5, p7);
l3 = gmsh.model.geo.addCircleArc(p4, p12, p7);
l4 = gmsh.model.geo.addCircleArc(p3, p9, p6);

#Middle Circle
l5 = gmsh.model.geo.addCircleArc(p8, p1, p10);
l6 = gmsh.model.geo.addCircleArc(p10, p1, p11);
l7 = gmsh.model.geo.addCircleArc(p11, p1, p13);
l8 = gmsh.model.geo.addCircleArc(p13, p1, p8);

#Outer Square
l9  = gmsh.model.geo.addLine(p14, p15)
l10 = gmsh.model.geo.addLine(p15, p16)
l11 = gmsh.model.geo.addLine(p16, p17)
l12 = gmsh.model.geo.addLine(p17, p14)

#4lines
l13 = gmsh.model.geo.addLine(p8, p14);
l14 = gmsh.model.geo.addLine(p10, p15);
l15 = gmsh.model.geo.addLine(p11, p16);
l16 = gmsh.model.geo.addLine(p13, p17);

l17 = gmsh.model.geo.addLine(p3, p13);
l18 = gmsh.model.geo.addLine(p4, p8);
l19 = gmsh.model.geo.addLine(p6, p11);
l20 = gmsh.model.geo.addLine(p7, p10);

# Closed loops

cl1 = gmsh.model.geo.addCurveLoop([l1, l18, -l8, -l17]);
s1 = gmsh.model.geo.addPlaneSurface([cl1]);

cl2 = gmsh.model.geo.addCurveLoop([l3, l20, -l5, -l18]);
s2 = gmsh.model.geo.addPlaneSurface([cl2]);

cl3 = gmsh.model.geo.addCurveLoop([l2, l20, l6, -l19]);
s3 = gmsh.model.geo.addPlaneSurface([cl3]);

cl4 = gmsh.model.geo.addCurveLoop([l4, l19, l7, -l17]);
s4 = gmsh.model.geo.addPlaneSurface([cl4]);

cl5 = gmsh.model.geo.addCurveLoop([l5, l14, -l9, -l13]);
s5 = gmsh.model.geo.addPlaneSurface([cl5]);

cl6 = gmsh.model.geo.addCurveLoop([l6, l15, -l10, -l14]);
s6 = gmsh.model.geo.addPlaneSurface([cl6]);

cl7 = gmsh.model.geo.addCurveLoop([l7, l16, -l11, -l15]);
s7 = gmsh.model.geo.addPlaneSurface([cl7]);

cl8 = gmsh.model.geo.addCurveLoop([l8, l13, -l12, -l16]);
s8 = gmsh.model.geo.addPlaneSurface([cl8]);



#Dividing Surfaces to 10pcs
gmsh.model.geo.mesh.setTransfiniteCurve(l1,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l3,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l2,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l4,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l5,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l6,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l7,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l8,  31)

#Out of circles
gmsh.model.geo.mesh.setTransfiniteCurve(l13, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l15, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l16, 31)

#Outside
gmsh.model.geo.mesh.setTransfiniteCurve(l9, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l11, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l12, 31)

gmsh.model.geo.mesh.setTransfiniteCurve(l17, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l18, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l19, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l20, 31)

gmsh.model.geo.mesh.setTransfiniteSurface(s1)
gmsh.model.geo.mesh.setTransfiniteSurface(s2)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setTransfiniteSurface(s4)
gmsh.model.geo.mesh.setTransfiniteSurface(s5)
gmsh.model.geo.mesh.setTransfiniteSurface(s6)
gmsh.model.geo.mesh.setTransfiniteSurface(s7)
gmsh.model.geo.mesh.setTransfiniteSurface(s8)

gmsh.model.geo.mesh.setRecombine(2, s1)
gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)
gmsh.model.geo.mesh.setRecombine(2, s5)
gmsh.model.geo.mesh.setRecombine(2, s6)
gmsh.model.geo.mesh.setRecombine(2, s7)
gmsh.model.geo.mesh.setRecombine(2, s8)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l9, l10, l11, l12],              11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1, l3, l2, l4],                 12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5, s6, s7, s8], 13 , name="fluid")

# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("gradak2.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()