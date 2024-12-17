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
p0 = gmsh.model.geo.addPoint(0,0,0, lc)
p1  = gmsh.model.geo.addPoint(2.0, 0.0, 0.0, lc);
p2  = gmsh.model.geo.addPoint(4*math.sqrt(2)/2 ,4*math.sqrt(2)/2, 0,lc);
p3  = gmsh.model.geo.addPoint( 0,2, 0,lc);
p4  = gmsh.model.geo.addPoint(-4*math.sqrt(2)/2 ,4*math.sqrt(2)/2, 0,lc);
p5  = gmsh.model.geo.addPoint( -2,0, 0,lc);
p6  = gmsh.model.geo.addPoint(-4*math.sqrt(2)/2 ,-4*math.sqrt(2)/2, 0,lc);
p7  = gmsh.model.geo.addPoint( 0,-2, 0,lc);
p8  = gmsh.model.geo.addPoint(4*math.sqrt(2)/2 ,-4*math.sqrt(2)/2, 0,lc);

p9  = gmsh.model.geo.addPoint( 4,0, 0,lc);
p10  = gmsh.model.geo.addPoint( 0,4, 0,lc);
p11  = gmsh.model.geo.addPoint( -4,0, 0,lc);
p12  = gmsh.model.geo.addPoint( 0,-4, 0,lc);

p13 = gmsh.model.geo.addPoint( 5,5, 0,lc);
p14 = gmsh.model.geo.addPoint( -5,5, 0,lc);
p15 = gmsh.model.geo.addPoint( -5,-5, 0,lc);
p16 = gmsh.model.geo.addPoint( 5,-5, 0,lc);
# Curves
#Inner Star
l1 = gmsh.model.geo.addLine(p1, p2);
l2 = gmsh.model.geo.addLine(p2, p3);
l3 = gmsh.model.geo.addLine(p3, p4);
l4 = gmsh.model.geo.addLine(p4, p5);
l5 = gmsh.model.geo.addLine(p5, p6);
l6 = gmsh.model.geo.addLine(p6, p7);
l7 = gmsh.model.geo.addLine(p7, p8);
l8 = gmsh.model.geo.addLine(p8, p1);

#Middle square
l9  = gmsh.model.geo.addCircleArc(p2, p0, p4)
l10 = gmsh.model.geo.addCircleArc(p4, p0, p6)
l11 = gmsh.model.geo.addCircleArc(p6, p0, p8)
l12 = gmsh.model.geo.addCircleArc(p8, p0, p2)

#4lines
l13 = gmsh.model.geo.addLine(p2, p13);
l14 = gmsh.model.geo.addLine(p4, p14);
l15 = gmsh.model.geo.addLine(p6, p15);
l16 = gmsh.model.geo.addLine(p8, p16);

#Outer Square
l17 = gmsh.model.geo.addLine(p13, p14);
l18 = gmsh.model.geo.addLine(p14, p15);
l19 = gmsh.model.geo.addLine(p15, p16);
l20 = gmsh.model.geo.addLine(p16, p13);

# Closed loops
cl1 = gmsh.model.geo.addCurveLoop([l2, l3, -l9]);
s1 = gmsh.model.geo.addPlaneSurface([cl1]);

cl2 = gmsh.model.geo.addCurveLoop([l4, l5, -l10]);
s2 = gmsh.model.geo.addPlaneSurface([cl2]);

cl3 = gmsh.model.geo.addCurveLoop([l6, l7, -l11]);
s3 = gmsh.model.geo.addPlaneSurface([cl3]);

cl4 = gmsh.model.geo.addCurveLoop([l8, l1, -l12]);
s4 = gmsh.model.geo.addPlaneSurface([cl4]);

cl5 = gmsh.model.geo.addCurveLoop([l14, -l17, -l13, l9]);
s5 = gmsh.model.geo.addPlaneSurface([cl5]);

cl6 = gmsh.model.geo.addCurveLoop([l15, -l18, -l14, l10]);
s6 = gmsh.model.geo.addPlaneSurface([cl6]);

cl7 = gmsh.model.geo.addCurveLoop([l16, -l19, -l15, l11]);
s7 = gmsh.model.geo.addPlaneSurface([cl7]);

cl8 = gmsh.model.geo.addCurveLoop([l16, l20, -l13, -l12]);
s8 = gmsh.model.geo.addPlaneSurface([cl8]);

#Dividing Surfaces to 10pcs
#Circles
gmsh.model.geo.mesh.setTransfiniteCurve(l1,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l8,  31)
gmsh.model.geo.mesh.setTransfiniteCurve(l12, 31)

gmsh.model.geo.mesh.setTransfiniteCurve(l2, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l9, 31)

gmsh.model.geo.mesh.setTransfiniteCurve(l4, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 31)

gmsh.model.geo.mesh.setTransfiniteCurve(l6, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l11, 31)

#Out of circles
gmsh.model.geo.mesh.setTransfiniteCurve(l13, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l14, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l15, 31)
gmsh.model.geo.mesh.setTransfiniteCurve(l16, 31)
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


gmsh.model.geo.mesh.setRecombine(2, s2)
gmsh.model.geo.mesh.setRecombine(2, s3)
gmsh.model.geo.mesh.setRecombine(2, s4)
gmsh.model.geo.mesh.setRecombine(2, s5)
gmsh.model.geo.mesh.setRecombine(2, s6)
gmsh.model.geo.mesh.setRecombine(2, s7)
gmsh.model.geo.mesh.setRecombine(2, s8)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l17, l18, l19, l20],            11 , name="outer")
gmsh.model.addPhysicalGroup(1, [l1,l2,l3,l4,l5,l6,l7,l8],       12 , name="inner")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5, s6, s7, s8],13 , name="fluid")

# Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("gradak1.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()