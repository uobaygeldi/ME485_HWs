import gmsh
import sys

# Initialize Gmsh
gmsh.initialize(sys.argv)

# Add a model named "mixed"
gmsh.model.add("mixed")
lc1 = 0.2
lc2 = 0.05
lc3 = 0.1
# Channel
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc1)
p2 = gmsh.model.geo.addPoint(3.0, 0.0, 0.0, lc3)
p3 = gmsh.model.geo.addPoint(8.0, 0.0, 0.0, lc3)
p4 = gmsh.model.geo.addPoint(8.0, 0.8, 0.0, lc3)
p5 = gmsh.model.geo.addPoint(4.0, 1.4, 0.0, lc3)
p6 = gmsh.model.geo.addPoint(0.0, 2.0, 0.0, lc1)

#obstacle
p7 = gmsh.model.geo.addPoint(4.0, 0.2, 0.0, lc2)
p8 = gmsh.model.geo.addPoint(7.0, 0.6, 0.0, lc2)
p9 = gmsh.model.geo.addPoint(6.0, 0.8, 0.0, lc2)
p10 = gmsh.model.geo.addPoint(2.0, 0.7, 0.0, lc2)

#refinement
p11 = gmsh.model.geo.addPoint(1.950847637 , 1.19757818, 0.0, lc3)
p12 = gmsh.model.geo.addPoint(1.9062124069 , 0.208874876, 0.0, lc3)


#lines
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p5)
l5 = gmsh.model.geo.addLine(p5, p6)
l6 = gmsh.model.geo.addLine(p6, p1)
l7 = gmsh.model.geo.addLine(p7, p8)
l8 = gmsh.model.geo.addLine(p8, p9)
l9 = gmsh.model.geo.addLine(p9, p10)
l10 = gmsh.model.geo.addLine(p10, p7)
l11 = gmsh.model.geo.addLine(p5, p11)
l12 = gmsh.model.geo.add_circle_arc(p11, p10, p12)
l13 = gmsh.model.geo.addLine(p12, p2)

#surfaces

cl1 = gmsh.model.geo.addCurveLoop([l1, -l13, -l12, -l11,l5, l6])
cl2 = gmsh.model.geo.addCurveLoop([l11, l12, l13, l2, l3, l4])
cl3 = gmsh.model.geo.addCurveLoop([l7, l8, l9, l10])

s1 = gmsh.model.geo.addPlaneSurface([cl1])
s2 = gmsh.model.geo.addPlaneSurface([cl2, cl3])


# Synchronize the model
gmsh.model.geo.synchronize()

# Define physical groups for solver compatibility
gmsh.model.addPhysicalGroup(1, [l6], 13, name="inflow")
gmsh.model.addPhysicalGroup(1, [l1, l2, l4, l5, l7, l8, l9, l10], 14, name="wall")
gmsh.model.addPhysicalGroup(1, [l3], 15, name="outflow")
gmsh.model.addPhysicalGroup(2, [s1,s2], 16, name="fluid")

# Mesh generation
gmsh.model.mesh.generate(2)
gmsh.write("V2.msh")

# Run GUI
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalize
gmsh.finalize()
