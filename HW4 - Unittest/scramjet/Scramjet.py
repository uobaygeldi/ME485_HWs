import gmsh
import sys

# Initialize Gmsh
gmsh.initialize(sys.argv)

# Add a model named "mixed"
gmsh.model.add("mixed")
lc1 = 0.1
lc2 = 0.1

#channel
p1  = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc1)
p2 = gmsh.model.geo.addPoint(8.0, 0.0, 0.0, lc1)
p3 = gmsh.model.geo.addPoint(8.0, 0.8, 0.0, lc1)
p4 = gmsh.model.geo.addPoint(0.0, 2.0, 0.0, lc1)

#obstacle
p5 = gmsh.model.geo.addPoint(4.0, 0.2, 0.0, lc2)
p6 = gmsh.model.geo.addPoint(7.0, 0.6, 0.0, lc2)
p7 = gmsh.model.geo.addPoint(6.0, 0.7, 0.0, lc2)
p8 = gmsh.model.geo.addPoint(2.0, 0.7, 0.0, lc2)


#lines
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
l5 = gmsh.model.geo.addLine(p5, p6)
l6 = gmsh.model.geo.addLine(p6, p7)
l7 = gmsh.model.geo.addLine(p7, p8)
l8 = gmsh.model.geo.addLine(p8, p5)



# Surfaces
cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])

s1 = gmsh.model.geo.addPlaneSurface([cl1, cl2])



# Synchronize the model
gmsh.model.geo.synchronize()

# Define physical groups for solver compatibility
gmsh.model.addPhysicalGroup(1, [l4], 13, name="inflow")
gmsh.model.addPhysicalGroup(1, [l1, l3, l5, l6, l7, l8], 14, name="wall")
gmsh.model.addPhysicalGroup(1, [l2], 15, name="outflow")
gmsh.model.addPhysicalGroup(2, [s1], 16, name="fluid")

# Mesh generation
gmsh.model.mesh.generate(2)
gmsh.model.mesh.optimize("Laplace2D")
gmsh.write("V1.msh")

# Run GUI
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalize
gmsh.finalize()
