import gmsh
import sys

# Initialize Gmsh
gmsh.initialize(sys.argv)

# Add a model named "mixed"
gmsh.model.add("mixed")
lc = 0.025

# Points
p1  = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
p2  = gmsh.model.geo.addPoint(0.25, 0.0, 0.0, lc)
p3  = gmsh.model.geo.addPoint(0.5, 0.0, 0.0, lc)
p4  = gmsh.model.geo.addPoint(0.75, 0.0, 0.0, lc)
p5  = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, lc)
p6  = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)
p7  = gmsh.model.geo.addPoint(0.75, 1.0, 0.0, lc)
p8  = gmsh.model.geo.addPoint(0.5, 1.0, 0.0, lc)
p9  = gmsh.model.geo.addPoint(0.25, 1.0, 0.0, lc)
p10 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)


# Lines
l1  = gmsh.model.geo.addLine(p1,  p2)
l2  = gmsh.model.geo.addLine(p2,  p3)
l3  = gmsh.model.geo.addLine(p3,  p4)
l4  = gmsh.model.geo.addLine(p4,  p5)
l5  = gmsh.model.geo.addLine(p5,  p6)
l6  = gmsh.model.geo.addLine(p6,  p7)
l7  = gmsh.model.geo.addLine(p7,  p8)
l8  = gmsh.model.geo.addLine(p8,  p9)
l9  = gmsh.model.geo.addLine(p9,  p10)
l10 = gmsh.model.geo.addLine(p10, p1)
l11 = gmsh.model.geo.addLine(p2, p9)
l12 = gmsh.model.geo.addLine(p3, p8)
l13 = gmsh.model.geo.addLine(p4, p7)


# Surfaces
cl1 = gmsh.model.geo.addCurveLoop([l1, l11, l9, l10])
s1 = gmsh.model.geo.addPlaneSurface([cl1])
cl2 = gmsh.model.geo.addCurveLoop([l2, l12, l8, -l11])
s2 = gmsh.model.geo.addPlaneSurface([cl2])
cl3 = gmsh.model.geo.addCurveLoop([l3, l13, l7, -l12])
s3 = gmsh.model.geo.addPlaneSurface([cl3])
cl4 = gmsh.model.geo.addCurveLoop([l4, l5, l6, -l13])
s4 = gmsh.model.geo.addPlaneSurface([cl4])

gmsh.model.geo.mesh.setTransfiniteCurve(l1, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l11, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(l9, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l10, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(l12, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 10)
gmsh.model.geo.mesh.setTransfiniteCurve(l13, 40)
gmsh.model.geo.mesh.setTransfiniteCurve(l7, 10)


# Transfinite and recombination
gmsh.model.geo.mesh.setTransfiniteSurface(s1)
gmsh.model.geo.mesh.setRecombine(2, s1)
gmsh.model.geo.mesh.setTransfiniteSurface(s3)
gmsh.model.geo.mesh.setRecombine(2, s3)

# Synchronize and physical groups
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10], 11, name="outer")
gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4], 13, name="fluid")

# Mesh generation
gmsh.model.mesh.generate(2)
gmsh.write("advection.msh")

# Run GUI
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalize
gmsh.finalize()
