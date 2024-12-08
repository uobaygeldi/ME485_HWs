import gmsh
import sys

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
p12 = gmsh.model.geo.addPoint(1, -1, 0, lc)
p13 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p14 = gmsh.model.geo.addPoint(-1, 1, 0, lc)
p15 = gmsh.model.geo.addPoint(-1, -1, 0, lc)

# Hexagon
p6 = gmsh.model.geo.addPoint(1, 0, 0, lc)
p7 = gmsh.model.geo.addPoint(0.5, 0.866, 0, lc)
p8 = gmsh.model.geo.addPoint(-0.5, 0.866, 0, lc)
p9 = gmsh.model.geo.addPoint(-1, 0, 0, lc)
p10 = gmsh.model.geo.addPoint(-0.5, -0.866, 0, lc)
p11 = gmsh.model.geo.addPoint(0.5, -0.866, 0, lc)

# Outer circle
l1 = gmsh.model.geo.addCircleArc(p12, p1, p13)
l2 = gmsh.model.geo.addCircleArc(p13, p1, p14)
l3 = gmsh.model.geo.addCircleArc(p14, p1, p15)
l4 = gmsh.model.geo.addCircleArc(p15, p1, p12)

# Inner circle
l5 = gmsh.model.geo.addCircleArc(p2, p1, p3)
l6 = gmsh.model.geo.addCircleArc(p3, p1, p4)
l7 = gmsh.model.geo.addCircleArc(p4, p1, p5)
l8 = gmsh.model.geo.addCircleArc(p5, p1, p2)

# Hexagon lines
l9 = gmsh.model.geo.addLine(p6, p7)
l10 = gmsh.model.geo.addLine(p7, p8)
l11 = gmsh.model.geo.addLine(p8, p9)
l12 = gmsh.model.geo.addLine(p9, p10)
l13 = gmsh.model.geo.addLine(p10, p11)
l14 = gmsh.model.geo.addLine(p11, p6)

# Hex to outer circle
l15 = gmsh.model.geo.addLine(p11, p12)
l16 = gmsh.model.geo.addLine(p12, p6)
l17 = gmsh.model.geo.addLine(p6, p13)
l18 = gmsh.model.geo.addLine(p7, p13)
l19 = gmsh.model.geo.addLine(p14, p8)
l20 = gmsh.model.geo.addLine(p14, p9)
l21 = gmsh.model.geo.addLine(p9, p15)
l22 = gmsh.model.geo.addLine(p15, p10)

# Hex to inner circle
l23 = gmsh.model.geo.addLine(p10, p5)
l24 = gmsh.model.geo.addLine(p11, p2)
l25 = gmsh.model.geo.addLine(p2, p6)
l26 = gmsh.model.geo.addLine(p3, p6)
l27 = gmsh.model.geo.addLine(p3, p7)
l28 = gmsh.model.geo.addLine(p8, p4)
l29 = gmsh.model.geo.addLine(p9, p4)
l30 = gmsh.model.geo.addLine(p9, p5)

# Surfaces
curve_loops = [
    gmsh.model.geo.addCurveLoop([l4, -l15, -l13, -l22]),
    gmsh.model.geo.addCurveLoop([l1, -l17, -l16]),
    gmsh.model.geo.addCurveLoop([l2, l19, -l10, l18]),
    gmsh.model.geo.addCurveLoop([l3, -l21, -l20]),
    gmsh.model.geo.addCurveLoop([l22, -l12, l21]),
    gmsh.model.geo.addCurveLoop([l15, l16, -l14]),
    gmsh.model.geo.addCurveLoop([l17, -l18, -l9]),
    gmsh.model.geo.addCurveLoop([-l19, l20, -l11]),
    gmsh.model.geo.addCurveLoop([l29, -l28, l11]),
    gmsh.model.geo.addCurveLoop([l29, l7, -l30]),
    gmsh.model.geo.addCurveLoop([l30, -l23, -l12]),
    gmsh.model.geo.addCurveLoop([l13, l24, -l8, -l23]),
    gmsh.model.geo.addCurveLoop([l14, -l25, -l24]),
    gmsh.model.geo.addCurveLoop([-l26, -l5, l25]),
    gmsh.model.geo.addCurveLoop([l9, -l27, l26]),
    gmsh.model.geo.addCurveLoop([l10, l28, -l6, l27])
]

surfaces = [gmsh.model.geo.addPlaneSurface([cl]) for cl in curve_loops]

# Meshing
for l in range(1, 31):
    gmsh.model.geo.mesh.setTransfiniteCurve(l, 31)

for s in [5, 6, 7, 8, 9, 11, 12, 13, 15, 16]:
    gmsh.model.geo.mesh.setTransfiniteSurface(surfaces[s - 1])
    gmsh.model.geo.mesh.setRecombine(2, surfaces[s - 1])

gmsh.model.geo.synchronize()

# Physical groups
gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4], 11, name="outer")
gmsh.model.addPhysicalGroup(1, [l5, l6, l7, l8], 12, name="inner")
gmsh.model.addPhysicalGroup(2, surfaces, 13, name="fluid")

# Generate mesh and save
gmsh.model.mesh.generate(2)
gmsh.write("gradok1.msh")

# Launch GUI
if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
