import gmsh
import sys
import math

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize(sys.argv)

# Next add a new model named "cavity" 
gmsh.model.add("mixed")
lc = 0.1

#Point 1 start
p1  = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc);


# BUMP COORDINATES
p2  = gmsh.model.geo.addPoint(2.0, 0.0, 0.0, 0.05);
p3 =  gmsh.model.geo.addPoint(3, -6.210, 0.0, 0.05); # The center point which ensures passing from  (3, 0.08, 0)
p4 =  gmsh.model.geo.addPoint(4.0, 0.0, 0.0, 0.05);


p5 = gmsh.model.geo.addPoint(6.0, 0.0, 0.0, lc);

p6  = gmsh.model.geo.addPoint(6.0, 1.1, 0.0, lc);

p7  = gmsh.model.geo.addPoint(0.0, 1.1, 0.0, lc);



l1  = gmsh.model.geo.addLine(p1,  p2) # Lower channel until  bump
l2  = gmsh.model.geo.addCircleArc(p2, p3, p4) #   Bump
l3  = gmsh.model.geo.addLine(p4,  p5) # The line after bump

l4  = gmsh.model.geo.addLine(p5,  p6) # Outlet
l5  = gmsh.model.geo.addLine(p6,  p7) # Upper Channel


l6  = gmsh.model.geo.addLine(p7,  p1) # Inlet


cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6]); # Based on outer bump


s1 = gmsh.model.geo.addPlaneSurface([cl1]) # The region



gmsh.model.geo.mesh.setTransfiniteCurve(l1, 60)

gmsh.model.geo.mesh.setTransfiniteCurve(l2, 60)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, 60)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, 60)
gmsh.model.geo.mesh.setTransfiniteCurve(l5, 60)
gmsh.model.geo.mesh.setTransfiniteCurve(l6, 60)



#gmsh.model.geo.mesh.setTransfiniteSurface(s1)

#gmsh.model.geo.mesh.setRecombine(2, s1)



gmsh.model.geo.synchronize()


gmsh.model.addPhysicalGroup(1, [l6], 13 , name="inflow")
gmsh.model.addPhysicalGroup(1, [l4], 14 , name="outflow")


gmsh.model.addPhysicalGroup(1, [l1,l2,l3,l5], 15 , name="wall")

gmsh.model.addPhysicalGroup(2, [s1],9 , name="fluid")



# # Save it to disk
gmsh.model.mesh.generate(2)
gmsh.write("bump_supersonic1.msh")
# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
# This should be called when you are done using the Gmsh Python API:
gmsh.finalize()