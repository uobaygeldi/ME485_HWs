import gmsh
import sys

# Gmsh'i başlat
gmsh.initialize(sys.argv)

# Yeni bir model oluştur
gmsh.model.add("mixed")
lc = 0.02  # Eleman boyutu

# Dış kare noktaları
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, lc)
p3 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)
p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)

# Dış kare çizgileri
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Dış kare yüzeyi
cl_outer = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
s_outer = gmsh.model.geo.addPlaneSurface([cl_outer])

# Geometriyi senkronize et
gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4], 11 , name="outer")
gmsh.model.addPhysicalGroup(2, [s_outer],13 , name="fluid")

gmsh.model.mesh.generate(2)  # 2D mesh (otomatik olarak üçgen elemanlarla mesh oluşturur)

# Mesh dosyasını kaydet
gmsh.write("advection.msh")

# Sonuçları görüntüle
if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

# Gmsh işlemini sonlandır
gmsh.finalize()