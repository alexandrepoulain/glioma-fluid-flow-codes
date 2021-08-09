import SVMTK as svmtk

def create_brain_mesh(stls, output,
                      resolution=32, remove_ventricles=True):

    # Load each of the Surfaces
    surfaces = [svmtk.Surface(stl) for stl in stls]
    
    # Take the union of the left (#3) and right (#4)
    # white surface and put the result into
    # the (former left) white surface
    surfaces[2].union(surfaces[3])

    # ... and drop the right white surface from the list
    surfaces.pop(3)

    # Define identifying tags for the different regions 
    tags = {"pial": 1, "tumor": 2, "white": 3, "ventricle": 4}

    # Label the different regions
    smap = svmtk.SubdomainMap()
    smap.add("10000", tags["pial"]) 
    smap.add("01000", tags["pial"]) 
    smap.add("10100", tags["white"])
    smap.add("01010", tags["white"])
    smap.add("11010", tags["white"])
    smap.add("10011", tags["ventricle"])
    smap.add("01011", tags["ventricle"])
    smap.add("11011", tags["ventricle"])
    smap.add("01100", tags["tumor"])

    # Generate mesh at given resolution
    domain = svmtk.Domain(surfaces, smap)
    domain.create_mesh(resolution)

    # Remove ventricles perhaps
    #if remove_ventricles:
    #    domain.remove_subdomain(tags["ventricle"])

    # Save mesh    
    domain.save(output)

stls = ("lh.pial.stl", "rh.pial.stl",
        "lh.white.stl", "rh.white.stl",
        "lh.ventricles.stl","lh.tumor-surface.stl")


create_brain_mesh(stls, "neil-brain.mesh")













