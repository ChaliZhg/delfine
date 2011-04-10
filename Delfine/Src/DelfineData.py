#############################################################
# File: DelfineData.py
# Function: Defines the data variables to be used along the program. Defines also the class that
# will allow data transfer among the several objects of the program.
# Author: Bruno Luna
# Date: 13/02/11
# Modifications Date:
# 16/02/11 - Data structure altered. Now the classes, methods and attributes correspond to the
# ones found in the ./XML/delfineGrammar.rnc file. This should make this file a whole more
# complicated, but should make the rest of the program easier to deal with.
#
# 31/03/11 - Added a 'rocks' list to the RockData class in order to store all the rocks created 
# with the rockType type
#
############################################################
from dolfin import *
# 5th level-------------------------------------------------------------------------------------------------------------------------------------
class Viscosity:
    """5th level - Class to define viscosity model and value"""
    def __init__(self):
        self.model = None
        self.value = None
class Density:
    """5th level - Class to define density model and value"""
    def __init__(self):
        self.model = None
        self.value = None
class CapillaryP:
    """5th level - Class to define capillary pressure model and value"""
    def __init__(self):
        self.model = None
        self.value = None
class Porosity:
    """5th level - Class to define porosity model and value"""
    def __init__(self):
        self.compressible = None
        self.value = None
class RockHeatCoeff:
    """5th level - Class to define rock conductive heat transfer coefficient"""
    def __init__(self):
        self.value = None
class Permeability:
    """5th level - Class to define rock flow permeability"""
    def __init__(self):
        self.type = None
        # if type == per-domain
        self.K = [] # 6th level
        # elif type == per-element-list
        self.filename = None
class Kr:
    """5th level - Class to define the phases relative permeability parameters"""
    def __init__(self):
        self.krEnd = None # 6th level
        self.Sr = None # 6th level
        self.n = None # 6th level

# 4th level-------------------------------------------------------------------------------------------------------------------------------------
class DolfinGen:
    """4th level - Class for dolfin generated mesh data"""
    def __init__(self):
        self.type = None
        self.n = []
        self.p0 = []
        self.p1 = []
class SubDom:
    """4th level - Class to define the subdomains"""
    def __init__(self):
        self.quantity = None
class Well:
    """4th level - Class to define the wells"""
    def __init__(self):
        self.id = None
        self.bctype = None
        self.function = None
        self.value = None
class Fluid:
    """4th level - Class to define fluid data(w,o,g)"""
    def __init__(self):
        self.use = None # E.g gas.Fluid.use = no or None for 2-phase o-w
        self.viscosity = Viscosity()
        self.density = Density()
        self.capillaryP = CapillaryP()
class RockType:
    """4th level - Class to define the rock type"""
    def __init__(self):
        self.id = None
        self.porosity = Porosity()
        self.rockHeatCoeff = RockHeatCoeff()
        self.permeability = Permeability()
class RelativePerm:
    """4th level - Class to define the relative permeability model and data"""
    def __init__(self):
        self.type = None
        self.krw = Kr()
        self.kro = Kr()
        self.krg = Kr()
class PreConditioning:
    """4th level - Class to define the preconditioner for the pressure solver"""
    def __init__(self):
        self.type = None
        self.numCoarseLevel = None # 5th level
        self.numRelaxIter = None # 5th level
# 3nd level-------------------------------------------------------------------------------------------------------------------------------------
class MeshData:
    """3rd level - Class for mesh information"""
    def __init__(self):
        self.dim = None
        self.type = None
        # if type = dolfin-generated
        self.dolfinGen = DolfinGen()
        # elif type = gmsh | xml
        self.filename = None
        self.subdom = SubDom()
class BcData:
    """3rd level - Class for boundary conditions information"""
    def __init__(self):
        self.wellType = Well()
        self.wells = []
class FluidData:
    """3rd level - Class for fluids properties"""
    def __init__(self):
        self.water = Fluid()
        self.oil = Fluid()
        self.gas = Fluid()
class RockData:
    """3rd level - Class for rock properties"""
    def __init__(self):
        self.rockType = RockType()
        self.rocks = []
class RockFluidData: 
    """3rd level - Class for rock-fluid properties"""
    def __init__(self):
        self.relativePerm = RelativePerm()
class PressSolvData:
    """3rd level - Class for pressure solver parameters"""
    def __init__(self):
        self.type = None
        self.tolerance = None #4th level
        self.maxNumSteps = None #4th level
        self.preConditioning = PreConditioning()
class SaturSolvData:
    """3rd level - Class for saturation solver parameters"""
    def __init__(self):
        self.totalTime = None # 4th level
        self.courant = None # 4th level
        self.limiterType = None # 4th level
# 2nd level-------------------------------------------------------------------------------------------------------------------------------------
class GeomData:
    """2nd level - Class for geometrical information"""
    def __init__(self):
        self.mesh = MeshData()
        self.bc = BcData()
class PhysData:
    """2nd level - Class for physical information"""
    def __init__(self):
        self.fluid = FluidData()
        self.rock = RockData()
        self.rockFluid = RockFluidData()
class NumData:
    """2nd level - Class for numerical information"""
    def __init__(self):
        self.pressSolv = PressSolvData()
        self.saturSolv = SaturSolvData()
# 1st level-------------------------------------------------------------------------------------------------------------------------------------
class DelfineData:
    """1st level - Class for information transfer among objects"""
     
    def __init__(self, type):
        if (type == "transfer"):
            self.A = None
            self.rhs = None
            self.x = None
            self.residuals = None
            self.V = None
            self.u0 = None
        elif (type == "input"):
            # Geometric
            self.geom = GeomData()
            # Physical
            self.phys = PhysData()
            # Numerical
            self.num = NumData()
#############################################################
