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
# 31/03/11 - Added a 'rocks' list to the RockData class in order to store all the rocks created 
# with the rockType type
# 28/05/11 - Permeability Tensor class moved in to this module for code organization
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
        self.order = None
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
        self.formulation = None
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
            self.K = None
            self.mob = None
            self.velocity = None
        elif (type == "input"):
            # Geometric
            self.geom = GeomData()
            # Physical
            self.phys = PhysData()
            # Numerical
            self.num = NumData()
 
 # Auxiliary class for elliptic assembly module------------------------------------------------------------------------------------------           
class PermeabilityTensor3D(Expression):
    """Defines the 3D permeability tensor for each cell"""
    def __init__(self, mesh, parameter):
        self.mesh = mesh
        self.rocks = parameter.phys.rock.rocks
        self.DomainID = []
        self.K =[]
        # Getting rocks IDs and respective K from parameter structure
        for i in range(len(self.rocks)):
            if (self.rocks[i].permeability.type == 'per-domain'):
                self.DomainID.append(int(self.rocks[i].id))
                self.K.append(self.rocks[i].permeability.K)
    def eval_cell(self, values, x, ufc_cell):
        # Get material indicator(which corresponds to rock ID) from mesh with mesh_function
        mf = self.mesh.data().mesh_function("material_indicators")
        i = ufc_cell.index
        j = 0
        if (type(i) is int):
            if (mf != None):
                for id in self.DomainID:
                    if (mf[i] == id):
                        values[0] = self.K[j][0] # Kxx
                        values[1] = self.K[j][1] # Kxy
                        values[2] = self.K[j][2] # Kxz            | Kxx  Kxy  Kxz |       | values[0] values[1] values[2] |
                        values[3] = values[1]   # Kyx     K = | Kyx  Kyy  Kyz  | => | values[3] values[4] values[5] |
                        values[4] = self.K[j][3] # Kyy            | Kzx  Kzy  Kzz  |       | values[6] values[7] values[8] |
                        values[5] = self.K[j][4] # Kyz
                        values[6] = values[2]   # Kzx
                        values[7] = values[5]   # Kzy
                        values[8] = self.K[j][5] # Kzz
                    j += 1
            else: 
                # For dolfin-generated meshes or meshes without "material_indicator"
                # This option consider just homogeneous cases
                values[0] = self.K[0][0] # Kxx
                values[1] = self.K[0][1] # Kxy
                values[2] = self.K[0][2] # Kxz            | Kxx  Kxy  Kxz |       | values[0] values[1] values[2] |
                values[3] = values[1]    # Kyx     K = | Kyx  Kyy  Kyz  | => | values[3] values[4] values[5] |
                values[4] = self.K[0][3] # Kyy            | Kzx  Kzy  Kzz  |       | values[6] values[7] values[8] |
                values[5] = self.K[0][4] # Kyz
                values[6] = values[2]    # Kzx
                values[7] = values[5]    # Kzy
                values[8] = self.K[0][5] # Kzz
        else:
            pass
    def value_shape(self):
        return (3, 3)
# Auxiliary class for elliptic assembly module------------------------------------------------------------------------------------------           
class InversePermeabilityTensor3D(Expression):
    """Defines the 3D inverse permeability tensor for each cell. Used for MixedFEM"""
    def __init__(self, mesh, parameter):
        self.mesh = mesh
        self.rocks = parameter.phys.rock.rocks
        self.DomainID = []
        self.K =[]
        self.detK =[]
        # Auxiliary variables for matrix inversion
        # (Used convention of en.wikipedia.org/wiki/Invertible_matrix)
        self.A = []
        self.B = []
        self.C = []
        self.D = []
        self.E = []
        self.F = []
        self.G = []
        self.H = []
        # Getting rocks IDs and respective K from parameter structure
        for i in range(len(self.rocks)):
            if (self.rocks[i].permeability.type == 'per-domain'):
                self.DomainID.append(int(self.rocks[i].id))
                self.K.append(self.rocks[i].permeability.K)
                # Rock permeabilities
                Kxx = self.K[i][0]
                Kxy = self.K[i][1]
                Kxz = self.K[i][2]
                Kyy = self.K[i][3]
                Kyz = self.K[i][4]
                Kzz = self.K[i][5]
                # Permability matrix determinant
                self.detK.append(Kxx*(Kyy*Kzz - Kyz*Kyz) + \
                Kxy*(Kyz*Kxz - Kzz*Kxy) + \
                Kxz*(Kxy*Kyz - Kyy*Kxz))
                # Cofactors term for inverse matrix (see wikipedia ref. above)
                self.A.append(Kyy*Kzz - Kyz*Kyz)
                self.B.append(Kyz*Kxz - Kzz*Kxy)
                self.C.append(Kxy*Kyz - Kyy*Kxz)
                self.D.append(Kxz*Kyz - Kxy*Kzz)
                self.E.append(Kxx*Kzz - Kxz*Kxz)
                self.F.append(Kxz*Kxy - Kxx*Kyz)
                self.G.append(Kxy*Kyz - Kxz*Kyy)
                self.H.append(Kxz*Kxy - Kxx*Kyz)
                self.K.append(Kxx*Kyy - Kxy*Kxy)
    def eval_cell(self, values, x, ufc_cell):
        # Get material indicator(which corresponds to rock ID) from mesh with mesh_function
        mf = self.mesh.data().mesh_function("material_indicators")
        i = ufc_cell.index
        j = 0
        if (type(i) is int):
            if (mf != None):
                for id in self.DomainID:
                    if (mf[i] == id):
                        invDetK = 1/(self.detK[j])
                        values[0] = invDetK*self.A[j]
                        values[1] = invDetK*self.D[j]
                        values[2] = invDetK*self.G[j]
                        values[3] = invDetK*self.B[j]    #    invK = (1/detK) * | A  D  G| =>  | values[0] values[1] values[2]| 
                        values[4] = invDetK*self.E[j]    #                                | B  E  H| =>  | values[3] values[4] values[5]| 
                        values[5] = invDetK*self.H[j]    #                                | C  F  K| =>  | values[6] values[7] values[8]| 
                        values[6] = invDetK*self.C[j]   
                        values[7] = invDetK*self.F[j]   
                        values[8] = invDetK*self.K[j] 
                    j += 1
            else: 
                # For dolfin-generated meshes or meshes without "material_indicator"
                # This option consider just homogeneous cases
                values[0] = self.K[0][0] 
                values[1] = self.K[0][1] 
                values[2] = self.K[0][2]     
                values[3] = values[1]    
                values[4] = self.K[0][3]        
                values[5] = self.K[0][4] 
                values[6] = values[2]    
                values[7] = values[5]    
                values[8] = self.K[0][5] 
        else:
            pass
    def value_shape(self):
        return (3, 3)
 # Auxiliary class for elliptic assembly module------------------------------------------------------------------------------------------  
class PermeabilityTensor2D(Expression):
    """Defines the 2D permeability tensor for each cell"""
    def __init__(self, mesh, parameter):
        self.mesh = mesh
        self.rocks = parameter.phys.rock.rocks
        self.DomainID = []
        self.K =[]
        # Getting rocks IDs and respective K from parameter structure
        for i in range(len(self.rocks)):
            if (self.rocks[i].permeability.type == 'per-domain'):
                self.DomainID.append(int(self.rocks[i].id))
                self.K.append(self.rocks[i].permeability.K)
    def eval_cell(self, values, x, ufc_cell):
        # Get material indicator(which corresponds to rock ID) from mesh with mesh_function
        mf = self.mesh.data().mesh_function("material_indicators")
        i = ufc_cell.index
        j = 0
        if (type(i) is int):
            if (mf != None):
                for id in self.DomainID:
                    if (mf[i] == id):
                        values[0] = self.K[j][0] # Kxx
                        values[1] = self.K[j][1] # Kxy     K = | Kxx  Kxy | => | values[0] values[1] |
                        values[2] = values[1]   # Kyx            | Kyx  Kyy  |       | values[2] values[3] |      
                        values[3] = self.K[j][2] # Kyy
                    j += 1
            else: 
                # For dolfin-generated meshes or meshes without "material_indicator"
                # This option consider just homogeneous cases
                values[0] = self.K[0][0] # Kxx
                values[1] = self.K[0][1] # Kxy
                values[2] = values[1]    # Kyx
                values[3] = self.K[0][2] # Kyy
        else:
            pass
    def value_shape(self):
        return (2, 2)
# Auxiliary class for elliptic assembly module------------------------------------------------------------------------------------------  
class InversePermeabilityTensor2D(Expression):
    """Defines the 2D inverse permeability tensor for each cell. Used for MixedFEM"""
    def __init__(self, mesh, parameter):
        self.mesh = mesh
        self.rocks = parameter.phys.rock.rocks
        self.DomainID = []
        self.K = []
        self.detK =[]
        # Getting rocks IDs and respective K from parameter structure
        for i in range(len(self.rocks)):
            if (self.rocks[i].permeability.type == 'per-domain'):
                self.DomainID.append(int(self.rocks[i].id))
                self.K.append(self.rocks[i].permeability.K)
                Kxx = self.K[i][0]
                Kxy = self.K[i][1]
                Kyy = self.K[i][2]
                self.detK.append(Kxx*Kyy - Kxy*Kxy)
    def eval_cell(self, values, x, ufc_cell):
        # Get material indicator(which corresponds to rock ID) from mesh with mesh_function
        mf = self.mesh.data().mesh_function("material_indicators")
        i = ufc_cell.index
        j = 0
        if (type(i) is int):
            if (mf != None):
                for id in self.DomainID:
                    if (mf[i] == id):
                        invDetK = 1/(self.detK[j])
                        values[0] = invDetK*self.K[j][2]
                        values[1] = invDetK*(-self.K[j][1])  #    invK = (1/detK) * | Kyy  -Kxy | =>  | values[0] values[1] |
                        values[2] = values[1]                     #                               | -Kyx  Kxx  |       | values[2] values[3] |      
                        values[3] = invDetK*self.K[j][0]      # 
                    j += 1
            else: 
                # For dolfin-generated meshes or meshes without "material_indicator"
                # This option consider just homogeneous cases
                        invDetK = 1/(self.detK[0])
                        values[0] = invDetK*self.K[0][2]
                        values[1] = invDetK*(-self.K[0][1])  #    invK = (1/detK) * | Kyy  -Kxy | =>  | values[0] values[1] |
                        values[2] = values[1]                      #                               | -Kyx  Kxx  |       | values[2] values[3] |      
                        values[3] = invDetK*self.K[0][0]      # 
        else:
            pass
    def value_shape(self):
        return (2, 2)
# Auxiliary class for elliptic assembly module------------------------------------------------------------------------------------------ 
# Do not confuse with the Viscosity class used also in this module just for data transfer.
class ViscScalar(Expression):
    """Defines the constant viscosity for each point"""
    def __init__(self, fluid, parameter):
        self.fluid = fluid
        if (self.fluid == "water"):
            self.mu =  parameter.phys.fluid.water.viscosity.value
        elif (self.fluid == "oil"):
            self.mu =  parameter.phys.fluid.oil.viscosity.value
    def eval (self, values, x):
        values[0] = self.mu
# Auxiliary class for elliptic assembly module------------------------------------------------------------------------------------------ 
class RelatPermScalar(Expression):
    """Defines the relative permeability for each point"""
    def __init__(self, fluid, parameter, Sw):
        self.fluid = fluid
        self.model = parameter.phys.rockFluid.relativePerm.type
        self.Sw = Sw
        # Read data for Corey model (Attention: other models(Stone's I, II) still pending)
        if (self.model == "corey"):
            if (self.fluid == "water"):
                self.Swr =  parameter.phys.rockFluid.relativePerm.krw.Sr
                self.krwEnd = parameter.phys.rockFluid.relativePerm.krw.krEnd
                self.nw = parameter.phys.rockFluid.relativePerm.krw.n
            elif (self.fluid == "oil"):
                self.Sor =  parameter.phys.rockFluid.relativePerm.kro.Sr
                self.no = parameter.phys.rockFluid.relativePerm.kro.n
    def eval (self, values, x):
        if (self.model == "corey"):
            if (self.fluid == "water"):
                values[0] = self.Sw**4 #FIXME: Fix corey's formulae
            elif (self.fluid == "oil"):
                values[0] = (1 - self.Sw)**4
#############################################################
