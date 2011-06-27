#############################################################
# File: DelfineInput.py
# Function: Reads the .xml input file with all the physical and numerical parameters. Parses the
# file using ElemTree and validates it using RNV and the pre-defined grammar. Calls also the 
# script which reads the geometrical file containing the mesh.
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 27/03/11 : Read the rest of the input .xml file in the data structure.
# 30/03/11 : Created a check.log file to check against file read errors
# 31/03/11 : Created a check (if condition) to test for the not compulsory data. If they are not
# present in the input file the program associates None to them and continues the workflow.
# These not compulsory data type are marked with '?' in the XML/delfineGrammar.rnc file.
# 27/04/11 : Read the Kxz permeability data that was missing. Also test to check if all the
# permeability tensor data are present for the given problem dimension
#
############################################################
import sys, os
from xml.etree import ElementTree as ET
from dolfin import *
from DelfineData import Well # Import this class locally for read purposes (see Dev. Diary for details)
from DelfineData import RockType # Import this class locally for read purposes (see Dev. Diary for details)

class DelfineInput:
     """Reads all the input data"""
     
     def __init__(self):
         pass
    
     def readData(self, parameter, argv):
        """Parses the 'casename'.xml file and stores data in delfineData"""
        inputFile = "./CaseFiles/" + argv[1] + ".xml"
        
        # Parses 'caseName'.xml File
        try:
            tree = ET.parse(inputFile)
        except Exception, inst:
            print "Error(2):"   
            print "Unexpected error opening %s: %s" % (inputFile, inst)
            print " "
            sys.exit(1)
        
        # Validate the input file against delfine grammar (very important to avoid
        # a lot of possible input errors)
        cmd = 'rnv -q ' + './XML/delfineGrammar.rnc ' + inputFile
        if (os.system(cmd) != 0):
            print "Error(3):"   
            print "Input file failed check against grammar ./XML/delfineGrammar.rnc"
            print "This may indicate that the file was bad generated "
            print " "
            sys.exit(1)
        
        # Creates check.log file. This file is used to check against reading errors
        log = open("Results/check.log",  mode='w')
        
        # 1st level(Root) xml element
        delfineXML = tree.getroot()
        
        # 2nd level xml elements
        geom = delfineXML.find('geometry')
        phys = delfineXML.find('physical')
        num = delfineXML.find('numerical')
        
        # 3rd level and on xml elements
        
        # # Geometrical parameters - geom
        mesh = geom.find('mesh')
        bc = geom.find('boundary-conditions')
        
        # Reads data from xml elements and copy it to parameters
        # attributes defined in DelfineData.py
        
        # ## Mesh - geom.mesh
        parameter.geom.mesh.type = mesh.attrib['type']
        parameter.geom.mesh.dim = int(mesh.attrib['dimension'])
        parameter.geom.mesh.order = int(mesh.attrib['order'])
        
        meshType = parameter.geom.mesh.type
        meshDim = parameter.geom.mesh.dim
        meshOrder = parameter.geom.mesh.order
        
        # Retrives mesh type and number of divisions, in case of a 
        # structured mesh generated by dolfin, or just take the
        # filename in case of a gmsh or xml mesh
        
        # ### Dolfin-Generated - geom.mesh.dolfinGen
        if (meshType == "dolfin-generated"):
            dolfinGen = mesh.find('dolfin-generated')
            parameter.geom.mesh.dolfinGen.type = dolfinGen.attrib['type']
            dolfinGenType = parameter.geom.mesh.dolfinGen.type
            if (dolfinGenType == "UnitInterval"):
                parameter.geom.mesh.dolfinGen.n = [int(dolfinGen.attrib['nx'])]
            elif (dolfinGenType == "UnitSquare"):
                parameter.geom.mesh.dolfinGen.n = [int(dolfinGen.attrib['nx']),  int(dolfinGen.attrib['ny'])]
            elif (dolfinGenType == "UnitCube"):
                parameter.geom.mesh.dolfinGen.n = [int(dolfinGen.attrib['nx']),  int(dolfinGen.attrib['ny']),  \
                int(dolfinGen.attrib['nz'])]
            elif (dolfinGenType == "UnitCircle"):
                parameter.geom.mesh.dolfinGen.n = [int(dolfinGen.attrib['nr'])]
            elif (dolfinGenType == "Rectangle"):
                parameter.geom.mesh.dolfinGen.p0 = [float(dolfinGen.attrib['x0']),  float(dolfinGen.attrib['y0'])]
                parameter.geom.mesh.dolfinGen.p1 = [float(dolfinGen.attrib['x1']),  float(dolfinGen.attrib['y1'])]
                parameter.geom.mesh.dolfinGen.n = [int(dolfinGen.attrib['nx']),  int(dolfinGen.attrib['ny'])]
        # Gmsh or Xml file      
        elif (meshType == "gmsh") | (meshType == "xml"):
            parameter.geom.mesh.filename = mesh.findtext('filename')
        else:
            print "Error(4):"
            print "No valid mesh format selected in the .xml input file"
            print " "
            sys.exit(1)
        
        # Writing in .log file for check purposes
        log.write("The Type of Mesh is: " + str(meshType) +  "\n")
        log.write("The Dimension of Mesh is: " + str(meshDim) +  "\n")
        log.write("The Order of Mesh is: " + str(meshOrder) +  "\n")
        
        # ## Boundary conditions - geom.bcgeom.bc
        etWells = bc.findall('well')
        log.write("The Number of Wells is: " + str(len(etWells)) +  "\n")
        w =[] # Dummy list to store well informations
        for i in range(len(etWells)):
            w.append(Well())
            w[i].id = etWells[i].attrib["id"]
            w[i].function = etWells[i].attrib["function"]
            w[i].bctype = etWells[i].attrib["bctype"]
            w[i].value = etWells[i].text
            parameter.geom.bc.wells.insert(i, w[i])
            # Writing in .log file for check purposes
            log.write("Well ID: " + str(parameter.geom.bc.wells[i].id) +"\n")

        # # Physical parameters - phys
        fluidProp = phys.find('fluid-properties')
        rockProp = phys.find('rock-properties')
        rockFluidProp = phys.find('rock-fluid-properties')
        
        # ## Fluid properties - phys.fluid
        water = fluidProp.find('water')
        oil = fluidProp.find('oil')
        gas = fluidProp.find('gas')
        
        # Read water properties
        parameter.phys.fluid.water.use = water.attrib["use"]
        
        parameter.phys.fluid.water.viscosity.model = water.find('viscosity').attrib["model"]
        parameter.phys.fluid.water.viscosity.value = water.find('viscosity').text
       
        parameter.phys.fluid.water.density.model = water.find('density').attrib["model"]
        parameter.phys.fluid.water.density.value = water.find('density').text
        
        dummyCapP = water.find('capillary-pressure') 
        if (dummyCapP != None): # Checking if not compulsory data is given in input file
            parameter.phys.fluid.water.capillaryP.model = dummyCapP.attrib["model"]
            parameter.phys.fluid.water.capillaryP.value = dummyCapP.text
        
        # Writing in .log file for check purposes
        log.write("Use Water? " + str(parameter.phys.fluid.water.use) + "\n")
        log.write("Water Viscosity (Model | Value): " + str(parameter.phys.fluid.water.viscosity.model) \
        + " | " + str(parameter.phys.fluid.water.viscosity.value) + "\n")
        log.write("Water Density (Model | Value): " + str(parameter.phys.fluid.water.density.model) \
        + " | " + str(parameter.phys.fluid.water.density.value) + "\n")
        log.write("Water CapillaryP (Model | Value): " + str(parameter.phys.fluid.water.capillaryP.model) \
        + " | " + str(parameter.phys.fluid.water.capillaryP.value) + "\n")
        
        # Read oil properties
        if (oil != None): # Checking if not compulsory data is given in input file
            parameter.phys.fluid.oil.use = oil.attrib["use"]
        
            parameter.phys.fluid.oil.viscosity.model = oil.find('viscosity').attrib["model"]
            parameter.phys.fluid.oil.viscosity.value = oil.find('viscosity').text
       
            parameter.phys.fluid.oil.density.model = oil.find('density').attrib["model"]
            parameter.phys.fluid.oil.density.value = oil.find('density').text
        
            dummyCapP = oil.find('capillary-pressure') 
            if (dummyCapP != None): # Checking if not compulsory data is given in input file
                parameter.phys.fluid.oil.capillaryP.model = dummyCapP.attrib["model"]
                parameter.phys.fluid.oil.capillaryP.value = dummyCapP.text
        
        # Writing in .log file for check purposes
        log.write("Use Oil? " + str(parameter.phys.fluid.oil.use) + "\n")
        log.write("Oil Viscosity (Model | Value): " + str(parameter.phys.fluid.oil.viscosity.model) \
        + " | " + str(parameter.phys.fluid.oil.viscosity.value) + "\n")
        log.write("Oil Density (Model | Value): " + str(parameter.phys.fluid.oil.density.model) \
        + " | " + str(parameter.phys.fluid.oil.density.value) + "\n")
        log.write("Oil CapillaryP (Model | Value): " + str(parameter.phys.fluid.oil.capillaryP.model) \
        + " | " + str(parameter.phys.fluid.oil.capillaryP.value) + "\n")
        
        # Read gas properties
        if (gas != None): # Checking if not compulsory data is given in input file
            parameter.phys.fluid.gas.use = gas.attrib["use"]
        
            parameter.phys.fluid.gas.viscosity.model = gas.find('viscosity').attrib["model"]
            parameter.phys.fluid.gas.viscosity.value = gas.find('viscosity').text
       
            parameter.phys.fluid.gas.density.model = gas.find('density').attrib["model"]
            parameter.phys.fluid.gas.density.value = gas.find('density').text
        
            dummyCapP = gas.find('capillary-pressure') 
            if (dummyCapP != None): # Checking if not compulsory data is given in input file
                parameter.phys.fluid.gas.capillaryP.model = dummyCapP.attrib["model"]
                parameter.phys.fluid.gas.capillaryP.value = dummyCapP.text
            
        # Writing in .log file for check purposes
        log.write("Use Gas? " + str(parameter.phys.fluid.gas.use) + "\n")
        log.write("Gas Viscosity (Model | Value): " + str(parameter.phys.fluid.gas.viscosity.model) \
        + " | " + str(parameter.phys.fluid.gas.viscosity.value) + "\n")
        log.write("Gas Density (Model | Value): " + str(parameter.phys.fluid.gas.density.model) \
        + " | " + str(parameter.phys.fluid.gas.density.value) + "\n")
        log.write("Gas CapillaryP (Model | Value): " + str(parameter.phys.fluid.gas.capillaryP.model) \
        + " | " + str(parameter.phys.fluid.gas.capillaryP.value) + "\n")
        
        # ## Rock roperties - phys.rock
        etRockTypes = rockProp.findall('rock-type')
        log.write("The Number of Rock types is: " + str(len(etRockTypes)) +  "\n")
        
        r =[] # Dummy list to store rock-types informations
        for i in range(len(etRockTypes)):
            r.append(RockType())
            r[i].id = etRockTypes[i].attrib["id"]
            r[i].porosity.compressible = etRockTypes[i].find('porosity').attrib["compressible"]
            r[i].porosity.value = etRockTypes[i].find('porosity').text
            dummyHCoef= etRockTypes[i].find('rock-heat-coefficient')
            if (dummyHCoef != None): # Checking if not compulsory data is given in input file
                r[i].rockHeatCoeff.value = dummyHCoef.text
            r[i].permeability.type = etRockTypes[i].find('permeability').attrib["type"]
            if (r[i].permeability.type == "per-domain"):
                Kxx = etRockTypes[i].find('permeability').find('Kxx')
                Kxy = etRockTypes[i].find('permeability').find('Kxy')
                Kxz = etRockTypes[i].find('permeability').find('Kxz')
                Kyy = etRockTypes[i].find('permeability').find('Kyy')
                Kyz = etRockTypes[i].find('permeability').find('Kyz')
                Kzz = etRockTypes[i].find('permeability').find('Kzz')
                
                r[i].permeability.K.append(float(Kxx.text))
                # Checking if not compulsory data is given in input file (as only Kxx is alway mandatory)
                # For 2-D:
                # K = | Kxx Kxy |   Kxy = Kyx
                #        | Kyx Kyy |
                if (parameter.geom.mesh.dim == 2): 
                    if (Kxy != None): r[i].permeability.K.append(float(Kxy.text)) 
                    if (Kyy != None): r[i].permeability.K.append(float(Kyy.text))
                    if ((Kxy == None) | (Kyy == None)):
                        print "Error(4.1):"
                        print "Missing permeability information for 2-D problem"
                        print "This may indicate that you forgot the Kxy or Kyy data in the .xml file"
                        print " "
                        sys.exit(1)
                # For 3-D:
                # K = | Kxx Kxy Kxz | Kxy = Kyx e Kyz = Kzy
                #        | Kyx Kyy Kyz | 
                #        | Kzx Kzy Kzz |
                elif(parameter.geom.mesh.dim == 3): # 3-D
                    if (Kxy != None): r[i].permeability.K.append(float(Kxy.text)) 
                    if (Kxz != None): r[i].permeability.K.append(float(Kxz.text))   
                    if (Kyy != None): r[i].permeability.K.append(float(Kyy.text)) 
                    if (Kyz != None): r[i].permeability.K.append(float(Kyz.text))  
                    if (Kzz != None): r[i].permeability.K.append(float(Kzz.text))
                    if ((Kxy == None) | (Kxz == None) | (Kyy == None) | (Kyz == None) | (Kzz ==None)):
                        print "Error(4.2):"
                        print "Missing permeability information for 3-D problem"
                        print "This may indicate that you forgot the Kxy or Kxz or Kyy \
                    or Kyz or Kzz data in the .xml file "
                        print " "
                        sys.exit(1)
            elif(r[i].permeability.type == "per-element-list"):
                r[i].permeability.filename = etRockTypes[i].find('permeability').findtext('filename')
            
            parameter.phys.rock.rocks.insert(i, r[i])
            # Writing in .log file for check purposes
            log.write("Rock ID: " + str(parameter.phys.rock.rocks[i].id) +"\n")
            log.write("Rock Porosity Compressibility?: " + str(parameter.phys.rock.rocks[i].porosity.compressible) +"\n")
            log.write("Rock Porosity Value: " + str(parameter.phys.rock.rocks[i].porosity.value) +"\n")
            log.write("Rock Heat Coefficient Value: " + str(parameter.phys.rock.rocks[i].rockHeatCoeff.value) +"\n")
            log.write("Rock Permeability Values: " + str(parameter.phys.rock.rocks[i].permeability.K) +"\n")
        
        # ## Rock-Fluid properties - phys.rockFluid
        kr = rockFluidProp.find('relative-permeability')
        parameter.phys.rockFluid.relativePerm.type = kr.attrib["model"]
        krw = kr.find("krw")
        kro = kr.find("kro")
        krg = kr.find("krg")
        
        # Read data for Corey model (Attention: other models(Stone's I, II) still pending)
        if (parameter.phys.rockFluid.relativePerm.type == "corey"):
            # Checking if not compulsory data is given in input file (as only krw is alway mandatory)
            # Important also to note that the 'krn_end' property doesnt exist for the oil phase (s. Aziz & Settari)
            parameter.phys.rockFluid.relativePerm.krw.krEnd = float(krw.findtext('krw_end'))
            if (krg != None): parameter.phys.rockFluid.relativePerm.krg.krEnd = float(krg.findtext('krg_end'))
            
            parameter.phys.rockFluid.relativePerm.krw.Sr = float(krw.findtext('Swc'))
            if (kro != None): parameter.phys.rockFluid.relativePerm.kro.Sr = float(kro.findtext('Sor'))
            if (krg != None): parameter.phys.rockFluid.relativePerm.krg.Sr = float(krg.findtext('Sgc'))
            
            parameter.phys.rockFluid.relativePerm.krw.n = float(krw.findtext('nw'))
            if (kro != None): parameter.phys.rockFluid.relativePerm.kro.n = float(kro.findtext('no'))
            if (krg != None): parameter.phys.rockFluid.relativePerm.krg.n = float(krg.findtext('ng'))
            
            # Writing in .log file for check purposes
            log.write("Rock Relative Permeability Model : " + str(parameter.phys.rockFluid.relativePerm.type) +"\n")
            log.write("Rock Water Relative Permeability (krwEnd | Swc | nw) : " + \
            str(parameter.phys.rockFluid.relativePerm.krw.krEnd) + " " + str(parameter.phys.rockFluid.relativePerm.krw.Sr) + \
            " " + str(parameter.phys.rockFluid.relativePerm.krw.n) + "\n")
            
            log.write("Rock Oil Relative Permeability (Soc | no) : " + \
            str(parameter.phys.rockFluid.relativePerm.kro.Sr) + " " + str(parameter.phys.rockFluid.relativePerm.kro.n) + "\n")
            
            log.write("Rock Gas Relative Permeability (krgEnd | Sgc | ng) : " + \
            str(parameter.phys.rockFluid.relativePerm.krg.krEnd) + " " + str(parameter.phys.rockFluid.relativePerm.krg.Sr) + \
            " " + str(parameter.phys.rockFluid.relativePerm.krg.n) + "\n")
        
        # # Numerical parameters - num
        pressSolv = num.find('pressure-solver')
        saturSolv = num.find('saturation-solver')
        
        # ## Pressure Solver parameters - num.pressSolv
        parameter.num.pressSolv.formulation = pressSolv.attrib["formulation"]
        parameter.num.pressSolv.type = pressSolv.attrib["type"]
        parameter.num.pressSolv.tolerance = float(pressSolv.findtext('tolerance'))
        parameter.num.pressSolv.maxNumSteps = int(pressSolv.findtext('max-number-steps'))
        parameter.num.pressSolv.preConditioning.type = pressSolv.find('pre-conditioning').attrib["type"]
        # Writing in .log file for check purposes
        log.write("Pressure Solver Formulation: " + str(parameter.num.pressSolv.formulation) +"\n")
        log.write("Pressure Solver Type: " + str(parameter.num.pressSolv.type) +"\n")
        log.write("Pressure Solver Tolerance: " + str(parameter.num.pressSolv.tolerance) +"\n")
        log.write("Pressure Solver Max. Num. Steps: " + str(parameter.num.pressSolv.maxNumSteps) +"\n")
        log.write("Pressure Solver Pre-Conditioner: " + str(parameter.num.pressSolv.preConditioning.type) +"\n")
        
        # Read data for Algebraic Multigrid(AMG) pre-conditioner (other types (e.g, LU) still pending)
        if (parameter.num.pressSolv.preConditioning.type == "amg"):
            # Checking if not compulsory data is given in input file (as nCoarse and nRelaxIter are mandatory
            # only when using algebraic multigrid)
            nCoarse = pressSolv.find('pre-conditioning').find('number-coarse-levels')
            if (nCoarse != None):
                parameter.num.pressSolv.preConditioning.numCoarseLevel = int(nCoarse.text)
            else: # Default value in case this data is not given
                parameter.num.pressSolv.preConditioning.numCoarseLevel = 2
            
            nRelaxIter = pressSolv.find('pre-conditioning').find('number-relax-iterations')
            if (nRelaxIter != None):
                parameter.num.pressSolv.preConditioning.numRelaxIter = int(nRelaxIter.text)
            else: # Default value in case this data is not given
                parameter.num.pressSolv.preConditioning.numRelaxIter = 2
            
            # Writing in .log file for check purposes
            log.write("AMG Number of Coarse Levels: " + str(parameter.num.pressSolv.preConditioning.numCoarseLevel) +"\n")
            log.write("AMG Number of Relaxation Iterations: " + str(parameter.num.pressSolv.preConditioning.numRelaxIter) +"\n")
        
        # ## Saturation Solver parameters - num.saturSolv
        parameter.num.saturSolv.totalTime = float(saturSolv.findtext('total-time-analysis'))
        parameter.num.saturSolv.courant = float(saturSolv.findtext('courant'))
        parameter.num.saturSolv.limiterType = saturSolv.find('limiter').attrib["type"]
        
        # Writing in .log file for check purposes
        log.write("Total Time of Analysis: " + str(parameter.num.saturSolv.totalTime) +"\n")
        log.write("Saturation Solver Courant Number: " + str(parameter.num.saturSolv.courant) +"\n")
        log.write("Saturation Solver Limiter Type: " + str(parameter.num.saturSolv.limiterType) +"\n")
        
        # Closing the check.log file
        log.close()
#############################################################
