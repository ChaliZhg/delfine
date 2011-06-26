# Copyright (C) 2006 Anders Logg
# Licensed under the GNU LGPL Version 2.1
#
# Modified by Garth N. Wells (gmsh function)
# Modified by Alexander H. Jarosch (gmsh fix)
# Modified by Angelo Simone (Gmsh and Medit fix)
# Modified by Andy R. Terrel (gmsh fix and triangle function)
# Modified by Magnus Vikstrom (metis and scotch function)
# Modified by Bartosz Sawicki (diffpack function)
# Modified by Gideon Simpson (Exodus II function)
# Modified by Arve Knudsen (make into module, abaqus support)
# Modified by Kent-Andre Mardal (Star-CD function)
# Modified by Nuno Lopes (fix for emc2 mesh format (medit version 0))
#
# This is the standard meshconvert module
# with the support for physical regions in
# from gmsh, and it invokes just the 
# gmsh2xml function of Neilen Marais
# This was done in order to not mess with
# the standard function and still use the
# functionality of Neilen's Patch
# Modified by Bruno Luna
""" Module for converting various mesh formats.
"""

import getopt
import sys
from dolfin_utils.commands import getoutput
import re
import warnings
import os.path

def format_from_suffix(suffix):
    "Return format for given suffix"
    if suffix == "xml":
        return "xml"
    elif suffix == "mesh":
        return "mesh"
    elif suffix == "gmsh":
        return "gmsh"
    elif suffix == "msh":
        return "gmsh"
    elif suffix == "gra":
        return "metis"
    elif suffix == "grf":
        return "scotch"
    elif suffix == "grid":
        return "diffpack"
    elif suffix == "inp":
        return "abaqus"
    elif suffix == "ncdf":
        return "NetCDF"
    elif suffix =="exo":
        return "ExodusII"
    elif suffix =="e":
        return "ExodusII"
    elif suffix == "vrt" or suffix == "cel":
        return "StarCD"
    else:
        _error("Sorry, unknown suffix %s." % suffix)

def gmsh2xml(ifilename, handler):
    """Convert between .gmsh v2.0 format (http://www.geuz.org/gmsh/) and .xml,
    parser implemented as a state machine:

        0 = read 'MeshFormat'
        1 = read  mesh format data
        2 = read 'EndMeshFormat'
        3 = read 'Nodes'
        4 = read  number of vertices
        5 = read  vertices
        6 = read 'EndNodes'
        7 = read 'Elements'
        8 = read  number of cells
        9 = read  cells
        10 = done

    """

    print "Converting from Gmsh format (.msh, .gmsh) to DOLFIN XML format"

    # Open files
    ifile = open(ifilename, "r")

    # Scan file for cell type
    cell_type = None
    dim = 0
    # Delfine-convert:------------------------------------------------
    pointTag = False
    # end-delfine-convert---------------------------------------------
    line = ifile.readline()
    while line:

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]
        
        # Delfine-convert:------------------------------------------------
        if line.find("$Nodes") == 0:
            line = ifile.readline()
            num_nodes = int(line)
        
        # end-delfine-convert---------------------------------------------
        # Read dimension
        if line.find("$Elements") == 0:

            line = ifile.readline()
            num_cells  = int(line)
            num_cells_counted = 0
            if num_cells == 0:
                _error("No cells found in gmsh file.")
            line = ifile.readline()

            # Now iterate through elements to find largest dimension.  Gmsh
            # format might include elements of lower dimensions in the element list.
            # We also need to count number of elements of correct dimensions.
            # Also determine which vertices are not used.
            dim_2_count = 0
            dim_3_count = 0
            # Delfine-convert:------------------------------------------------
            # Array used to store gmsh tags for 1D (type 1/point) elements
            dim_1_count = 0
            tags_1 = []
            vertices_1_used = []
            # end-delfine-convert---------------------------------------------
            # Array used to store gmsh tags for 2D (type 2/triangular) elements
            tags_2 = []
            vertices_2_used = []
            # Array used to store gmsh tags for 3D (type 4/tet) elements
            tags_3 = []
            vertices_3_used = []
            while line.find("$EndElements") == -1:
                element = line.split()
                elem_type = int(element[1])
                num_tags = int(element[2])
                # Delfine-convert:------------------------------------------------
                if elem_type == 15:
                    if pointTag == False:                    
                        pointTag = True
                    node_num_list = [int(node) for node in element[3 + num_tags:]]
                    vertices_1_used.extend(node_num_list)
                    if num_tags > 0:
                        tags_1.append(tuple(int(tag) for tag in element[3:3+num_tags]))
                    dim_1_count += 1
                # end-delfine-convert---------------------------------------------
                if elem_type == 2:
                    if dim < 2:
                        cell_type = "triangle"
                        dim = 2
                    node_num_list = [int(node) for node in element[3 + num_tags:]]
                    vertices_2_used.extend(node_num_list)
                    if num_tags > 0:
                        tags_2.append(tuple(int(tag) for tag in element[3:3+num_tags]))
                    dim_2_count += 1
                elif elem_type == 4:
                    if dim < 3:
                        cell_type = "tetrahedron"
                        dim = 3
                        vertices_2_used = None
                    node_num_list = [int(node) for node in element[3 + num_tags:]]
                    vertices_3_used.extend(node_num_list)
                    if num_tags > 0:
                        tags_3.append(tuple(int(tag) for tag in element[3:3+num_tags]))
                    dim_3_count += 1
                line = ifile.readline()
        else:
            # Read next line
            line = ifile.readline()

    # Check that we got the cell type and set num_cells_counted
    if cell_type == None:
        _error("Unable to find cell type.")
    if dim == 3:
        num_cells_counted = dim_3_count
        vertex_set = set(vertices_3_used)
        vertices_3_used = None
    elif dim == 2:
        num_cells_counted = dim_2_count
        vertex_set = set(vertices_2_used)
        vertices_2_used = None

    vertex_dict = {}
    for n,v in enumerate(vertex_set):
        vertex_dict[v] = n

    # Step to beginning of file
    ifile.seek(0)
    
    # Set mesh type
    handler.set_mesh_type(cell_type, dim)

    # Initialise node list (gmsh does not export all vertexes in order)
    nodelist = {}

    # Current state
    state = 0

    # Write data
    num_vertices_read = 0
    num_cells_read = 0

    while state != 10:

        # Read next line
        line = ifile.readline()
        if not line: break

        # Skip comments
        if line[0] == '#':
            continue

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        if state == 0:
            if line == "$MeshFormat":
                state = 1
        elif state == 1:
            (version, file_type, data_size) = line.split()
            state = 2
        elif state == 2:
            if line == "$EndMeshFormat":
                state = 3
        elif state == 3:
            if line == "$Nodes":
                state = 4
        elif state == 4:
            num_vertices = len(vertex_dict)
            handler.start_vertices(num_vertices)
            state = 5
        elif state == 5:
            (node_no, x, y, z) = line.split()
            if vertex_dict.has_key(int(node_no)):
                node_no = vertex_dict[int(node_no)]
            else:
                continue
            nodelist[int(node_no)] = num_vertices_read
            handler.add_vertex(num_vertices_read, [x, y, z])
            num_vertices_read +=1

            if num_vertices == num_vertices_read:
                handler.end_vertices()
                state = 6
        elif state == 6:
            if line == "$EndNodes":
                state = 7
        elif state == 7:
            if line == "$Elements":
                state = 8
        elif state == 8:
            handler.start_cells(num_cells_counted)
            state = 9
        elif state == 9:
            element = line.split()
            elem_type = int(element[1])
            num_tags  = int(element[2])
            if elem_type == 2 and dim == 2:
                node_num_list = [vertex_dict[int(node)] for node in element[3 + num_tags:]]
                for node in node_num_list:
                    if not node in nodelist:
                        _error("Vertex %d of triangle %d not previously defined." %
                              (node, num_cells_read))
                cell_nodes = [nodelist[n] for n in node_num_list]
                handler.add_cell(num_cells_read, cell_nodes)
                num_cells_read +=1
            elif elem_type == 4 and dim == 3:
                node_num_list = [vertex_dict[int(node)] for node in element[3 + num_tags:]]
                for node in node_num_list:
                    if not node in nodelist:
                        _error("Vertex %d of tetrahedron %d not previously defined." %
                              (node, num_cells_read))
                cell_nodes = [nodelist[n] for n in node_num_list]
                handler.add_cell(num_cells_read, cell_nodes)
                num_cells_read +=1

            if num_cells_counted == num_cells_read:
                handler.end_cells()
                state = 10
        elif state == 10:
            break

    # Write mesh function based on the Physical Regions defined by
    # gmsh, but only if they are not all zero. All zero physical
    # regions indicate that no physical regions were defined.
    if dim == 2:
        tags = tags_2
    elif dim == 3:
        tags = tags_3
    else:
        _error("Gmsh tags not supported for dimension %i. Probably a bug" % dim)

    # Delfine-convert---------------------------------------------------
    # Initialize general meshfunction header in case there are any tags
    if not all(tag == 0 for tag in tags):
        handler.start_meshfunctionfile()
    elif not all(tagW == 0 for tagW in tags_1):
        handler.start_meshfunctionfile()
    # end - delfine-convert---------------------------------------------

    # Delfine-convert---------------------------------------------------
    # This snippet is executed only if well(physical) points were found
    if (dim_1_count > 0):
        physical_points = tuple(tag[0] for tag in tags_1)
        point_index = tuple(tag[1] for tag in tags_1)
        if not all(tag == 0  for tag in tags_1):
            handler.start_meshfunction("well_indicators", 0, num_nodes)
            wellType = [0]*num_nodes
            for i in xrange(num_nodes):
                for j in xrange(len(point_index)):
                    if (i == (point_index[j] - 1)):
                        wellType[i] = physical_points[j]
                handler.add_entity_meshfunction((i), wellType[i])
            handler.end_meshfunction()
    # end - delfine-convert---------------------------------------------

    physical_regions = tuple(tag[0] for tag in tags)
    if not all(tag == 0 for tag in tags):
        handler.start_meshfunction("material_indicators", dim, num_cells_counted)
        for i, physical_region in enumerate(physical_regions):
            handler.add_entity_meshfunction(i, physical_region)
        handler.end_meshfunction()
       
    # Delfine-convert---------------------------------------------------
    # Finalize general meshfunction footer in case there were any tags
    if not all(tag == 0 for tag in tags):
        handler.end_meshfunctionfile()
    elif not all(tagW == 0 for tagW in tags_1):
        handler.end_meshfunctionfile()
    # end - delfine-convert---------------------------------------------

    # Check that we got all data
    if state == 10:
        print "Conversion done"
    else:
       _error("Missing data, unable to convert \n\ Did you use version 2.0 of the gmsh file format?")

    # Close files
    ifile.close()

def write_header_meshfunction(ofile, name, dimensions, size):
    # Delfine-convert: Alteration to ensure that all the information
    # is saved in just one file
    header = """        <data_entry name="%s">
          <meshfunction type="uint" dim="%d" size="%d">\n""" % (name, dimensions, size)
    ofile.write(header)

def write_entity_meshfunction(ofile, index, value):
    ofile.write("""            <entity index=\"%d\" value=\"%d\"/>
""" % (index, value))

def write_footer_meshfunction(ofile):
    # Delfine-convert: Alteration to ensure that all the information
    # is saved in just one file
    ofile.write("""          </meshfunction>
        </data_entry>\n""")

# Delfine-convert---------------------------------------------------
# General header and footer for meshfunction file
def write_header_meshfunctionfile(ofile):
    ofile.write("""      <data>\n""")

def write_footer_meshfunctionfile(ofile):
    ofile.write("""      </data>\n""")
# end - delfine-convert---------------------------------------------

class ParseError(Exception):
    """ Error encountered in source file.
    """

class DataHandler(object):
    """ Baseclass for handlers of mesh data.

    The actual handling of mesh data encountered in the source file is
    delegated to a polymorfic object. Typically, the delegate will write the
    data to XML.
    @ivar _state: the state which the handler is in, one of State_*.
    @ivar _cell_type: cell type in mesh. One of CellType_*.
    @ivar _dim: mesh dimensions.
    """
    State_Invalid, State_Init, State_Vertices, State_Cells, State_MeshFunction = range(5)
    CellType_Tetrahedron, CellType_Triangle = range(2)

    def __init__(self):
        self._state = self.State_Invalid

    def set_mesh_type(self, cell_type, dim):
        assert self._state == self.State_Invalid
        self._state = self.State_Init
        if cell_type == "tetrahedron":
            self._cell_type = self.CellType_Tetrahedron
        elif cell_type == "triangle":
            self._cell_type = self.CellType_Triangle
        self._dim = dim

    def start_vertices(self, num_vertices):
        assert self._state == self.State_Init
        self._state = self.State_Vertices

    def add_vertex(self, vertex, coords):
        assert self._state == self.State_Vertices

    def end_vertices(self):
        assert self._state == self.State_Vertices
        self._state = self.State_Init

    def start_cells(self, num_cells):
        assert self._state == self.State_Init
        self._state = self.State_Cells

    def add_cell(self, cell, nodes):
        assert self._state == self.State_Cells

    def end_cells(self):
        assert self._state == self.State_Cells
        self._state = self.State_Init

    def start_meshfunction(self, name, dim, size):
        assert self._state == self.State_Init
        self._state = self.State_MeshFunction

    def add_entity_meshfunction(self, index, value):
        assert self._state == self.State_MeshFunction

    def end_meshfunction(self):
        assert self._state == self.State_MeshFunction
        self._state = self.State_Init
     
    # Delfine-convert:-----------------------------------------------------------
    def start_meshfunctionfile(self):
        assert self._state == self.State_Init
        self._state = self.State_Init
    def end_meshfunctionfile(self):
        assert self._state == self.State_Init
        self._state = self.State_Init
    # end-delfine-convert:-----------------------------------------------------------
    
    def warn(self, msg):
        """ Issue warning during parse.
        """
        warnings.warn(msg)

    def error(self, msg):
        """ Raise error during parse.

        This method is expected to raise ParseError.
        """
        raise ParseError(msg)

    def close(self):
        self._state = self.State_Invalid

class XmlHandler(DataHandler):
    """ Data handler class which writes to Dolfin XML.
    """
    def __init__(self, ofilename):
        DataHandler.__init__(self)
        self._ofilename = ofilename
        self.__ofile = file(ofilename, "wb")
        self.__ofile_meshfunc = None

    def set_mesh_type(self, cell_type, dim):
        DataHandler.set_mesh_type(self, cell_type, dim)
        write_header_mesh(self.__ofile, cell_type, dim)

    def start_vertices(self, num_vertices):
        DataHandler.start_vertices(self, num_vertices)
        write_header_vertices(self.__ofile, num_vertices)

    def add_vertex(self, vertex, coords):
        DataHandler.add_vertex(self, vertex, coords)
        write_vertex(self.__ofile, vertex, *coords)

    def end_vertices(self):
        DataHandler.end_vertices(self)
        write_footer_vertices(self.__ofile)

    def start_cells(self, num_cells):
        DataHandler.start_cells(self, num_cells)
        write_header_cells(self.__ofile, num_cells)

    def add_cell(self, cell, nodes):
        DataHandler.add_cell(self, cell, nodes)
        if self._cell_type == self.CellType_Tetrahedron:
            func = write_cell_tetrahedron
        if self._cell_type == self.CellType_Triangle:
            func = write_cell_triangle
        func(self.__ofile, cell, *nodes)

    def end_cells(self):
        DataHandler.end_cells(self)
        write_footer_cells(self.__ofile)

    def start_meshfunction(self, name, dim, size):
        DataHandler.start_meshfunction(self, name, dim, size)
        #fname = os.path.splitext(self.__ofile.name)[0]
        #self.__ofile_meshfunc = file("%s_%s.xml" % (fname, name), "wb")
        #write_header_meshfunction(self.__ofile_meshfunc, dim, size)
        # Delfine-convert: Alteration to ensure that all the information
        # is saved in just one file
        write_header_meshfunction(self.__ofile, name, dim, size)

    def add_entity_meshfunction(self, index, value):
        DataHandler.add_entity_meshfunction(self, index, value)
        #write_entity_meshfunction(self.__ofile_meshfunc, index, value)
        # Delfine-convert: Alteration to ensure that all the information
        # is saved in just one file
        write_entity_meshfunction(self.__ofile, index, value)

    def end_meshfunction(self):
        DataHandler.end_meshfunction(self)
        #write_footer_meshfunction(self.__ofile_meshfunc)
        #self.__ofile_meshfunc.close()
        #self.__ofile_meshfunc = None
        # Delfine-convert: Alteration to ensure that all the information
        # is saved in just one file
        write_footer_meshfunction(self.__ofile)
        self.__ofile_meshfunc = None
    
    # Delfine-convert:-----------------------------------------------------------
    def start_meshfunctionfile(self):
        DataHandler.start_meshfunctionfile(self)
        write_header_meshfunctionfile(self.__ofile)
    def end_meshfunctionfile(self):
        DataHandler.end_meshfunctionfile(self)
        write_footer_meshfunctionfile(self.__ofile)
    # end-delfine-convert:-----------------------------------------------------------
    
    def close(self):
        DataHandler.close(self)
        if self.__ofile.closed:
            return
        write_footer_mesh(self.__ofile)
        self.__ofile.close()
        # Delfine-convert: Alteration to ensure that all the information
        # is saved in just one file
        #if self.__ofile_meshfunc is not None:
        #    self.__ofile_meshfunc.close()

# Write mesh header
def write_header_mesh(ofile, cell_type, dim):
    ofile.write("""\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>

<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">
  <mesh celltype="%s" dim="%d">
""" % (cell_type, dim))

# Write graph header
def write_header_graph(ofile, graph_type):
    ofile.write("""\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>

<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">
  <graph type="%s">
""" % (graph_type))

# Write mesh footer
def write_footer_mesh(ofile):
    ofile.write("""\
  </mesh>
</dolfin>
""")

# Write graph footer
def write_footer_graph(ofile):
    ofile.write("""\
  </graph>
</dolfin>
""")

def write_header_vertices(ofile, num_vertices):
    "Write vertices header"
    print "Expecting %d vertices" % num_vertices
    ofile.write("    <vertices size=\"%d\">\n" % num_vertices)

def write_footer_vertices(ofile):
    "Write vertices footer"
    ofile.write("    </vertices>\n")
    print "Found all vertices"

def write_header_edges(ofile, num_edges):
    "Write edges header"
    print "Expecting %d edges" % num_edges
    ofile.write("    <edges size=\"%d\">\n" % num_edges)

def write_footer_edges(ofile):
    "Write edges footer"
    ofile.write("    </edges>\n")
    print "Found all edges"

def write_vertex(ofile, vertex, x, y, z):
    "Write vertex"
    ofile.write("      <vertex index=\"%d\" x=\"%s\" y=\"%s\" z=\"%s\"/>\n" % \
        (vertex, x, y, z))

def write_graph_vertex(ofile, vertex, num_edges, weight = 1):
    "Write graph vertex"
    ofile.write("      <vertex index=\"%d\" num_edges=\"%d\" weight=\"%d\"/>\n" % \
        (vertex, num_edges, weight))

def write_graph_edge(ofile, v1, v2, weight = 1):
	 "Write graph edge"
	 ofile.write("      <edge v1=\"%d\" v2=\"%d\" weight=\"%d\"/>\n" % \
        (v1, v2, weight))

def write_header_cells(ofile, num_cells):
    "Write cells header"
    ofile.write("    <cells size=\"%d\">\n" % num_cells)
    print "Expecting %d cells" % num_cells

def write_footer_cells(ofile):
    "Write cells footer"
    ofile.write("    </cells>\n")
    print "Found all cells"

def write_cell_triangle(ofile, cell, n0, n1, n2):
    "Write cell (triangle)"
    ofile.write("      <triangle index=\"%d\" v0=\"%d\" v1=\"%d\" v2=\"%d\"/>\n" % \
        (cell, n0, n1, n2))

def write_cell_tetrahedron(ofile, cell, n0, n1, n2, n3):
    "Write cell (tetrahedron)"
    ofile.write("      <tetrahedron index=\"%d\" v0=\"%d\" v1=\"%d\" v2=\"%d\" v3=\"%d\"/>\n" % \
        (cell, n0, n1, n2, n3))

def _error(message):
    "Write an error message"
    for line in message.split("\n"):
        print "*** %s" % line
    sys.exit(2)

def convert2xml(ifilename, ofilename, iformat=None):
    """ Convert a file to the DOLFIN XML format.
    """
    convert(ifilename, XmlHandler(ofilename), iformat=iformat)

def convert(ifilename, handler, iformat=None):
    """ Convert a file using a provided data handler.

    Note that handler.close is called when this function finishes.
    @param ifilename: Name of input file.
    @param handler: The data handler (instance of L{DataHandler}).
    @param iformat: Format of input file.
    """
    if iformat is None:
        iformat = format_from_suffix(os.path.splitext(ifilename)[1][1:])
    # XXX: Backwards-compat
    if hasattr(handler, "_ofilename"):
        ofilename = handler._ofilename
    # Choose conversion
    #if iformat == "mesh":
    #    # Convert from mesh to xml format
    #    mesh2xml(ifilename, ofilename)
    #elif iformat == "gmsh":
    #########################################################################
    # Delfine-convert only!! All other file formats set as unable.
    # In order to use then with delfine, it is necessary to run the
    # standard dolfin-convert and then use the resultant .xml file.
    if iformat == "gmsh":
        # Convert from gmsh to xml format
        gmsh2xml(ifilename, handler)
    #########################################################################
    #elif iformat == "Triangle":
    #    # Convert from Triangle to xml format
    #    triangle2xml(ifilename, ofilename)
    #elif iformat == "xml-old":
    #   # Convert from old to new xml format
    #    xml_old2xml(ifilename, ofilename)
    #elif iformat == "metis":
    #    # Convert from metis graph to dolfin graph xml format
    #    metis_graph2graph_xml(ifilename, ofilename)
    #elif iformat == "scotch":
    #    # Convert from scotch graph to dolfin graph xml format
    #    scotch_graph2graph_xml(ifilename, ofilename)
    #elif iformat == "diffpack":
    #    # Convert from Diffpack tetrahedral grid format to xml format
    #    diffpack2xml(ifilename, ofilename)
    #elif iformat == "abaqus":
    #    # Convert from abaqus to xml format
    #    _abaqus(ifilename, handler)
    #elif iformat == "NetCDF":
    #    # Convert from NetCDF generated from ExodusII format to xml format
    #    netcdf2xml(ifilename, ofilename)
    #elif iformat =="ExodusII":
    #    # Convert from ExodusII format to xml format via NetCDF
    #    exodus2xml(ifilename, ofilename)
    #elif iformat == "StarCD":
    #    # Convert from Star-CD tetrahedral grid format to xml format
    #    starcd2xml(ifilename, ofilename)
    else:
        _error("Sorry, delfine-convert cannot convert between %s and DOLFIN xml file formats. \
Try to use the standard dolfin-convert and then run delfine with the resultant .xml file" % iformat)

    # XXX: handler.close messes things for other input formats than abaqus or gmsh
    if iformat in ("abaqus", "gmsh"):
        handler.close()
