/*!
 * @page o3d_xml_manual The OPTIMET-3D XML Input File Manual
 *
 * @section o3d_xml_introduction Introduction
 * OPTIMET-3D uses an eXtended Markup Language (XML) file to setup
 * its simulation. The .xml file describes the entire simulation and its
 * name will be appended to all output files. When launching a simulation,
 * the .xml file name needs to be passed as a parameter to the Optimet3D
 * executable WITHOUT the .xml extension.
 *
 * The .xml file consists of several blocks describing different aspects
 * of the simulation. The order of the blocks is irrelevant, but this manual describes
 * a typical setup which goes through all stages of the simulation.
 *
 * @section o3d_xml_simulation The Simulation Block
 * The simulation block describes the overall properties of the simulation.
 * The only property currently supported is the maximum number of harmonics (nMax).
 * The block syntax (for 20 harmonics) is:
 * @code{.xml}
 * <simulation>
 *  <harmonics nmax="20" />
 * </simulation>
 * @endcode
 *
 * @section o3d_xml_geometry The Geometry Block
 * The geometry block defines the geometry of the simulation. The block supports
 * two types of geometries, an object by object geometry and a structured geometry.
 * Once the geometry is created, OPTIMET-3D will validate it to ensure that no
 * two objects intersect. At the present time, only vacuum backgrounds are
 * supported with variable permitivity and permeability backgrounds to be added.
 *
 * @subsection o3d_xml_sub_object The Object Sub-block
 * The object sub-block defines a scattered in terms of position, geometry and
 * electromagnetic properties. Several object sub-blocks can be added to the
 * geometry one after the other. The actual order they will be stored in the
 * Geometry object is given by the order they appear in the .xml file. Currently,
 * only the sphere object is supported.
 *
 * @subsubsection o3d_xml_subsub_sphere The Sphere Object
 * A sphere object is defined by its coordinates (Cartesian or Spherical), its radius and
 * its electromagnetic properties. A typical sphere object, using static electromagnetic
 * properties and Cartesian coordinates, is defined as (all lengths in nanometers):
 *
 * @code{.xml}
 * <object type="sphere">
 *  <cartesian x="100.0" y="50.0" z="-200.0" />
 *  <properties radius="300.0" />
 *  <epsilon type="relative" value="13" />
 *  <mu type="relative" value="1" />
 * </object>
 * @endcode
 *
 * Note that the leading .0 in double values is optional. However, if an integer value contains a .0,
 * the code will trigger an error message before rounding, so for clarity it is preferred.
 *
 * The object sub-block also supports Sellmeier type equations for the relative permitivity. The following
 * block defines a Si sphere using the Sellmeier model and Spherical coordinates in nanometers and radians
 * (note that Cartesian work as defined above for this case):
 *
 * @code{.xml}
 * <object type="sphere">
 *  <spherical rrr="100.0" the="1.8" phi="3.14" />
 *  <properties radius="200.0" />
 *  <epsilon type="sellmeier"/>
 *    <parameters B1="10.6684293" C1="0.301516485" B2="0.003043475" C2="1.13475115" B3="1.54133408" C3="1104.0" />
 *  </epsilon>
 *  <mu type="relative" value="1" />
 * </object>
 * @endcode
 *
 * @subsection o3d_xml_sub_structure The Structure Sub-Block
 * The structure block defines a pre-made structure. It allows for quick writing of simulation files
 * with pre-determined geometric structures, without having to manually add each scatterer. A structure
 * consists of an object sub-block (see @ref o3d_xml_sub_object) and various
 * geometrical parameters. The structure will be formed of objects identical to the one defined in the
 * object sub-block.
 *
 * OPTIMET-3D currently supports one type of structure, the spiral. A spiral is defined by its number of
 * arms, number of scatterers (points) on each arm and the separation distance between the objects on an arm, which
 * is kept constant (the distance is defined as center to center). The block below defines a typical three arm
 * spiral with five points (including the central one which is shared) separated by 100 nanometers (surface to surface):
 *
 * @code {.xml}
 * <structure type="spiral" arms="3">
 *  <object type="sphere">
 *    <properties radius="300.0" />
 *    <epsilon type="relative" value="13" />
 *    <mu type="relative" value="1" />
 *  </object>
 *  <properties points="5" distance="700.0" />
 * </structure>
 * @endcode
 *
 * Note that the Cartesian and Spherical coordinates of the object are not given as these are calculated by the
 * structure builder. If any values are given, they will be ignored. At the current time, if a structure sub-block
 * is present in the geometry block, ALL other sub-blocks (object and/or any other structure) will be ignored.
 *
 * @section o3d_xml_source The Source Block
 *
 * The source block defines the excitation (or incoming wave) that acts upon the geometry. A source is defined
 * by its type (planewaves currently supported, Laguerre-Gaussian beams to be included in the near future), wavelength,
 * incidence vector (or angle) and its polarization. A typical block used to define a planewave is:
 *
 * @code{.xml}
 * <source type="planewave">
 *  <wavelength value="1460" />
 *  <propagation theta="90" phi="90" />
 *  <polarization Etheta.real="0.70710678118" Etheta.imag="0.0" Ephi.real="0.0" Ephi.imag="0.70710678118"/>
 * </source>
 * @endcode
 *
 * Here, the wavelength is given in nanometers, the propagation is defined with respect to the center of the cluster
 * and given in terms of polar (theta) and azimuthal (phi) angles and the polarization is defined in a plane perpendicular
 * to the propagation vector. Both real and imaginary components are accepted for the Etheta and Ephi componens of
 * the incoming plane wave. Note that in the case above, this will lead to a circular polarization. The intensity of the
 * incoming plane wave is directly defined in the polarization section by the value of the components (in V/m^2).
 *
 * @section o3d_xml_output The Output Block
 *
 * The output block is used to define the type of output that the user desires. Depending on the request, several output
 * files may be created. At the moment, a single output is possible, although future versions will allow for more than
 * one output block. There are several types of output blocks, explained in what follows.
 *
 * @subsection o3d_xml_sub_fields Field Profile Output
 *
 * The field profile output creates a profile of the electric and magnetic field (for all three vectorial components)
 * on a grid defined by the user. Based on the information required, the user can opt for a full field profile or
 * a profile of just one of the harmonic components of the Fourier-Bessel expansion.
 *
 * @subsubsection o3d_xml_subsub_fullfields Full Field Profile
 *
 * A full field profile will take into account ALL harmonic components of the Fourier-Bessel expansion. A full field profile
 * is requested using the following block:
 *
 * @code{.xml}
 * <output type="field">
 *  <grid type="cartesian">
 *    <x min="-1000.0" max="1000.0" steps="101" />
 *    <y min="-1000.0" max="1000.0" steps="101" />
 *    <z min="-1000.0" max="1000.0" steps="101" />
 *  </grid>
 * </output>
 * @endcode
 *
 * Here, the grid is defined of type "cartesian" (the only available type at the moment) and grid min and max values for
 * each coordinate are in nanometers and relative to the center of the cluster.
 * Note that for any simulation, at least 2 (TWO) steps are required for each coordinate.
 *
 * @subsubsection o3d_xml_subsub_singlefields Single Harmonic Field Profile
 *
 * The single harmonic field profile is defined identical to the full field profile, but a special sub-block is added
 * to request only one harmonic to be used. Also, it is possible to output only the TE or TM part of the harmonic. The
 * block bellow requests a single harmonic (with n=1 and m=0) for the TM part:
 *
 * @code{.xml}
 * <output type="field">
 *  <grid type="cartesian">
 *    <x min="-1000.0" max="1000.0" steps="101" />
 *    <y min="-1000.0" max="1000.0" steps="101" />
 *    <z min="-1000.0" max="1000.0" steps="101" />
 *  </grid>
 *  <singlemode n="1" m="0" component="TM" />
 * </output>
 * @endcode
 *
 * It is also possible to let OPTIMET-3D calculate the dominant harmonic for each simulation. The code below
 * outputs the dominant harmonic profile for both TE and TM (the latter part is done simply by not adding a component="" property
 * to the singlemode sub-block).
 *
 * @code{.xml}
 * <output type="field">
 *  <grid type="cartesian">
 *    <x min="-1000.0" max="1000.0" steps="101" />
 *    <y min="-1000.0" max="1000.0" steps="101" />
 *    <z min="-1000.0" max="1000.0" steps="101" />
 *  </grid>
 *  <singlemode dominant="auto" />
 * </output>
 * @endcode
 *
 * All output will be directed to an HDF5 file named after the project name (identical to the name of .xml file, without the extension).
 *
 * @subsection o3d_xml_response Spectral Response Output
 *
 * The spectral response output is used to perform a wavelength scan of the structure in order to find the absorption
 * and extinction cross sections. To request a wavelength scan, the output block needs to be similar to:
 * 
 * @code{.xml}
 * <output type="response">
 *  <scan type="A+E">
 *    <wavelength initial="1400" final="2400" steps="101" />
 *  </scan>
 * </output>
 * @endcode
 * 
 * Here, a wavelength scan is defined starting from 1400 nanometers, up to 2400 nanometers and suing 101 steps. 
 * It is worth noting that, in the case of a scan, any wavelength defined in the source block will be
 * ignored. Output will occur in two separate .dat files named after the simulation name and containing the
 * keywords AbsCrossSection and ExtCrossSection.
 * 
 * It is also possible to request a radius scan using the following block:
 * 
 * @code{.xml}
 * <output type="response">
 *  <scan type="A+E">
 *    <radius initial="1000" final="2000" steps="101" />
 *  </scan>
 * </output>
 * @endcode
 * 
 * The format is identical to a wavelength scan but the radius of the spheres is varied at each
 * simulation step. Furthermore, a double scan of both radius and wavelength is possible using the
 * following combination:
 * 
 * @code{.xml}
 * <output type="response">
 *  <scan type="A+E">
 *    <wavelength initial="1400" final="2400" steps="101" />
 *    <radius initial="1000" final="2000" steps="101" />
 *  </scan>
 * </output>
 * @endcode
 * 
 * Here, the simulation will run for each possible (radius, wavelength) combination, with the output
 * being a 2D matrix for both the Absorption and Extinction cross sections.
 */
