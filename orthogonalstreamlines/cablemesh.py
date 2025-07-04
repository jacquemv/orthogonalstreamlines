import time
from contextlib import contextmanager
from collections import namedtuple, Counter
import numpy as np
from surfacetopology import surface_topology
from evenlyspacedstreamlines import evenly_spaced_streamlines
from . import intersection
from . import tessellation
from . import triangulation

SPACES = '    '
UNIT = 'cm'
RATIO_DX_RADIUS = 1.38

OrthogonalStreamlinesMesh = namedtuple('OrthogonalStreamlinesMesh', [
    'cables', 'nlc', 'ntc', 'dx', 'vertices', 'triangles', 
    'facets', 'boundaries', 'ver_to_orig_tri', 'tri_to_facet',
    'neighbors', 'sign', 'info', 'random_seed', 'is_node'
])

#-----------------------------------------------------------------------------
def create_orthogonal_streamlines_mesh(
        vertices, triangles, orientation, dx, *,
        nb_seeds=1024, options=None,
        random_seed=None, verbose=True, add_ghost_nodes=False,
        unit='cm') -> OrthogonalStreamlinesMesh:
    """Generate an interconnected cable mesh from evenly spaced orthogonal 
    streamlines

    Args:
        vertices (nv-by-3 array): x, y, z coordinates of the nv vertices 
            of type numpy.float64
        triangles (nt-by-3 int array): indices of the vertices of the nt 
            triangles of type numpy.int32
        orientation (nt-by-3 array): orientation vector in each triangle
        dx (float or tuple): target mesh resolution in the same unit as 
            'vertices'; if dx is a tuple, resolution is different in the 
            longitudinal (first value) and transverse direction (second 
            value)
        nb_seeds (int): number of seed points for streamline generation
            (default: 1024)
        options (dict): additional arguments passed to the function
            evenly_spaced_streamlines
        random_seed (tuple of int): set random seed for streameline 
            generation; there is one seed for longitudinal streamlines and 
            one for transverse streamlines
        add_ghost_nodes (bool): add the vertices of the streamlines to the 
            cables; these "ghost" nodes are intended to facilitate 
            visualization and create curves instead of segments between 
            consecutive nodes of the cable, but they are not evenly spaced
            (default: False)
        verbose (bool): print informations during computations (default: True)
        unit (str): unit of distance in the input mesh (used only for
            displaying information); default: 'cm'
    
    Returns:
        namedtuple containing the following fields:
            cables (list of int arrays): each cable is represented by an 
                array of vertex indices; longitudinal cables are listed first
            nlc (int): number of longitudinal cables
            ntc (int): number of transverse cables
            dx (float array of size ne): length of each edge
            vertices (nv-by-3 float array): cable vertex positions
            triangles (nt-by-3 int array): new triangulation of the surface
            facets (list of int arrays): facets[n] is a k-by-n array whose 
                k rows are facets with n sides
            boundaries (list of int arrays): boundaries[n] is the array of
                vertex indices of the n-th boundary or hole
            ver_to_orig_tri (int array of size nv): maps cable vertices to 
                triangle indices from the original triangulated surface
            tri_to_facet (int array of size nt): maps triangle indices to
                facet indices
            neighbors (nv-by-4 int array): indices of up to 4 neighboring 
                vertices of each vertex; the first two columns are for 
                neighbors in the longitudinal direction, the next two columns
                for the transverse direction; index is -1 where there is no 
                neighbor
            sign (int array of size nv): sign (with respect to the vector 
                normal to the surface) of the cross product of the 
                tangent vectors of the streamlines at their intersection
            random_seed (tuple of int): random seeds used for streamlines 
                generation
            is_node (int array of size nv): = 0 if the vertex is a ghost node,
                and 1 otherwise
            info (dict): internal information about the calculations
    """
    UNIT = unit
    info = {}

    # check matrix format

    vertices = np.ascontiguousarray(vertices, dtype=np.float64)
    triangles = np.ascontiguousarray(triangles, dtype=np.int32)
    orientation = np.ascontiguousarray(orientation, dtype=np.float64)

    # check geometry

    with timer('original surface', verbose):
        topo = surface_topology(triangles)
        if len(topo) != 1:
            raise ValueError('The triangular mesh must have one single '+
                             f'connected component ({len(topo)} identified)')
        if not topo[0].manifold:
            raise ValueError('The triangular mesh is not a manifold')
        if not topo[0].oriented:
            raise ValueError('The triangular mesh is not consistently '+
                             'oriented')
        mesh = OrthogonalStreamlines(vertices, triangles, orientation, 
                                     topo[0].n_boundaries)
    if verbose:
        print(SPACES+f'{vertices.shape[0]} vertices, '
              f'{triangles.shape[0]} triangles\n')
    info['nb_boundaries'] = topo[0].n_boundaries

    # generate longitudinal streamlines

    with timer('longitudinal streamlines', verbose):
        if options is None:
            options= {} 
        if random_seed is not None:
            options['random_seed'] = random_seed[0]
        mesh.generate_streamlines(dx, direction='long', nb_seeds=nb_seeds, 
                                  options=options)
    L = mesh.infos_long.lengths
    info['streamline_info_long'] = mesh.infos_long
    info['stats_streamlines_long'] = calc_stats(L)
    if verbose:
        print(SPACES+f'{L.size} curves')
        print_stats('length', info['stats_streamlines_long'], UNIT)
        print(SPACES+f'random seed: {mesh.infos_long.random_seed}\n')
    
    # generate transverse streamlines

    with timer('transversal streamlines', verbose):
        if random_seed is not None:
            options['random_seed'] = random_seed[1]
        mesh.generate_streamlines(dx, direction='trans', nb_seeds=nb_seeds,
                                  options=options)
    L = mesh.infos_trans.lengths
    info['streamline_info_trans'] = mesh.infos_trans
    info['stats_streamlines_trans'] = calc_stats(L)
    if verbose:
        L = mesh.infos_trans.lengths
        print(SPACES+f'{L.size} curves')
        print_stats('length', info['stats_streamlines_trans'] , UNIT)
        print(SPACES+f'random seed: {mesh.infos_trans.random_seed}\n')
    
    # compute streamline intersections

    with timer('streamline intersections', verbose):
        mesh.create_cable_mesh(add_ghost_nodes)
        dx_long, dx_trans = mesh.space_steps()
    info['stats_dx_long'] = calc_stats(dx_long)
    info['stats_dx_trans'] = calc_stats(dx_trans)
    if verbose:
        print(SPACES+f'{mesh.cablenet.vertices.shape[0]} vertices')
        print(SPACES+f'{mesh.cablenet.nc_long} long. '+
              f'and {mesh.cablenet.nc_trans} trans. cables')
        print_stats('dx(long)', info['stats_dx_long'] , UNIT, more_space=1)
        print_stats('dx(trans)', info['stats_dx_trans'] , UNIT)
        print()
    info['cnt_loose_ends'] = mesh.cablenet.cnt_loose_ends
    info['cnt_empty_cables'] = mesh.cablenet.cnt_empty_cables
    info['cnt_duplicates'] = mesh.cablenet.cnt_duplicates
    info['cnt_isolated_vertices'] = mesh.cablenet.cnt_isolated_vertices
    info['nb_connected_components'] = mesh.cablenet.nb_connected_components

    # tessellate the cable mesh

    with timer('tessellation', verbose):
        mesh.tessellate()
        hist = mesh.histogram_facets()
    if verbose:
        print(SPACES+f'{len(mesh.facets)} facets')
        print(SPACES+'# of n-gon facets for each n:')
        print(SPACES[:-1], hist, '\n')
    info['histogram_facets'] = hist

    # triangulate the facets

    with timer('triangulation', verbose):
        mesh.triangulate()
        wrong_orient = mesh.fix_orientation()
    if verbose:
        print(SPACES+f'{mesh.triangles.shape[0]} triangles')
        failed_loops = mesh.triangulation_failures
        if failed_loops:
            count = 0
            for loop in failed_loops:
                count += len(np.unique(loop)) < len(loop)
            n1, n2 = len(failed_loops), count
            s1, s2 = 's' * (n1>1), 's' * (n2>1)
            print(SPACES+f'{n1} facet triangulation failure{s1}'
                  f' ({n2} non-simple polygon{s2})')
        if wrong_orient > 0:
            print(SPACES+f'orientation fixed in {wrong_orient} triangles')
        print()
    info['triangulation_failure'] = mesh.triangulation_failures

    # process boundaries

    with timer('boundaries', verbose):
        facets, boundaries, permutation = tessellation.group_facets_by_size(
            mesh.facets, mesh.cutoff)
        assert np.unique(permutation).size == permutation.size
        inv_permutation = np.zeros_like(permutation)
        inv_permutation[permutation] = np.arange(permutation.size)
        tri_to_facet = inv_permutation[mesh.facetid]
    if verbose:
        print(SPACES+f'{mesh.nb_boundaries} boundaries\n')
    
    # create and output the data structure

    return OrthogonalStreamlinesMesh(
        cables=intersection.unpack_cables(mesh.cablenet.cables, 
                                          mesh.cablenet.cables_len),
        nlc=mesh.cablenet.nc_long, 
        ntc=mesh.cablenet.nc_trans,
        dx=np.concatenate((dx_long, dx_trans)),
        vertices=mesh.cablenet.vertices,
        triangles=mesh.triangles,
        facets=facets,
        boundaries=boundaries,
        ver_to_orig_tri=mesh.cablenet.indices_tri,
        tri_to_facet=tri_to_facet,
        neighbors=mesh.vneigh,
        sign=mesh.cablenet.sign,
        is_node=mesh.cablenet.is_node,
        random_seed=(mesh.infos_long.random_seed, 
                     mesh.infos_trans.random_seed),
        info=info
    )

#-----------------------------------------------------------------------------
def calc_stats(x):
    return dict(min=x.min(), max=x.max(), mean=x.mean(), std=x.std(),
                quartiles=np.percentile(x, [25, 50, 75]))

#-----------------------------------------------------------------------------
def print_stats(name, stat, unit, more_space=0):
    spaces = SPACES + ' ' * len(name) + '  ' + ' ' * more_space
    print(SPACES+name+': ' + ' ' * more_space + 
          f"mean = {stat['mean']:.4f}, std = {stat['std']:.4f} " + unit)
    q = stat['quartiles']
    print(spaces + 
          f'q1 =   {q[0]:.4f}, med = {q[1]:.4f}, q3 = {q[2]:.4f} ' + unit)
    print(spaces + 
          f"min =  {stat['min']:.4f}, max = {stat['max']:.4f} " + unit)

#-----------------------------------------------------------------------------
@contextmanager
def timer(block_name, verbose=True):
    if not verbose:
        yield
        return
    s = block_name + ': ' + '-' * (68 - len(block_name)) + ' '
    print(s, end='')
    t0 = time.perf_counter()
    yield
    print(f'{time.perf_counter() - t0:.3f} s')

#-----------------------------------------------------------------------------
OrientedSurface = namedtuple('OrientedSurface', ['ver', 'tri', 'orient'])

#-----------------------------------------------------------------------------
class OrthogonalStreamlines:

    def __init__(self, vertices, triangles, orientation, nb_boundaries):
        self.trisurf = OrientedSurface(ver=vertices, tri=triangles, 
                                       orient=orientation)
        self.nb_boundaries = nb_boundaries
    
    #-------------------------------------------------------------------------
    def generate_streamlines(self, dx, direction, nb_seeds=1024, 
                             options=None):
        if isinstance(dx, (tuple, list)):
            # the distance between longitudinal cables determines transverse 
            # resolution, so the longitudinal radius is based on the 
            # transverse dx (when direction == 'long' is equal to 1)
            # and reciprocally
            radius = dx[direction == 'long'] / RATIO_DX_RADIUS
        else:
            radius = dx / RATIO_DX_RADIUS
        if options is None:
            options = {}

        output = evenly_spaced_streamlines(
            self.trisurf.ver, self.trisurf.tri, self.trisurf.orient, 
            radius=radius, orthogonal=(direction != 'long'),
            seed_points=nb_seeds, **options
        )
        if direction == 'long':
            self.lines_long, self.faces_long, self.infos_long = output
        else:
            self.lines_trans, self.faces_trans, self.infos_trans = output
        
    #-------------------------------------------------------------------------
    def create_cable_mesh(self, add_ghost_nodes=False):
        normals = intersection.tri_normals(self.trisurf.ver, self.trisurf.tri)
        self.cablenet = intersection.create_cable_network(
            normals,
            self.lines_long, self.faces_long,
            self.lines_trans, self.faces_trans, 
            add_ghost_nodes=add_ghost_nodes
        )
    
    #-------------------------------------------------------------------------
    def space_steps(self):
        dx = np.linalg.norm(self.cablenet.vertices[self.cablenet.cables[1:]]
                           -self.cablenet.vertices[self.cablenet.cables[:-1]], 
                            axis=-1)
        cable_end = np.cumsum(self.cablenet.cables_len)
        dx = np.delete(dx, cable_end[:-1]-1)
        sep = np.sum(self.cablenet.cables_len[:self.cablenet.nc_long]-1)
        return dx[:sep], dx[sep:]
    
    #-------------------------------------------------------------------------
    def tessellate(self):
        self.facets, self.vneigh = tessellation.tessellate(
            self.cablenet.cables, self.cablenet.cables_len, 
            self.cablenet.nc_long, self.cablenet.sign, 
            return_moves=True, return_neigh=True)
    
    #-------------------------------------------------------------------------
    def triangulate(self):
        largest, self.cutoff = tessellation.find_largest_facets(
            self.facets, self.nb_boundaries)
        self.triangles, self.facetid, self.triangulation_failures = \
            triangulation.triangulate_facets(self.cablenet.vertices, 
                                             self.facets, self.cutoff)
    
    #-------------------------------------------------------------------------
    def fix_orientation(self):
        calc_normals = intersection.tri_normals
        normals = calc_normals(self.cablenet.vertices, self.triangles)
        orig_normals = calc_normals(self.trisurf.ver, self.trisurf.tri)
        idx_orig = self.cablenet.indices_tri[self.triangles[:, 0]]
        I =  np.sum(normals * orig_normals[idx_orig], axis=1) <= 0
        self.triangles[I] = self.triangles[I][:, [0, 2, 1]]
        return I.sum()
    
    #-------------------------------------------------------------------------
    def histogram_facets(self):
        sizes = tessellation.facet_size(self.facets)
        S = dict(Counter(sizes))
        return {n: S[n] for n in sorted(S.keys())[:-self.nb_boundaries]}