import platform
from setuptools import setup, Extension
from Cython.Distutils import build_ext

NAME = "orthogonalstreamlines"
VERSION = "0.1.2"
DESCR = "Create an interconnected cable mesh from an orientation field on a triangulated 3D surface"
KEYWORDS = "mesh,cables,surface,streamlines"
URL = "http://github.com/jacquemv/orthogonalstreamlines"
REQUIRES = ['numpy', 'cython', 'surfacetopology', 'evenlyspacedstreamlines']
AUTHOR = "Vincent Jacquemet"
EMAIL = "vincent.jacquemet@umontreal.ca"
LICENSE = "MIT"
SRC_DIR = "orthogonalstreamlines"
PACKAGES = [SRC_DIR, 
            SRC_DIR+'/intersection', 
            SRC_DIR+'/tessellation', 
            SRC_DIR+'/triangulation']

if platform.system() == 'Windows':
    compiler_args = ['/openmp', '/O2']
    linker_args = []
else:
    compiler_args = ['-fopenmp', '-O3', '-Wall']
    linker_args = ['-fopenmp']

FILES = {
    "intersection.runengine": ("/intersection/runengine.pyx",
                               "/intersection/intersection.cpp"),
    "tessellation.tessellation": ("/tessellation/tessellation.pyx",),
    "triangulation.triangulation": ("/triangulation/triangulation.pyx",
                                    "/triangulation/triangulatefacets.cpp")
}
EXTENSIONS = []
for target, src_files in FILES.items():
    ext = Extension(SRC_DIR + "." + target,
        sources=[SRC_DIR + file for file in src_files],
        libraries=[],
        extra_compile_args=compiler_args,
        extra_link_args=linker_args,
        language="c++",
        include_dirs=[SRC_DIR, SRC_DIR + "/common"]
    )
    ext.cython_directives = {'language_level': "3"}
    EXTENSIONS.append(ext)

setup(install_requires=REQUIRES,
      packages=PACKAGES,
      zip_safe=False,
      name=NAME,
      version=VERSION,
      description=DESCR,
      keywords=KEYWORDS,
      long_description=open('README.md', 'r').read(),
      long_description_content_type='text/markdown',
      author=AUTHOR,
      author_email=EMAIL,
      url=URL,
      license=LICENSE,
      cmdclass={"build_ext": build_ext},
      ext_modules=EXTENSIONS
)
