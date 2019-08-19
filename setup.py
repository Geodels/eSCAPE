import numpy
from numpy.distutils.core import setup, Extension
try:
    from numpy.distutils.fcompiler import FCompiler
    def runtime_library_dir_option(self, dir):
        return self.c_compiler.runtime_library_dir_option(dir)
    FCompiler.runtime_library_dir_option = \
        runtime_library_dir_option
except Exception:
    pass

def configuration(parent_package='',top_path=None):
    INCLUDE_DIRS = []
    LIBRARY_DIRS = []
    LIBRARIES    = []

    # PETSc
    import os
    PETSC_DIR  = os.environ['PETSC_DIR']
    PETSC_ARCH = os.environ.get('PETSC_ARCH', '')
    from os.path import join, isdir
    if PETSC_ARCH and isdir(join(PETSC_DIR, PETSC_ARCH)):
        INCLUDE_DIRS += [join(PETSC_DIR, PETSC_ARCH, 'include'),
                         join(PETSC_DIR, 'include')]
        LIBRARY_DIRS += [join(PETSC_DIR, PETSC_ARCH, 'lib')]
    else:
        if PETSC_ARCH: pass
        INCLUDE_DIRS += [join(PETSC_DIR, 'include')]
        LIBRARY_DIRS += [join(PETSC_DIR, 'lib')]
    LIBRARIES += [#'petscts', 'petscsnes', 'petscksp',
                  #'petscdm', 'petscmat',  'petscvec',
                  'petsc']

    import os
    import petsc4py
    INCLUDE_DIRS += [petsc4py.get_include()]

    # Configuration
    from numpy.distutils.misc_util import Configuration

    config = Configuration('', parent_package, top_path)

    config.add_extension('eSCAPE._fortran',
                        sources = ['fortran/functions.pyf',
                                     'fortran/functions.F90'],
                        depends = ['fortran/functionsmodule.h'],
                        # f2py_options=['--quiet'],
                        define_macros=[], #[('F2PY_REPORT_ON_ARRAY_COPY',0)],
                        include_dirs=INCLUDE_DIRS + [os.curdir],
                        libraries=LIBRARIES,
                        library_dirs=LIBRARY_DIRS,
                        extra_f90_compile_args=['-fPIC','-O3'],
                        #extra_f90_compile_args = ['-fPIC', '-O0', '-g', '-fbacktrace','-fcheck=all'],
                        #extra_link_args = ['-shared '],
                        runtime_library_dirs=LIBRARY_DIRS)
                        #define_macros=[],
                        #include_dirs=[os.curdir],
                        # libraries=[],
                        # library_dirs=[],
                        # extra_f90_compile_args=['-fPIC','-O3'],
                        # extra_link_args = ['-shared '],
                        # runtime_library_dirs=[])

    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup
    setup(name = 'eSCAPE',
          author            = "Tristan Salles  ",
          author_email      = "tristan.salles@sydney.edu.au",
          url               = "https://github.com/Geodels/eSCAPE",
          version           = "0.1",
          description       = "Scalable Parallelised Landscape Evolution Model",
          # ext_modules       = [ext],
          configuration     = configuration,
          packages          = ['eSCAPE', 'eSCAPE.tools','eSCAPE.mesher','eSCAPE.pit','eSCAPE.flow'],
          # package_data      = {'eSCAPE': ['Examples/data',
          #                                   'Examples/Notebooks/*.ipynb',
          #                                   'Examples/Python/*.py']},
          classifiers       = ['Programming Language :: Python :: 2',
                               'Programming Language :: Python :: 2.6',
                               'Programming Language :: Python :: 2.7',
                               'Programming Language :: Python :: 3',
                               'Programming Language :: Python :: 3.3',
                               'Programming Language :: Python :: 3.4',
                               'Programming Language :: Python :: 3.5',
                               'Programming Language :: Python :: 3.6']
          )
