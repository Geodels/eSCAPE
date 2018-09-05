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

    import os
    import petsc4py
    INCLUDE_DIRS += [petsc4py.get_include()]

    # Configuration
    from numpy.distutils.misc_util import Configuration

    config = Configuration('', parent_package, top_path)

    config.add_extension('gSCAPE._fortran',
                        sources = ['fortran/functions.pyf',
                                     'fortran/functions.F90'],
                        depends = ['fortran/functionsmodule.h'],
                        # f2py_options=['--quiet'],
                        define_macros=[],
                        include_dirs=[os.curdir],
                        libraries=[],
                        library_dirs=[],
                        extra_f90_compile_args=['-fPIC','-O3'],
                        extra_link_args = ['-shared '],
                        runtime_library_dirs=[])

    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup
    setup(name = 'gSCAPE',
          author            = "Tristan Salles  ",
          author_email      = "tristan.salles@sydney.edu.au",
          url               = "https://github.com/Geodels/gSCAPE",
          version           = "0.1",
          description       = "Scalable Parallelised Landscape Evolution Model",
          # ext_modules       = [ext],
          configuration     = configuration,
          packages          = ['gSCAPE', 'gSCAPE.tools','gSCAPE.mesher','gSCAPE.pit','gSCAPE.flow'],
          # package_data      = {'gSCAPE': ['Examples/data',
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
