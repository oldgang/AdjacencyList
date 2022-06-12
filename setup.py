from distutils.core import setup, Extension

module = Extension('simple_graphs', sources=['simple_graphs.c'])

setup(
    name = 'simple_graphs',
    version = '0.0.1',
    ext_modules = [module]
)
