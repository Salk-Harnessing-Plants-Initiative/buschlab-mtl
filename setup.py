from distutils.core import setup

setup(
    name='mtl',
    version='0.2.0',
    packages=['mtl', 'mtl.pygwas_modules', 'mtl.core'],
    url='',
    license='',
    author='christian.goeschl',
    author_email='',
    description='multi trait gwas using limix', requires=['h5py', 'limix']
)
