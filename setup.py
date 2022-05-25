from setuptools import setup

setup(
    name='genespector',
    version='0.1.0',
    description='Package for interactive viewing of AnnData objects/files',
    license='MIT',
    packages=['genespector'],
    author='Iwo Kucinski',
    author_email='',
    url='https://github.com/Iwo-K/genespector',
    include_package_data=True,
    package_data={'': ['data/tiny_example2.h5ad']}
)
