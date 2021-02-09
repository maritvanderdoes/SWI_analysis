from setuptools import setup, find_packages

contrib = ['Marit van der Does']

# setup.
setup(
    name='swi-analysis',
    version='0.0.0',
    description='SWI analysis in python',
    author=', '.join(contrib),
    license='MIT',
    packages=find_packages(exclude=[
        'tests',
    ]),
    install_requires=[
        'numpy>=1.19',
        'matplotlib>=3.3',
        'scikit-image>=0.17',
        'scipy>=1.5',
        'seaborn>=0.11',
        'pandas>=1.1',
        'luigi>=3.0'
    ],
)
