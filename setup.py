from setuptools import setup
import os

setup(
    name='exocomp',
    version='0.1.0',
    py_modules=['exocomp'],
    author='Joshua D. Lothringer',
    author_email='jlothringer@stsci.edu',
    description='A tool for exoplanet atmospheric abundance computations.',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    long_description_content_type='text/markdown',
    url='https://github.com/jlothringer/exocomp',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy',
        'pandas',
        'easychem',
        'scipy',
        'matplotlib',
        'random'
    ],
)
