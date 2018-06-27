from setuptools import setup, find_packages
from distutils.util import convert_path

# Get version
main_ns = {}
version_path = convert_path("ctg/version.py")
with open(version_path) as handle:
	exec(handle.read(), main_ns)

# Set up
setup(
	name="ctg",
	version=main_ns['__version__'],
	description="Composition and Time-course aware Genetic interaction analysis.",
	url="https://github.com/bpmunson/ctg",
	author="Brenton Munson",
	author_email="bpmunson@eng.ucsd.edu",
	liscense="MIT",
	classifiers=[
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 3 - Alpha',

		# Indicate who your project is intended for
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',

		# Pick your license as you wish (should match "license" above)
		'License :: OSI Approved :: MIT License',

		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
	], 
	entry_points={
		'console_scripts': [
			'ctg=ctg:main',
		],

	},
	keywords='genetic interaction crispr',
	#packages=['ctg'],
	packages=find_packages(),
	install_requires=['numpy','scipy','pysam','statsmodels','ConfigArgParse'],
	include_package_data=True,
)
