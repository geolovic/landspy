#   Licensed under the MIT License

import setuptools

with open("README.md", "r", encoding = "utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "landspy",
    version = "1.0",
    author = "J. Vicente Perez Pena",
    author_email = "geolovic@gmail.com",
    description = "Python library for topographic analysis from DEMs",
    long_description = long_description,
    url =  "https://github.com/geolovic/landspy",
    project_urls = {
            "Bug Tracker" : "https://github.com/geolovic/landspylandspy/issues"
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
    ],

    package_dir = {"": "src"},
    packages = setuptools.find_packages(where="src"),
    python_requires = ">=3.6",
    install_requires=["GDAL>=3", "numpy", "matplotlib", "scikit-image", "scipy"],
    tests_require=['unittest']
)
