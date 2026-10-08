#   Licensed under the MIT License

import setuptools

with open("README.md", "r", encoding = "utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "landspy",
    version = "1.4.1",
    author = "J. Vicente Perez Pena",
    author_email = "geolovic@gmail.com",
    description = "Python library for landscape analysis with DEMs",
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url =  "https://github.com/geolovic/landspy",
    project_urls = {
            "Bug Tracker" : "https://github.com/geolovic/landspy/issues"
    },
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
    ],

    package_dir = {"": "src"},
    packages = setuptools.find_packages(where="src"),
    python_requires = ">=3.10",
    install_requires=["GDAL", "numpy", "matplotlib", "scikit-image", "scipy", "shapely>=2", "numba>=0.53"],
    include_package_data=False,
)
