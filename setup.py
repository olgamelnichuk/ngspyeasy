try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'a description',
    'author': 'an author',
    'url': 'an url',
    'download_url': 'a download url',
    'author_email': 'an email',
    'version': '0.1',
    'install_requires': [],
    'packages': ['ngspyeasy'],
    'scripts': [],
    'name': 'NGSPyEasy'
}

setup(**config)