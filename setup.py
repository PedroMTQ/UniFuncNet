import pathlib
from setuptools import setup,find_packages


package_name='DRAX'
def get_package_version():
    import os
    from sys import platform
    if platform.startswith('win'):
        SPLITTER = '\\'
    else:
        SPLITTER = '/'


    dir_name=os.path.dirname(os.path.abspath(__file__))
    init_path=f'{dir_name}{SPLITTER}{package_name.lower()}{SPLITTER}__init__.py'
    package_version=None
    with open(init_path) as file:
        for line in file:
            if '__version__' in line:
                package_version=line.replace('__version__','')
                package_version=package_version.strip('\n')
                package_version=package_version.strip()
                package_version=package_version.strip('=')
                package_version=package_version.strip()
                package_version=package_version.strip('"')
                package_version=package_version.strip('"')
    return package_version



# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text(encoding='utf-8')

long_description='DRAX is network annotation tool that integrates multiple biological databases into one composite output.'

setup(
    name=package_name,
    version=get_package_version(),
    author="Pedro Queirós",
    author_email="pdqueiros@gmail.com",
    description="Biological network annotation tool.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/PedroMTQ/drax",
    project_urls={
        "Bug Tracker": "https://github.com/PedroMTQ/DRAX/issues",
    },
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license="MIT",
    #install_requires=['requests','python>=3.7'],
    entry_points={
        "console_scripts": [
            "drax=drax.__main__:main",
        ],
    },
)
