#!/usr/bin/env python

import sys
import os
import io

import distutils.sysconfig as dsc
# from distutils import core, dir_util
# from distutils.core import setup,Extension
# from distutils.command import build_ext
from distutils.sysconfig import get_config_var
# from distutils.ccompiler import show_compilers
# use setuptools
from setuptools import setup, Extension
try:
    from Cython.Distutils import build_ext
    # from Cython.Build import cythonize
except ImportError:
    from distutils.command import build_ext

pyincdir = dsc.get_python_inc(plat_specific=1)
pylibdir = os.path.join('/', *pyincdir.split('/')[:-2] + ['lib'])
try:
    import numpy as np
    np_include_dir = np.get_include()
    # include_dirs.append(np_include_dir)
except ImportError:
    np_include_dir = os.path.join(
                pylibdir.replace('lib/python', 'local/lib/python'),
                'numpy', 'core', 'include')
    print("Unable to import numpy, trying header %s".format(np_include_dir))
    sys.exit(1)

from glob import glob
pj = os.path.join

# from distutils.command.install_data import install_data

try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

PY3K = sys.version_info[0] > 2

# Remove unwanted options
_UNWANTED_OPTS = frozenset(['-Wstrict-prototypes'])
os.environ['OPT'] = ' '.join(
                    _ for _ in get_config_var('OPT').strip().split()
                    if _ not in _UNWANTED_OPTS)

hierarchy = Extension(
    "biTK.ngs.cluster.hierarchy",
    include_dirs=[np_include_dir],
    sources=[
        'biTK/ngs/cluster/_hierarchy.pyx'
        ])

gibbs_sampler = Extension(
    "biTK.ngs.statistics.sampling._gibbs_sampler",
    include_dirs=[np_include_dir],
    sources=[
        'biTK/ngs/statistics/sampling/_gibbs_sampler.pyx'
        ])


def split_multiline(value):
    """Split a multiline string into a list, excluding blank lines."""
    val = value
    if not isinstance(value, (str)):
        val = value.encode('utf-8')
    return [element for element in
            (line.strip() for line in val.split('\n'))
            if element]


def cfg_to_args(path='setup.cfg'):
    """Compatibility helper to use setup.cfg in setup.py.

    This functions uses an existing setup.cfg to generate a dictionnary of
    keywords that can be used by distutils.core.setup(**kwargs).  It is used
    by generate_setup_py.

    *file* is the path to the setup.cfg file.  If it doesn't exist,
    PackagingFileError is raised.
    """

    # XXX ** == needs testing
    D1_D2_SETUP_ARGS = {"name": ("metadata",),
                        "version": ("metadata",),
                        "author": ("metadata",),
                        "author_email": ("metadata",),
                        "maintainer": ("metadata",),
                        "maintainer_email": ("metadata",),
                        "url": ("metadata", "home_page"),
                        "description": ("metadata", "summary"),
                        "long_description": ("metadata", "description"),
                        "download-url": ("metadata",),
                        "classifiers": ("metadata", "classifier"),
                        "platforms": ("metadata", "platform"),  # **
                        "license": ("metadata",),
                        "requires": ("metadata", "requires_dist"),
                        "provides": ("metadata", "provides_dist"),  # **
                        "obsoletes": ("metadata", "obsoletes_dist"),  # **
                        "package_dir": ("files", 'packages_root'),
                        "packages": ("files",),
                        "scripts": ("files",),
                        "resources": ("files",),
                        "py_modules": ("files", "modules"),  # **
                        "package_data": ("files", "package_data"),
                        }

    MULTI_FIELDS = ("classifiers",
                    "platforms",
                    "requires",
                    "provides",
                    "obsoletes",
                    "packages",
                    "scripts",
                    "py_modules",
                    "extension",
                    "package_data",
                    )

    def has_get_option(config, section, option):
        if config.has_option(section, option):
            return config.get(section, option)
        elif config.has_option(section, option.replace('_', '-')):
            return config.get(section, option.replace('_', '-'))
        else:
            return False

    # The real code starts here
    config = RawConfigParser()
    # f = codecs.open(path, encoding='utf-8')
    f = io.open(path, encoding="utf-8")
    try:
        config.readfp(f)
    finally:
        f.close()

    kwargs = {}
    for arg in D1_D2_SETUP_ARGS:
        if len(D1_D2_SETUP_ARGS[arg]) == 2:
            # The distutils field name is different than distutils2's
            section, option = D1_D2_SETUP_ARGS[arg]

        else:
            # The distutils field name is the same thant distutils2's
            section = D1_D2_SETUP_ARGS[arg][0]
            option = arg

        in_cfg_value = has_get_option(config, section, option)
        if not in_cfg_value:
            # There is no such option in the setup.cfg
            if arg == 'long_description':
                filenames = has_get_option(config, section, 'description-file')
                if filenames:
                    filenames = split_multiline(filenames)
                    in_cfg_value = []
                    for filename in filenames:
                        # fp = codecs.open(filename, encoding='utf-8')
                        fp = io.open(filename, encoding="utf-8")
                        try:
                            in_cfg_value.append(fp.read())
                        finally:
                            fp.close()
                    in_cfg_value = '\n\n'.join(in_cfg_value)
            else:
                continue

        if arg == 'package_dir' and in_cfg_value:
            in_cfg_value = {'': in_cfg_value}

        if arg in MULTI_FIELDS:
            # support multiline options
            in_cfg_value = split_multiline(in_cfg_value)
            if arg == 'packages' and in_cfg_value:
                if 'package_dir' in kwargs:
                    if kwargs['package_dir']['']:
                        in_cfg_value = [
                            kwargs['package_dir']['']+'.'+pack
                            for pack in in_cfg_value]
            if arg == 'package_data' and in_cfg_value:
                datafiles = {}
                for line in in_cfg_value:
                    split_path = line.split('/')
                    package_strcut_str = ".".join(split_path[:-1]).strip()
                    files = [f.split('/')[-1] for f in glob(line)]
                    datafiles[package_strcut_str] = files
                in_cfg_value = datafiles
                kwargs['install_package_data'] = True

            if arg == 'resources' and in_cfg_value:
                datafiles = []
                for line in in_cfg_value:
                    file_str = line.split('=')[0].strip()
                    datafiles.extend(
                        [f.split('/')[-1] for f in glob(file_str)])
                # kwargs['data_files'] = [('share/bilab',
                #    [line.split('=')[0].strip() for line in in_cfg_value])]
                # kwargs['data_files'] = [('biTK/ngs',datafiles)]

        kwargs[arg] = in_cfg_value
    return kwargs

# setup
general_settings = cfg_to_args()
# Key: resources has to be removed
# general_settings.pop('resources')
general_settings['zip-safe'] = False
general_settings['ext_modules'] = [hierarchy, gibbs_sampler]
# setup requires for external pacakges
general_settings['install_requires'] = [
    'matplotlib>1.4.0',
    'joblib>=0.9.0b4',
    'statsmodels>=0.6']
# print(general_settings['install_requires'])
setup(**general_settings)
# setup(**cfg_to_args())
