from setuptools import setup, find_packages

setup(
    name='py_mmseqs',
    version='0.0.1',
    url='https://github.com/CopenhagenCenterForGlycomics/py-mmseqs.git',
    author='Hiren Joshi',
    author_email='hirenj@gmail.com',
    description='Wrapper around mmseqs2',
    packages=["py_mmseqs"],
    package_dir={"": "src"},
    package_data={"py_mmseqs": ["**/*.json"]},
    install_requires=[],
)
