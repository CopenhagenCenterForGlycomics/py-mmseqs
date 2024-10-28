from setuptools import setup, find_packages

setup(
    name='py-mmseqs',
    version='0.0.1',
    url='https://github.com/CopenhagenCenterForGlycomics/py-mmseqs.git',
    author='Hiren Joshi',
    author_email='hirenj@gmail.com',
    description='Wrapper around mmseqs2',
    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    package_dir={"": "src"},
    install_requires=[],
)
