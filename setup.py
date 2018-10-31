from setuptools import setup

setup(name='miRBaseMiner',
      version='0.2',
      description='Mining the miRNA annotation in miRBase for comprehensive understanding in miRNA annotation reference before implementing in miRNA study.',
      long_description='Mining the miRNA annotation in miRBase for comprehensive understanding in miRNA annotation reference before implementing in miRNA study.',
      url='https://github.com/joey0214/miRBaseMiner',
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2.7',
      ],
      keywords='miRBase microRNA miRNA',
      author='Xiangfu (Joey) Zhong',
      author_email='joey.zhong.cn@gmail.com',
      license='GNU GPLv3',
      packages=['miRBaseMiner'],
      install_requires=[
          'biopython',
          'editdistance'
      ],
      python_requires='>=2.6, !=3.0.*, !=3.1.*, !=3.2.*, <4',
      zip_safe=True)


