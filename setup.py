from setuptools import setup

setup(name='miRBaseMiner',
      version='0.1',
      description='Mining the miRNA annotation in miRBase for comprehensive understanding in miRNA annotation reference before implementing in miRNA study.',
      url='https://github.com/joey0214/miRBaseMiner',
      classifiers=[
        'Development Status :: 0.1 - Alpha',
        'License :: OSI Approved :: GNU GPLv3 License',
        'Programming Language :: Python :: 2.7',
      ],
      keywords='miRBase microRNA miRNA',
      author='Xiangfu (Joey) Zhong',
      author_email='joey.zhong.cn@gmail.com',
      license='GNU GPLv3',
      packages=['miRBaseMiner'],
      install_requires=[
          'biopython',
          'editdistance',
          'ftplib'
      ],
      zip_safe=True)


