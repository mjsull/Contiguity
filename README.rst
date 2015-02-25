Contiguity
==========

Contiguity is a tool for constructing and visualising assembly graphs.
It uses a linear layout so that the assembly graph can be directly compared to a reference.

.. image:: https://pypip.in/version/Contiguity/badge.svg
        :target: https://pypi.python.org/pypi/Contiguity/
        :alt: Latest Version

.. image:: https://pypip.in/download/Contiguity/badge.svg
        :target: https://pypi.python.org/pypi/Contiguity/
        :alt: Downloads

.. image:: https://travis-ci.org/BeatsonLab-MicrobialGenomics/Contiguity.svg?branch=master
        :target: https://travis-ci.org/BeatsonLab-MicrobialGenomics/Contiguity
        :alt: Build status


.. image:: https://github.com/BeatsonLab-MicrobialGenomics/Contiguity/blob/master/docs/manual/Contiguity_SS.png
    :alt: Contiguity Screen shot
    :align: center


Requirements:
    * Python 2.7+
    * NCBI-BLAST+ (needed for overlap edge creation and automatic comparison 
      generation)
    * Bowtie 2 (needed for paired end edge creation)


Installation
------------

If you're not familiar with the command-line we recommend you ask local IT 
support to help you install Contiguity.


Checking requirements are installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will need to install/have installed:
    * ncbiblast >= 2.2.28
    * python >= 2.7 (**Python 3 is not supported**)
    * bowtie2 >= 2.1.0

You can check these are installed by::
    
    $ python --version
    $ blastn -version
    $ bowtie2 --version

Installation of python, blastn or bowtie2 (without a package manager) is 
beyond the scope of this document.

If you have both python, blastn and bowtie2 you need to (if not already 
present) install pip_.

You can check if pip_ exists with::

    $ which pip

If you get a "not found", please read the `pip installation instructions`_. 

**If you already have pip we do suggest you upgrade it.** We are using version 
1.5.6 at the time of writing this document. 

You can upgrade pip_ like this::

    $ pip install --upgrade pip


pip based installation of Contiguity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have root/admin something like::

    $ pip install Contiguity

Otherwise (not root/admin or permission denied errors running above)::

    $ pip install --user Contiguity

If you installed using the --user option of pip_, Contiguity will typically 
end up in: /home/$USER/.local/bin/ 
You need to add this location to you ~/.bash_profile. 

Add Contiguity to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile
    $ source !$


Testing the installation of Contiguity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run (in the Terminal)::
    
    $ Contiguity


Upgrading Contiguity
~~~~~~~~~~~~~~~~~~~~

You can upgrade like this::
    
    pip install --upgrade Contiguity


**Please regularly check back to make sure you're running the most recent 
Contiguity version.**


Usage/Docs
----------

For detailed information on how to use Contiguity please see the manual_
otherwise see Quick Start below.


Citation
--------

If you use Contiguity in your work, please cite it using::

    Mitchell J Sullivan, Nouri Ben Zakour, Brian Forde, Mitchell Stanton-Cook & Scott A Beatson*
    Contiguity: Contig adjacency graph construction and visualisation
    https://github.com/BeatsonLab-MicrobialGenomics/Contiguity

Quick Start
-----------

Constructing an assembly graph
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Contiguity works with Abyss (.dot), Velvet (LastGraph), newbler (.ace) and Spades (FASTG) formats.
For all other assemblies, a assembly graph (.cag) can be created from the GUI (file->create cag file)
or using the command line. We recommend both Velvet and Spades graphs be reconstructed using Contiguity.

USAGE:

    $ Contiguity.py -cl -c <contig_file.fa> -fq <read_file.fq> -o <output_folder>

REQUIREMENTS: With default settings Contigutiy requires at least 8gb of free memory (RAM)

contig file: FASTA file of contigs or scaffolds

read file: Interleaved fastq file - read1_left, read1_right, read2_left etc... orientated as such --> <--

output folder: folder to put output files in, can and will overwrite files in this folder.
Will create folder if folder doesn't exist. Default folder (GUI) is .contiguity_wd

The CAG file (output_folder/assembly.cag) can then be loaded into Contiguity.

The Contiguity GUI has three main functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visualising an assembly graph:
    * Load FASTG/LastGraph/CAG etc. using "File->Load assembly"
    * View assembly graph using "View->View Assembly"

Compare assembly graph to a reference:
    * Load FASTG/LastGraph/CAG etc. using "File->Load assembly"
    * Create comparison to a reference by selecting "File->Create Comparison"
    * Select a reference file and click ok, when asked if you want to generate a comparison, click "yes".
    * View assembly graph using "View->View Assembly"

Self Comparison (Long read assemblies):
    * Load assembly (FASTA) using "File->Load assembly"
    * Create a self comparison using "View->Self Comparison"
    * Select "OK" and when asked if you want to generate a comparison, click "yes"

For more in depth description of functionality and workflows please see the manual.

.. _manual: https://github.com/BeatsonLab-MicrobialGenomics/Contiguity/blob/master/docs/manual/Contiguity_manual.pdf
.. _pip: http://www.pip-installer.org/en/latest/
.. _pip installation instructions: http://pip.readthedocs.org/en/latest/installing.html
