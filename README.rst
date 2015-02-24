Contiguity
==========

A tool for visualising assemblies.

.. image:: https://github.com/BeatsonLab-MicrobialGenomics/Contiguity/blob/master/docs/manual/Contiguity_SS.png
    :alt: Contiguity Screen shot
    :align: center


Requirements:
    * Python 2.7+
    * NCBI-BLAST (needed for overlap edge creation and automatic comparison 
      generation)
    * Bowtie 2 (needed for paired end edge creation)


Installation
------------

If you're not familiar with the command-line we recommend you ask local IT 
support to help you install Contiguity.


Checking requirements are installed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will need to install/have installed:
    * ncbiblast >= ???
    * python >= 2.7 (**Python 3 is not supported**)
    * bowtie2 >= ???

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

For information on how to use Contiguity please see the manual_


Citation
--------

If you use Contiguity in your work, please cite it using::

    Mitchell J Sullivan, Nouri Ben Zakour, Brian Forde, Mitchell Stanton-Cook & Scott A Beatson*
    Contiguity: Contig adjacency graph construction and visualisation
    https://github.com/BeatsonLab-MicrobialGenomics/Contiguity


.. _manual: https://github.com/BeatsonLab-MicrobialGenomics/Contiguity/blob/master/docs/manual/Contiguity_manual_0.3.pdf
.. _pip: http://www.pip-installer.org/en/latest/
.. _pip installation instructions: http://pip.readthedocs.org/en/latest/installing.html
