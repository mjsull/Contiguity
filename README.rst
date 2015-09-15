Contiguity
==========

Contiguity is a tool for constructing and visualising assembly graphs.
It uses a linear layout so that the assembly graph can be directly compared 
to a reference.

The main website for contiguity can be found [here](http://mjsull.github.io/Contiguity)

Contiguity can be downloaded from [here](http://mjsull.github.io/Contiguity/files.html)

For detailed information on how to use Contiguity please see the manual_
otherwise see Quick Start below.

Quick Start
-----------

Supported formats & the CAG
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Contiguity works with ABySS_ (**.dot**), Velvet_ (**LastGraph**), Newbler_ 
(**.ace**), GFA_ and SPAdes_ (**FASTG**) formats.

For all other assemblies, an assembly graph (.cag) can be created from the 
Contiguity GUI (file->create cag file) or using the command line. 

**We recommend that both Velvet and SPAdes graphs be reconstructed using 
the Contiguity CAG format.**

You can read more about the CAG in the manual_.


Generation of the CAG file
~~~~~~~~~~~~~~~~~~~~~~~~~~

You can generate a CAG from the command line like::

    $ Contiguity -cl -c <contig_file.fa> -fq <read_file.fq> -o <output_folder>

Or by selecting File > Create CAG in the GUI and providing a read file  and a contig file

This assumes:
    * (~8GB of free memory)
    * contig_file.fa: is in FASTA file of contigs or scaffolds
    * read_file.fq: Interleaved fastq file - read1_left, read1_right, read2_left 
      etc... orientated as such --> <--
    * output_folder: folder to put output files in, can and will overwrite 
      files in this folder.

The resultant CAG file (output_folder/assembly.cag) can then be loaded into 
Contiguity.


The Contiguity GUI
~~~~~~~~~~~~~~~~~~

There are 3 main functions of the Contiguity GUI -

**Visualising an assembly graph**:
    * Load FASTG/LastGraph/CAG etc. using "File->Load assembly"
    * View assembly graph using "View->View Assembly"

**Compare assembly graph to a reference**:
    * Load FASTG/LastGraph/CAG etc. using "File->Load assembly"
    * Create comparison to a reference by selecting "File->Create Comparison"
    * Select a reference file and click ok, when asked if you want to 
      generate a comparison, click "yes".
    * View assembly graph using "View->View Assembly"

**Self Comparison (Long read assemblies)**:
    * Load assembly (FASTA) using "File->Load assembly"
    * Create a self comparison using "View->Self Comparison"
    * Select "OK" and when asked if you want to generate a comparison, click 
      "yes"

For more in depth description of functionality and work flows please see the 
manual_.


Citation
--------

If you use Contiguity in your work, please cite it using::

    Sullivan MJ, Ben Zakour NL, Forde BM, Stanton-Cook M, Beatson SA. (2015)
    Contiguity: Contig adjacency graph construction and visualisation.
    PeerJ PrePrints 3:e1273 https://dx.doi.org/10.7287/peerj.preprints.1037v1



.. _manual: https://github.com/mjsull/Contiguity/wiki
.. _pip installation instructions: http://pip.readthedocs.org/en/latest/installing.html
.. _ABySS: http://www.bcgsc.ca/platform/bioinfo/software/abyss 
.. _Velvet: https://www.ebi.ac.uk/~zerbino/velvet/
.. _Newbler: http://www.454.com/products/analysis-software/
.. _SPAdes: http://bioinf.spbau.ru/spades
.. _GFA: https://github.com/pmelsted/GFA-spec
