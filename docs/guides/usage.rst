Usage
=====

This section describes the general usage of `HTSinfer`.

General Usage
-------------

.. code-block:: bash

   htsinfer [--output-directory PATH]
            [--temporary-directory PATH]
            [--cleanup-regime {DEFAULT,KEEP_ALL,KEEP_NONE,KEEP_RESULTS}]
            [--records INT]
            [--threads INT]
            [--transcripts FASTA]
            [--read-layout-adapters PATH]
            [--read-layout-min-match-percentage FLOAT]
            [--read-layout-min-frequency-ratio FLOAT]
            [--library-source-min-match-percentage FLOAT]
            [--library-source-min-frequency-ratio FLOAT]
            [--library-type-max-distance INT]
            [--library-type-mates-cutoff FLOAT]
            [--read-orientation-min-mapped-reads INT]
            [--read-orientation-min-fraction FLOAT]
            [--tax-id INT]
            [--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}]
            [-h] [--version]
            PATH [PATH]

The above command allows the user to infer metadata for single- or paired-ended RNA-Seq libraries by specifying file paths and relevant parameters. The tool outputs metadata in JSON format to :code:`STDOUT` and logs to :code:`STDERR`.

Command-line Options
---------------------

Available command-line parameters are categorized as follows:

.. list-table:: General Options
   :widths: 40 60

   * - :code:`--output-directory`
     - Path where output data will be saved. 
   * - :code:`--temporary-directory`
     - Path for storing temporary files generated during execution.
   * - :code:`--cleanup-regime`
     - Specifies which data should be kept after completion. Options are: :code:`DEFAULT, KEEP_ALL, KEEP_NONE, KEEP_RESULTS`
   * - :code:`--verbosity`
     - Controls the verbosity level of log output. Options are: :code:`DEBUG, INFO, WARN, ERROR, CRITICAL`
   * - :code:`-h, --help`
     - Show help screen and exit.
   * - :code:`-v, --version`
     - Show version information and exit.


.. list-table:: Processing and Performance Options
   :widths: 40 60

   * - :code:`--records`
     - Limits the number of input records to process; setting this to 0 will process all records.
   * - :code:`--threads`
     - Specifies the number of threads for STAR to optimize performance.
   * - :code:`--tax-id`
     - Taxonomy ID for the sample source, aiding in organism-specific analyses.

.. list-table:: Library-specific Options
   :widths: 40 60

   * - :code:`PATH [PATH]`
     - Path(s) to the RNA-Seq input data. For paired-end libraries, provide paths to both mate files.
   * - :code:`--transcripts`
     - Path to the FASTA file containing transcript sequences for reference.
   * - :code:`--read-layout-adapters`
     - Path to a file with 3' adapter sequences (one sequence per line) used to identify adapter content.
   * - :code:`--read-layout-min-match-percentage`
     - Minimum percentage of reads containing an adapter for\n it to be considered as the library’s 3’-end adapter.
   * - :code:`--read-layout-min-frequency-ratio`
     - Minimum frequency ratio between the most and second most frequent adapters to select the 3’-end adapter.
   * - :code:`--library-source-min-match-percentage`
     - Minimum percentage of reads aligning with a library source for it to be considered representative of the library.
   * - :code:`--library-source-min-frequency-ratio`
     - Minimum frequency ratio between primary and secondary library sources,
       ensuring only the most prominent source is identified.
   * - :code:`--library-type-max-distance`
     - Maximum allowable distance between read pairs to classify the library type.
   * - :code:`--library-type-mates-cutoff`
     - Ratio cutoff to determine the consistency of mate orientation in paired-end reads.
   * - :code:`--read-orientation-min-mapped-reads`
     - Minimum number of mapped reads to ensure reliable inference of read orientation.
   * - :code:`--read-orientation-min-fraction`
     - Minimum fraction (must exceed 0.5) of reads supporting a given orientation to confirm its accuracy.
