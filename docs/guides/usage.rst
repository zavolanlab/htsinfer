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

.. table:: General Options
    :widths: 25, 75

.. table:: General Options
   :widths: 25, 75

   | Option | Description |
   |---|---|
   | :code:`--output-directory PATH` | Path where output data will be saved. |
   | :code:`--temporary-directory PATH` | Path for storing temporary files generated during execution. |
   | :code:`--cleanup-regime {DEFAULT,KEEP_ALL,KEEP_NONE,KEEP_RESULTS}` | Specifies which data should be kept after completion. |
   | :code:`--verbosity {DEBUG,INFO,WARN,ERROR,CRITICAL}` | Controls the verbosity level of log output. |
   | :code:`-h, --help` | Show help screen and exit. |
   | :code:`--version` | Show version information and exit. |

.. table:: Library-specific Options
   :widths: 25, 75

   | Option | Description |
   |---|---|
   | `PATH [PATH]` | Path(s) to the RNA-Seq input data. For paired-end libraries, provide paths to both mate files. |
   | `--transcripts FASTA` | Path to the FASTA file containing transcript sequences for reference. |
   | `--read-layout-adapters PATH` | Path to a file with 3' adapter sequences (one sequence per line) used to identify adapter content. |
   | `--read-layout-min-match-percentage FLOAT` | Minimum percentage of reads containing an adapter for it to be considered as the library’s 3’-end adapter. |
   | `--read-layout-min-frequency-ratio FLOAT` | Minimum frequency ratio between the most and second most frequent adapters to select the 3’-end adapter. |
   | `--library-source-min-match-percentage FLOAT` | Minimum percentage of reads aligning with a library source for it to be considered representative of the library. |
   | `--library-source-min-frequency-ratio FLOAT` | Minimum frequency ratio between primary and secondary library sources, ensuring only the most prominent source is identified. |
   | `--library-type-max-distance INT` | Maximum allowable distance between read pairs to classify the library type. |
   | `--library-type-mates-cutoff FLOAT` | Ratio cutoff to determine the consistency of mate orientation in paired-end reads. |
   | `--read-orientation-min-mapped-reads INT` | Minimum number of mapped reads to ensure reliable inference of read orientation. |
   | `--read-orientation-min-fraction FLOAT` | Minimum fraction (must exceed 0.5) of reads supporting a given orientation to confirm its accuracy. |

.. table:: Processing and Performance Options
   :widths: 25, 75

   | Option | Description |
   |---|---|
   | `--records INT` | Limits the number of input records to process; setting this to 0 will process all records. |
   | `--threads INT` | Specifies the number of threads for concurrent processing to optimize performance. |
   | `--tax-id INT` | Taxonomy ID for the sample source, aiding in organism-specific analyses. |
