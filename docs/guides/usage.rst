Usage
=====

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

The above command allows the user to infer metadata for single- or paired-ended RNA-Seq libraries by specifying file paths and relevant parameters. The tool outputs metadata in JSON format to `STDOUT` and logs to `STDERR`.

Command-line Options
---------------------

Available command-line parameters are categorized as follows:

- **General Options**: These include specifying directories, verbosity level, and other global settings.
- **Library-specific Options**: These parameters allow the user to modify settings related to the input data, such as transcript references, adapter sequences, and match thresholds.
- **Output Options**: These settings control the output format, including the number of records and the output destination.
- **Meta Options**: The user can also control the behavior of the tool with meta options such as cleanup regimes, thread count, and version information.

For a complete list of all available options, use the following command:

.. code-block:: bash

   htsinfer --help
