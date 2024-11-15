Examples
========

`HTSinfer` provides easy-to-use commands for analyzing single- and paired-ended RNA-Seq libraries.


Single-ended library example
----------------------------

To run `HTSinfer` on a single-ended RNA-Seq library, use the following command:

.. code-block:: bash

   htsinfer tests/files/adapter_single.fastq

Paired-ended library example
----------------------------

To run `HTSinfer` on a paired-ended RNA-Seq library, use the following command:

.. code-block:: bash

   htsinfer tests/files/adapter_1.fastq tests/files/adapter_2.fastq

Both commands will output the results in JSON format to :code:`STDOUT` and the log to :code:`STDERR`.

Example output
--------------

Here is a sample output for the paired-ended library:

.. code-block:: json

   {
      "library_stats": {
         "file_1": {
            "read_length": {
               "min": 75,
               "max": 75,
               "mean": 75.0,
               "median": 75,
               "mode": 75
            }
         },
         "file_2": {
            "read_length": {
               "min": 75,
               "max": 75,
               "mean": 75.0,
               "median": 75,
               "mode": 75
            }
         }
      },
      "library_source": {
         "file_1": {
            "short_name": "hsapiens",
            "taxon_id": "9606"
         },
         "file_2": {
            "short_name": "hsapiens",
            "taxon_id": "9606"
         }
      },
      "library_type": {
         "file_1": "first_mate",
         "file_2": "second_mate",
         "relationship": "split_mates"
      },
      "read_orientation": {
         "file_1": "SF",
         "file_2": "SR",
         "relationship": "ISF"
      },
      "read_layout": {
         "file_1": {
            "adapt_3": "AATGATACGGCGACC",
            "polyA_frac": 10.0
         },
         "file_2": {
            "adapt_3": "AATGATACGGCGACC",
            "polyA_frac": 10.0
         }
      }
   }

Results
-------

For more details on the output structure, refer to the :code:`Results` model in the API documentation.

- :code:`library_stats`: Read length statistics results, including the minimum, maximum, mean, median and mode.
- :code:`library_source`: Library source inference results, with the short name and NCBI taxonomy ID of the source organism.
- :code:`library_type`: Library type inference results, single- or paired-end. In case of paired-end samples, the mate designation.
- :code:`read_orientation`: Read orientation inference results, based on the fragment library types notation from `Salmon <https://salmon.readthedocs.io/en/latest/library_type.html>`_.
- :code:`read_layout`: Read layout inference results, including the 3'-adapter sequence and the poly(A) fraction. 