Examples
========

HTSinfer provides easy-to-use commands for analyzing single- and paired-ended RNA-Seq libraries.

Single-ended Library Example
----------------------------

To run HTSinfer on a single-ended RNA-Seq library, use the following command:

.. code-block:: bash

   htsinfer tests/files/adapter_single.fastq

Paired-ended Library Example
----------------------------

To run HTSinfer on a paired-ended RNA-Seq library, use the following command:

.. code-block:: bash

   htsinfer tests/files/adapter_1.fastq tests/files/adapter_2.fastq

Both commands will output the results in JSON format to `STDOUT` and the log to `STDERR`.

Example Output
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

For more details on the output structure, refer to the `Results` model in the API documentation.
