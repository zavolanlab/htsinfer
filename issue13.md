## Task:
Given a set of sequences, calculate the number of occurrences of all motifs within a range of 
length that occur in the input sequences. Return the results as a dictionary.

## Function description:

1) for every sequence: determine all motifs/k-mers (k in range(a,b))  in the sequences with string slicing, save as list 	
2) make one big list combining all the smaller lists 
3) construct an output dictionary with motif sequences (key) and corresponding frequencies (value)
4) output sorted dictionary

## Function outline: 

def count_motifs(input_sequences: list, min_motif_length: int, max_motif_length: int) -> dict
""" Function that calculates the occurrence of all possible motifs in one or multiple sequences and 
    returns a dictionary with all motifs within the specified length and their occurrence

Args: 
    input_sequences (list): list of sequences
    min_motif_length (int): minimal length of motif
    max_motif_length (int): maximal length of motif
   
Returns:
    motif frequency (dict) with paired data {"motif_seq": frequency }
    
Raises:
    ValueError: input_sequences contains non-standard (ACTG) characters
    TypeError: input_sequences is not a list
"""

# Pull request

### Description

This function calculates all the occuring motifs in a set of sequences in a specified length and returns them as a dict.

Fixes # 13

### Type of change

Please delete options that are not relevant.

- [ ] Bug fix (non-breaking change which fixes an issue)
- [x] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality
      to not work as expected)
- [ ] Documentation update

### Checklist

Please carefully read these items and tick them off if the statements are true
_or do not apply_.

- [x] I have performed a self-review of my own code
- [x] My code follows the existing coding style, lints and generates no new
      warnings
- [x] I have added type annotations to all function/method signatures, and I
      have added type annotations for any local variables that are non-trivial,
      potentially ambiguous or might otherwise benefit from explicit typing.
- [x] I have commented my code in hard-to-understand areas
- [x] I have added ["Google-style docstrings"] to all new modules, classes,
      methods/functions or updated previously existing ones
- [x] I have added tests that prove my fix is effective or that my feature
      works
- [x] New and existing unit tests pass locally with my changes and I have not
      reduced the code coverage relative to the previous state
- [x] I have updated any sections of the app's documentation that are affected
      by the proposed changes

If for some reason you are unable to tick off all boxes, please leave a
comment explaining the issue you are facing so that we can work on it
together.
