"""Infers single or paired end library for a given file.

This module parses fastq files and depending on whether they are
first mate or second mate, it returns tuple of values accordingly.
Tuple description :
    First member of the tuple corresponds to the result for the
    first file with possible values- 0,1,2,3,None.
    (None) if no file provided.
    Second member of the tuple corresponds to the result for the
    second file with possible values- 0,1,2,3,None.
    (None) if no file provided.
    Third member of the tuple indicates 4 if the files have split paired
    end library and 5 if the read id's don't match among the two files.

The modules parses sucessfully for fastq files that are strictly
as per as wiki convention.

Examples:
    @HWUSI-EAS100R:6:73:941:1973#0/1
    @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

"""
import sys
import gzip
import re
from functools import partial
from enum import Enum
from Bio.SeqIO.QualityIO import FastqGeneralIterator
sys.tracebacklimit = 0


class Messages(Enum):
    invalid_file = 0
    first_mate = 1
    second_mate = 2
    mixed_paired_end = 3
    split_paired_end = 4
    different_read_ids = 5


class EndParser:
    def __init__(self, file1_name: str, file2_name: str):
        """Constructor to assign filenames and initialize file pointers.

        Args:
            file1_name (str) : First file name provided by the user.
            file2_name (str) : Second file name provide by the user.

        Attributes:
            file1 : File pointer to first file.
            file2 : File pointer to second file.
        """
        self.file1_name = file1_name
        self.file2_name = file2_name
        self.file1 = None
        self.file2 = None

    def validate(self, seq_id: str):
        """Validates the type of identifier for the seq_id.

        Args:
            seq_id (str) : Sequence identifier.

        Returns:
            int : 0 , 1 , 2
            Depending on the sequence id whether it matches the
            wiki conventions it returns 0 for invalid identifier
            and 1 and 2 for first mate read and second mate read
            respectively.

        """
        first_convention = re.compile(r'^[^:#/ ]+:\d+:\d+:\d+:\d+#\d+/[1-2]\s*')
        second_convention = re.compile(r'^[^:#/ ]+:\d+:[^:#/ ]+:\d+:\d+:\d+:\d+ [1-2]:(Y|N):\d*[02468]:([ATCG]+|\d+)\s*')
        match_1 = first_convention.finditer(seq_id)
        match_2 = second_convention.finditer(seq_id)
        len_match_1 = None
        len_match_2 = None
        len_seq_id = len(seq_id)
        for match in match_2:
            len_match_2 = match.end() - match.start()
            break
        for match in match_1:
            len_match_1 = match.end() - match.start()
            break
        if len_match_1 is not None:
            if len_seq_id == len_match_1:
                return 1
            else:
                return 0
        else:
            if len_match_2 is None:
                return 0
            else:
                if len_seq_id == len_match_2:
                    return 2
                else:
                    return 0

    def check_file1(self):
        """Opening and checking errors for first file by the user.

        Attribute Error:
            Function can raise error if file doesn't exist or
            something wrong with opening of file.

        """
        try:
            _open = partial(
                gzip.open, mode='rt') if self.file1_name.endswith(".gz") else open
            self.file1 = _open(self.file1_name)
        except OSError:
            print("Could not open/read file:", self.file1_name)
            sys.exit()

    def check_file2(self):
        # Opening and checking errors for second file by the user
        try:
            _open = partial(
                gzip.open, mode='rt') if self.file2_name.endswith(".gz") else open
            self.file2 = _open(self.file2_name)
        except OSError:
            print("Could not open/read file:", self.file2_name)
            sys.exit()

    def read_line_in_fastq(self, seq_id: str, seq_id_dict: dict):
        """This function processes for every seq_id in a fastq file.

        Args:
            seq_id (str): Sequence identifier.
            seq_id_dict (dict) : Dictionary to maintain frequency
            of every seq_id.

        Atrributes:
            id_type (int): Determines to which wiki convention it corresponds.

        """
        seq_id.strip(' ')
        id_type = self.validate(seq_id)
        if id_type == 0:
            return 0
        else:
            if id_type == 1:
                split_id = list(map(str, seq_id.split('/')))
                if split_id[0] in seq_id_dict:
                    seq_id_dict[split_id[0]] += 1
                else:
                    seq_id_dict[split_id[0]] = 1
                return int(split_id[1])
            elif id_type == 2:
                split_id = list(map(str, seq_id.split()))
                if split_id[0] in seq_id_dict:
                    seq_id_dict[split_id[0]] += 1
                else:
                    seq_id_dict[split_id[0]] = 1
                return int(split_id[1][0])

    def result_from_count(self, count: list):
        """This function on basis of count values returns final result.

        Returns:
            int : 0 , 1 , 2 , 3 where each points to it's corresponding
            memeber in the enumerator.

        """
        if sum(count) == 0:
            return 0
        if count[0]*count[1] == 0:
            if count[1] == 0:
                return 1
            elif count[0] == 0:
                return 2
        else:
            return 3

    def read_file1(self, seq_id_dict: dict):
        """Function reads file1 and calls other helper functions.

        Args:
            seq_id_dict (dict): Dictionary to maintain frequency of every
            seq_id.

        Attributes:
            count (list) : List to maintain count of number of 1's and
                2's that will decide further whether first mate or
                second mate.
            result (int) : Stores final result among the enumerator values
                whether invalid, first mate or second mate.

        Function calls self.read_line_in_fast() to further process
        the seq_id, and calls self.result_from_count to get result.

        """

        self.check_file1()
        count = [0]*2
        for seq_id, seq, qual in FastqGeneralIterator(self.file1):
            x = self.read_line_in_fastq(seq_id, seq_id_dict)
            if x == 0:
                return 0
            else:
                count[x-1] += 1
        result = self.result_from_count(count)
        return result

    def read_file2(self, seq_id_dict: dict):
        # Function reads file2 and calls other functions.
        self.check_file2()
        count = [0]*2
        for seq_id, seq, qual in FastqGeneralIterator(self.file2):
            x = self.read_line_in_fastq(seq_id, seq_id_dict)
            if x == 0:
                return 0
            else:
                count[x-1] += 1
        result = self.result_from_count(count)
        return result

    def fastq(self):
        """This function calls self.read_file to read fastq files.

        Attributes:
            seq_id_dict1 (dict) : Dictionary to maintain frequency of every
            seq_id in file1.
            seq_id_dict2 (dict) : Dictionary to maintain frequency of every
            seq_id in file2.

            result1 (int) : Stores final result among the enumerator
            values for the first file.
            result2 (int) : Stores final result among the enumerator
            values for the second file.

        Returns:
            tuple : a tuple of length 3 is returned.

        Examples: for return type
            (1,2,4): first mate , second mate , split paired end.
            (1,2,5): first mate , second mate , read id's don't
            match among the two files.

        """
        seq_id_dict1 = dict()
        seq_id_dict2 = dict()
        result1 = None
        result2 = None
        tuple_third_element = None
        files = 0
        if self.file1_name is not None:
            files += 1
            result1 = self.read_file1(seq_id_dict1)
            self.file1.close()

        if self.file2_name is not None:
            files += 1
            result2 = self.read_file2(seq_id_dict2)
            self.file2.close()

        if files == 2 and result1*result2 == 2:
            if seq_id_dict1 == seq_id_dict2:
                tuple_third_element = 4
            else:
                tuple_third_element = 5
        return (result1, result2, tuple_third_element)
