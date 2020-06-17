import os
import sys
import gzip
import subprocess as sp
from itertools import islice
sys.tracebacklimit = 0

messages = {
    1: "Single End or First_mate \n",
    2: "Second_mate \n",
    3: "Mixed \n",
    5: "First_mate or Second_mate \n"
}


class End_parser:

    """ Infer Single or Paired """

    def __init__(self):
        print("\nRunning Parser\n")

    def SRA_contents(self, name, filename):

        """ Runs fastq-dump """

        try:
            with open(os.devnull, "w") as fnull:
                contents = sp.check_output(
                    ["fastq-dump", "-X", "1", "-Z", "--split-spot", name],
                    stderr=fnull)
        except sp.CalledProcessError:
            raise Exception("Error running fastq-dump on", filename)
        contents = contents.decode('utf-8')
        return contents

    def isPairedSRA(self, name, filename):

        """ Checks if SRA file is Single end or Paired """

        string_from_function = self.SRA_contents(name, filename)
        string_from_file_4 = ""
        string_from_file_8 = ""

        if filename.endswith(".gz"):
            file_in = gzip.open(filename, 'rt')
        else:
            file_in = open(filename)

        i = 1
        for lin in islice(file_in, 8):
            string_from_file_8 += lin
            if i == 4:
                string_from_file_4 = string_from_file_8
            i += 1
        file_in.close()
        count = 0
        for i in string_from_function:
            if ord(i) == 10:
                count += 1
        if count == 4:
            if string_from_function == string_from_file_4:
                return 1
            else:
                return 0

        elif count == 8:
            inp = list(map(str, string_from_function.split("\n")))
            f_h = inp[0] + "\n" + inp[1] + "\n" + inp[2] + "\n" + inp[3] + "\n"
            s_h = inp[4] + "\n" + inp[5] + "\n" + inp[6] + "\n" + inp[7] + "\n"
            if string_from_function == string_from_file_8:
                return 3
            elif f_h == string_from_file_4:
                return 1
            elif s_h == string_from_file_4:
                return 2
            else:
                return 0

    def validate_1(self, line):

        """Validate for First type Wiki Convention:
            @HWUSI-EAS100R:6:73:941:1973#0/1"""

        if line[0] == '@':
            line = line[1:]
            c = True
            inp = list(map(str, line.split(':')))
            if len(inp) == 5:
                p = list(map(str, inp[4].split('#')))
                inp[4] = p[0]

                for i in range(1, 5, 1):
                    if not inp[i].isnumeric():
                        return False

                if len(p) == 2:
                    c = list(map(str, p[1].split('/')))
                    if len(c) == 1:
                        if c[0].isnumeric():
                            return True
                        else:
                            return False
                    elif len(c) == 2:
                        if c[0].isnumeric():
                            if c[1].isnumeric():
                                if c[1] == '1' or c[1] == '2':
                                    return True
                                else:
                                    return False
                            else:
                                return False
                        else:
                            return False
                    else:
                        return False
                elif len(p) == 1:
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False

    def validate_2(self, line):

        """ Validate for Second Type Wiki Convention :
             @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG """

        if line[0] == '@':
            line = line[1:]
            inp = list(map(str, line.split()))
            if len(inp) == 2:
                c = list(map(str, inp[0].split(':')))
                d = list(map(str, inp[1].split(':')))
                if len(c) == 7 and len(d) == 4:
                    e = []
                    e.append(c[1])
                    e.append(c[3])
                    e.append(c[4])
                    e.append(c[5])
                    e.append(c[6])
                    for i in e:
                        if not i.isnumeric():
                            return False
                    if d[0].isnumeric():
                        if d[0] == '1' or d[0] == '2':
                            if d[1] == 'Y' or d[1] == 'N':
                                if d[2].isnumeric():
                                    if int(d[2]) % 2 == 0:
                                        return True
                                    else:
                                        return False
                                else:
                                    return False
                            else:
                                return False
                        else:
                            return False
                    else:
                        return False
                else:
                    return False
        else:
            return False

    def read_line_in_fastq(self, line, count, type, Seq_count):

        """ Update Count whether First Mate or Second Mate """

        inp = list(map(str, line.split()))

        if type == 1:
            if inp[0][len(inp[0])-2] == '/':
                count[int(line[len(line) - 1]) - 1] += 1
                if line[:len(line) - 2] in Seq_count:
                    Seq_count[line[:len(line) - 2]] += 1
                else:
                    Seq_count[line[:len(line) - 2]] = 1
            else:
                if line in Seq_count:
                    Seq_count[line] += 1
                else:
                    Seq_count[line] = 1
        else:
            if inp[0] in Seq_count:
                Seq_count[inp[0]] += 1
            else:
                Seq_count[inp[0]] = 1
            count[int(inp[1][0]) - 1] += 1

    def read_fastq(self, file, Seq_count, filename, emp_str):

        """ Reads FastQ file """

        count = [0]*2
        for i, line in enumerate(file):
            line = line[:-1]
            if i % 4 == 0 and line:
                x = line
                inp = list(map(str, x.split()))
                x = ""
                for p in inp:
                    x += p
                    x += " "
                x = x[:len(x) - 1]
                if len(inp) == 1:
                    if self.validate_1(x):
                        self.read_line_in_fastq(x, count, 1, Seq_count)
                    else:
                        return 0
                elif len(inp) == 2:
                    if self.validate_2(x):
                        self.read_line_in_fastq(x, count, 2, Seq_count)
                    else:
                        return 0
                elif len(inp) == 3:
                    q = list(map(str, x.split('.')))
                    name = q[0][1:]
                    # emp_str = name
                    return self.isPairedSRA(name, filename)
        # 5 - first/second mate
        # 3 - mixed
        # 2 - second mate
        # 1 - first mate
        # 0 - Invalid file
        c1 = 0
        c2 = 0
        if count[0]*count[1] == 0:
            if count[0] + count[1] == 0:
                for i in Seq_count:
                    if Seq_count[i] == 1:
                        c1 += 1
                    if Seq_count[i] == 2:
                        c2 += 1
                n = len(Seq_count)
                if n == 0:
                    return 0
                if n == c1:
                    return 5
                elif n == c2:
                    return 3
                else:
                    return 3
            else:
                if count[0] == 0:
                    return 2
                else:
                    return 1
        else:
            return 3

    def check_file(self, file1_name, file2_name):

        """ Checks if file exists and has proper extensions """

        if file1_name is None and file2_name is None:
            raise SystemError("Specify at least one file name.")

        if file1_name is not None:
            if os.path.isfile(file1_name) is False:
                raise FileNotFoundError(file1_name)
            if (file1_name.endswith(".fastq.gz") or
                    file1_name.endswith(".fastq")) is False:
                raise TypeError("Invalid file type : {}".format(file1_name))

        if file2_name is not None:
            if os.path.isfile(file2_name) is False:
                raise FileNotFoundError(file2_name)
            if (file2_name.endswith(".fastq.gz") or
                    file2_name.endswith(".fastq")) is False:
                raise TypeError("Invalid file type : {}".format(file2_name))

    def print_mssg(
                   self, parse1, parse2, Seq_count1,
                   Seq_count2, emp_str1, emp_str2, file1_name, file2_name):

        """Prints whether Single or Paired"""
        file_info = 0
        if parse2 == -1:
            if parse1 == 0:
                print(
                    "\nRead IDs do not adhere to Illumina or SRA format",
                    " : {} \n".format(file1_name))
            else:
                print(messages[parse1])

        else:
            if parse1*parse2 == 0:
                if parse1 == 0:
                    print("\nRead IDs do not adhere to Illumina or SRA format",
                          " : {} \n".format(file1_name))
                if parse2 == 0:
                    print("\nRead IDs do not adhere to Illumina or SRA format",
                          " : {} \n".format(file2_name))
                return

            elif parse1 == 1 or parse1 == 5:
                if parse2 == 2 or parse2 == 5:
                    if emp_str1 == "" and emp_str2 == "":
                        if Seq_count1 == Seq_count1:
                            print("\nSplit Paired End Library\n:")
                            print("{}: First Mate\n".format(file1_name))
                            print("{}: Second Mate\n".format(file2_name))
                            file_info = 1
                        else:
                            print("\nError: Read ID's don't match.",
                                  "Files are of different Experiments \n")
                    else:
                        if emp_str1 == emp_str2:
                            print("\nSplit Paired End Library:\n")
                            print("{}: First Mate\n".format(file1_name))
                            print("{}: Second Mate\n".format(file2_name))
                            file_info = 1
                        else:
                            print("\nError: Read ID's don't match.",
                                  "Files are of different Experiments \n")

            elif parse1 == 2 or parse1 == 5:
                if parse2 == 1 or parse2 == 5:
                    if emp_str1 == "" and emp_str2 == "":
                        if Seq_count1 == Seq_count1:
                            print("\nSplit Paired End Library:\n")
                            print("{}: First Mate\n".format(file2_name))
                            print("{}: Second Mate\n".format(file1_name))
                            file_info = 1
                        else:
                            print("\nError: Read ID's don't match.",
                                  "Files are of different Experiments \n")
                    else:
                        if emp_str1 == emp_str2:
                            print("\nSplit Paired End Library:\n")
                            print("{}: First Mate\n".format(file2_name))
                            print("{}: Second Mate\n".format(file1_name))
                            file_info = 1
                        else:
                            print("\nError: Read ID's don't match.",
                                  "Files are of different Experiments \n")

            if file_info == 0:
                print(
                    "\n{} : {} \n{} : {}\n".format(
                        file1_name, messages[parse1],
                        file2_name, messages[parse2]))

    def fastq(self, file1_name=None, file2_name=None):

        """ Open Fastq files """

        self.check_file(file1_name, file2_name)

        parse1 = -1
        parse2 = -1
        Seq_count1 = dict()
        Seq_count2 = dict()
        emp_str1 = ""
        emp_str2 = ""

        if file1_name is not None:
            if file1_name.endswith(".gz"):
                file = gzip.open(file1_name, 'rt')
            else:
                file = open(file1_name)
            parse1 = self.read_fastq(file, Seq_count1, file1_name, emp_str1)
            file.close()

        if file2_name is not None:
            if file2_name.endswith(".gz"):
                file = gzip.open(file2_name, 'rt')
            else:
                file = open(file2_name)
            parse2 = self.read_fastq(file, Seq_count2, file2_name, emp_str2)
            file.close()

        self.print_mssg(parse1, parse2, Seq_count1, Seq_count2, emp_str1,
                        emp_str2, file1_name, file2_name)
        return parse1, parse2


def main():
    pass
