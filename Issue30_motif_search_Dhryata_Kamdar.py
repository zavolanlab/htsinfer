#!/usr/bin/env python
# coding: utf-8

# ## 1. EXTRACT READ SEQUENCES FROM FASTA FILE

# 1.1 FUNCTION TO EXTRACT READ SEQUENCES

# In[ ]:


def fasta_extraction(fasta_file_name):
        """
        Description: takes fasta file and extracts the read sequence and its sequence ID
        
        Arguments:
            fasta_file: the name of the fasta file. 
            eg. 'transcripts.fasta'
        
        Returns:
           seq:  the list of read sequences.
           seq_id: the list of sequence IDs corresponding to each read in seq.
           
        Raises:
           TypeError: If the input is not a string
           Exception: If the input is not in fasta format.
           i.e. >header
                ATGCATGCATGC
                ATGCATGCATGC
                ...
        """
        #input needs to be a str
        if not isinstance(fasta_file_name,str):
            raise TypeError('Input is not a string')
            
        seq=[] #empty list for sequences of reads
        seq_id=[] #empty list for sequence IDs of reads
        
        #opens the file from argument input
        with open ('%s' %fasta_file_name,'r') as fasta_file:
            #reads the content of the fasta file 
            fasta_file=fasta_file.read()
            
            #splits the file into individual reads
            fasta_file=fasta_file.split('>')
            
            #for every read it extracts the sequence and sequence ID
            for i in range(1,len(fasta_file)):
                line=fasta_file[i].split('\n')
                read_id=line[0]
                sequence=line[1:len(line)]
                sequence_list=''.join(sequence)
                seq.append(sequence_list)
                seq_id.append(read_id)
                
                print(str(i),"DONE!",read_id)
                
            #if input is not in fasta format, no output expected 
            if (len(seq)==0):
                raise Exception('Input not in Fasta format')
                
            return seq, seq_id


# 1.2 USER INPUT OF FASTA FILE

# In[ ]:


#user input for fasta file
user_input= False
while not user_input in ['q','quit'] :
    #get current directory information
    import os
    path = os.getcwd() 
    
    #input name of the fasta file or change directory
    name_file=input("Prompt: Provide the name of the fasta file in %s \n or type '0' to change directory: " %path)
    
    #change to another directory, if not in the correct one and then extract sequence information
    if (name_file==str(0)):
        try:
            new_dir=input("Prompt: Input new directory: ")
            os.chdir('%s' %new_dir)
            path=os.getcwd()
            name_file=input("Prompt: Name of fasta file in %s \nor type 0 to change directory: " %path) 
            sequence=fasta_extraction(name_file)
            
        #raises error if the input directory is incorrect
        except FileNotFoundError:
            print('No such file or directory: %s' %new_dir)
            
    #if in the correct directory, extracts sequence information
    else:
        try:
            sequence=fasta_extraction(name_file)
            
        #raises error if the input directory is incorrect
        except FileNotFoundError:
            print('No such file or directory: %s' %name_file)
            
    #ask to quit or restart the function
    user_input=input("Enter 'q'/'quit' to exit operation \nor 'any integer' to restart: ") 


# ## 2. MOTIF INPUT

# 2.1 FUNCTION TO DETECT MOTIF AND ITS OVERLAP

# In[ ]:


def motif_overlap(seq,motif_input,overlap):
    """
    Description: searches a motif through a list of sequences and checks for overlap
    using as it slides through the sequence.
    
    Arguments: 
        seq: a list of sequences in string type.
        motif_input: a string of the motif sequence.
        overlap: the minimum overlap that counts the alignment of the motif to the
        sequence as a match.
    
    Returns:
        position_match: a list of tuples giving information of the match detected.
        eg. (a,b,c,d,e) a-> sequence position in argument seq
                        b-> start position where the motif aligns on the sequence 
                        c-> end position where the motif aligns on the sequence
                        d-> overlap % for the match
                        e-> total match
    
    Raises:
        TypeError: Incorrect argument
        ValueError: Input of motif_input and seq is empty or overlap % is below 0%
    """
    if not isinstance(motif_input, str):
        raise TypeError('Incorrect argument type: motif_input')
    if not isinstance(overlap, int):
        raise TypeError('Incorrect argument type: overlap')
    if len(motif_input) == 0 or len(seq) == 0:
        raise ValueError('length of motif and/or read must be bigger than 0')
    if overlap < 0:
        raise ValueError('overlap must be longer than 0')
    
    #dictionary of IUPAC code for consensus sequences 
    IUPAC_code= {'A':'A', 'C':'C', 'G':'G', 'T':'T',
             'R': ['A','G'], 'Y': ['C','T'], 'S': ['G','C'], 
             'W': ['A','T'], 'K': ['G','T'], 'M': ['A','C'],
             'B': ['C','G','T'], 'D': ['A','G','T'],
             'H': ['A','C','T'], 'V': ['A','C','G'],
             'N': ['A','C','G','T'], '.':['']
             }
    
    #calculates the minimum number of bases that need to match for min. overlap criteria
    min_overlap=round((overlap/100)*len(motif_input))  
    position_match=[] #empty list to store the tuples output. eg. (position,overlap)
    match_start=0 #count of the overlaps
    
    #for every read sequence in the seq argument
    for each_read in range(0,len(seq)):
        
        #for every position in the motif 
        for every_pos in range(0,(len(motif_input)-min_overlap+1)):
            
            #motif aligned at the read sequence begins with the minimum overlaping bases required
            #till the whole motif sequence is aligned with the start point of the read sequence
            motif_input_start=motif_input[abs((len(motif_input)-min_overlap)-every_pos):len(motif_input)]
            
            #align the motif to the read sequence
            sequence_start=seq[each_read][0:min_overlap+every_pos]
            
            #only proceeds if the motif length is greater than or equal to min overlap 
            #and is less than and equal to full length of motif 
            #so that it only slides till the whole motif is aligned with the sequence
            if len(motif_input_start)>=min_overlap and len(motif_input_start)<=len(motif_input):
                
                #count for the particular motif length aligned
                count=0
                
                #check for match at every motif position
                for position in range(0,len(motif_input_start)):
                    
                    #retrieve information on the nucleotides permissable at the position
                    IUPAC_motif=IUPAC_code[motif_input_start[position]]
                    
                    #increases count if match
                    if (motif_input_start[position]==sequence_start[position]) or (sequence_start[position] in IUPAC_motif):
                        count=count+1      
                    else:
                        count=count #no change
                        
                #for motif length aligned at the sequence read, if the base matches is greater than
                #min. overlap requirement, then a match is confirmed and the position and percent overlap stored
                if count>=min_overlap:
                    #counts the number of matches
                    match_start=match_start+1
                    
                    #end position
                    position_match_start=min_overlap+every_pos
                    
                    #% overlap
                    overlap_per_start=(count/len(motif_input))*100
                    
                    #read sequence position, start position, end position, % overlap, total matches till that read
                    info_match_start= each_read,position_match_start-len(motif_input_start),position_match_start,overlap_per_start,match_start
                    position_match.append(info_match_start)
                    
                #print(motif_input_start) in case to visualise the alignment
                #print(sequence_start)
                    
            else:
                break
                
        #for every position in the read sequence
        for every_pos in range(1,len(seq[each_read])):
            
            #the whole motif aligned to the read sequence and slides
            #till it reaches the end of the read sequence from where it left off
            motif_start_point=len(motif_input)
            sequence_center=seq[each_read][0+every_pos:motif_start_point+every_pos]
            
            #proceed only if the whole motif is considered
            if len(sequence_center)==len(motif_input):
                
                #count for the particular motif length aligned
                count=0
                
                #check for match at every motif position
                for position in range(0,len(motif_input)):
                    
                    #retrieve information on the nucleotides permissable at the position
                    IUPAC_motif=IUPAC_code[motif_input[position]]
                    if (motif_input[position]==sequence_center[position]) or (sequence_center[position] in IUPAC_motif):
                        
                        #increases count if match
                        count=count+1
                    else:
                        count=count 
                
                #for motif length aligned at the sequence read, if the base matches is greater than
                #min. overlap requirement, then a match is confirmed and the position and percent overlap stored
                if count>=min_overlap:
                    #counts the number of matches
                    match_start=match_start+1
                    
                    #end position
                    position_match_centre=motif_start_point+every_pos
                    
                    #% overlap
                    overlap_per_centre=(count/len(motif_input))*100
                    
                    #read sequence position, start position, end position, % overlap, total matches till that read
                    info_match_centre= each_read,0+every_pos,position_match_centre,overlap_per_centre,match_start
                    position_match.append(info_match_centre)
                    
                    #print("Match! Number:",match)
                #print(motif_input) in case to visualise the alignment
                #print(sequence_center)
                    
            else:
                break 
                
        #for every position in motif        
        for every_pos in range(1,min_overlap):
            
            
            motif_input_end=motif_input[0:len(motif_input)-every_pos]
            sequence_end=seq[each_read][len(seq[each_read])-len(motif_input)+every_pos:len(seq[each_read])]
            #only proceeds if the motif length is greater than or equal to min overlap 
            #and is less than full length of motif 
            if len(motif_input_end)>=min_overlap and len(motif_input_end)<=len(motif_input):
                #count for the particular motif length aligned
                count=0
                
                #check for match at every motif position
                for position in range(0,len(motif_input_end)):
                    
                    #retrieve information on the nucleotides permissable at the position
                    IUPAC_motif=IUPAC_code[motif_input_end[position]]
                    if (motif_input_end[position]==sequence_end[position]) or (sequence_end[position] in IUPAC_motif):
                        
                        #increases count if match
                        count=count+1
                    else:
                        count=count 
                        
                #for motif length aligned at the sequence read, if the base matches is greater than
                #min. overlap requirement, then a match is confirmed and the position and percent overlap stored
                if count>=min_overlap:
                    #counts the number of matches
                    match_start=match_start+1
                    
                    #end position
                    position_match_end=len(seq[each_read])
                    
                    #% overlap
                    overlap_end=(count/len(motif_input))*100
                    
                    #read sequence position, start position, end position, % overlap, total matches till that read
                    info_match_end=each_read,len(seq[each_read])-len(motif_input)+every_pos,position_match_end,overlap_end,match_start
                    position_match.append(info_match_end)
                    
            
                    #print("Match! Number:",match)
                #print(motif_input_end) in case to visualise the alignment
                #print(sequence_end)       
            else:
                break
        print(position_match)
    return position_match


# 2.2 USER INPUT FOR MOTIF AND OVERLAP

# In[ ]:


#user input for MOTIF SEARCH 
user_input= False
while not user_input in ['q','quit'] :
    start_user=input("Prompt: Start motif search type '1':")
    
    if (start_user==str(1)):
        #user input of motif
        motif_input=input("Prompt: please enter motif to be searched \n(consensus seq. should be according to IUPAC_code):")

        #user input of minimum overlap desired 
        overlap=int(input("Prompt: please enter min. overlap value (%):"))
        
        motif_overlap_search=motif_overlap(sequence[0],motif_input,overlap)
            
    else:
         print("Type 1 to start search")   
    #ask to quit or restart the function
    user_input=input("Enter 'q'/'quit' to exit operation \nor 'any integer' to restart: ") 


# ## 3.0 TEST

# 1. Download the transcripts.fasta file
# 2. Run the code and follow 
# 3. Done

# the end.
