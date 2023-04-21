#This is the solution for q1. 
#It will print Hello Bioinformatics
def HelloBI():
    print("Hello Bioinformatics")
#HelloBI()

#This is the answer for q3
#It will print out each header of sequences
with open('datafile.txt', 'r') as f:
    seq = ""
    for line in f:
        #The header line always starts with >
        if line.startswith('>'):
            if seq != "":
                print(header + " Length: " + str(len(seq)))
            header = line.strip()
            seq = ""
        else:
            seq += line.strip()
    print(header + " Length: " + str(len(seq)))

with open('datafile.txt', 'r') as f, open('mouse_rat_sequences.fasta', 'w') as out:
    seq = ""
    for line in f:
        if line.startswith('>'):
            if "mus" in line.lower() or "rat" in line.lower():
                if seq != "":
                    out.write(header + "\n")
                    out.write("\n".join([seq[i:i+60] for i in range(0, len(seq), 60)]) + "\n")
                header = line.strip()
                seq = ""
        else:
            seq += line.strip()
    if "mus" in header.lower() or "rat" in header.lower():
        out.write(header + "\n")
        out.write("\n".join([seq[i:i+60] for i in range(0, len(seq), 60)]) + "\n")




with open('datafile.txt', 'r') as f, open('data.seq', 'w') as seq_out:
    if_firstline= False
   
    
    for line in f:
        if line.startswith('>'):
            if if_firstline:
                seq_out.write("@")
           
            if_firstline=True
            
        else:
            seq_out.write(line.strip())

with open('datafile.txt', 'r') as f, open('data.in', 'w') as in_out:
  
    offset = 0
    
    for line in f:
        if line.startswith('>'):
            
            gi_number = line.strip().split("|")[1]
            in_out.write(gi_number + " " + str(offset) + "\n")
        else:
           
            
            offset += len(line.strip())


offsets = {}
with open('data.in', 'r') as f:
    for line in f:
        gi_number, offset = line.strip().split()
        offsets[gi_number] = int(offset)


query = "MHIQITDFGTAKVLSPDS"


with open('data.seq', 'r') as f:
    seq = f.read()
sorted_offsets = sorted(offsets.items(), key=lambda x: x[1])
##sorted_offsets+=["init",0]
for ind, (gi_number, offset) in enumerate(sorted_offsets[1:]):
    if query in seq[sorted_offsets[ind][1]:sorted_offsets[ind+1][1]]:
        print(sorted_offsets[ind][0])
       