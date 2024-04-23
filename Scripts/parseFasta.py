

from Bio import SeqIO

fileid=open("../Input/id_list.txt", "w")

for record in SeqIO.parse("../Input/input.txt", "fasta"):
       
        pid=record.id.strip()
        fasta=record.seq
        # print("Protein ID: ", pid)
        # print("Protein Sequence: ", fasta)
        # Write them into file
        filefasts=open("../Input/FASTA/"+pid+".fasta", "w")
        fileid.write(pid+"\n")
        filefasts.write(">"+pid+"\n")
        filefasts.write(str(fasta)+"\n")
        filefasts.close()
fileid.close()
