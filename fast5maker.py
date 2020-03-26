from ont_fast5_api import fast5_file
from ont_fast5_api.analysis_tools import basecall_1d,event_detection
from Bio import SeqIO
import sys
import datetime

a=[i for i in SeqIO.parse(sys.argv[1],'fasta')][:2]

f5_file = fast5_file.Fast5File('test.fast5' , 'w')
f5_file.add_read(1,'qwerdfsdv1',1,1,1,1)
latest_analysis = f5_file.get_latest_analysis('Basecall_1D' , increment=True)
latest_analysis2 = f5_file.get_latest_analysis('EventDetection' , increment=True)
f5_file.add_analysis('basecall_1d',latest_analysis,{})
f5_file.add_analysis('event_detection',latest_analysis2,{'Reads':'1'})

b1d = basecall_1d.Basecall1DTools(f5_file,mode='w',group_name=latest_analysis)

b1ev = event_detection.EventDetectionTools(f5_file,mode='w',group_name=latest_analysis2)
b1ev.set_event_data([1],{'read_number':1})

if a[0]:
    seq=a[0].seq
    b1d.add_called_sequence('template', 'read1', seq, len(seq)*'#')
#if a[1]:
#    seq=a[1].seq
#    b1d.add_called_sequence('complement', 'read1', seq, len(seq)*'#')
