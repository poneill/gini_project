from Header import Header
from SequenceData import SequenceData
from RecognizerFactory import RecognizerFactory
from estremo_utils import wc

import string

class Paragraph(object):
    def __init__(self,paragraph_lines,ignore_recognizer=False):
        self.paragraph_lines = paragraph_lines
        self.header = Header(paragraph_lines[0])
        #Paragraph can't end with empty line!
        seq_lines = [line for line in paragraph_lines[1:]
                     if line[0] in string.uppercase]
        recognizer_lines = [line for line in paragraph_lines[1:]
                                if not line[0] in string.uppercase]
        self.seq_data = map(SequenceData,seq_lines)
        if not ignore_recognizer:
            recognizer_data = [[float(val) for val in line.split(",")]
                               for line in recognizer_lines]
        else:
            recognizer_data = [[0 for i in range(len(self.seq_data[0].site))]
                               for j in range(4)]
        self.recognizer = RecognizerFactory(recognizer_data)
                                
    def toFasta(self):
        return "\n".join([">Site %s\n%s" % (i,seq.site)
                          for i,seq in enumerate(self.seq_data)])
    def motif(self):
        return [seq.site for seq in self.seq_data]

    def adjusted_motif(self):
        def adjust(site):
            if self.recognizer.score(site) < self.recognizer.score(wc(site)):
                return site
            else:
                return wc(site)
        return [adjust(site) for site in self.motif()]

