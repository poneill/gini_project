# I am actually writing a Factory function.  What have I become?
# https://twitter.com/jleedev/status/224875765588701184
from SLP import SLP
from MLP import MLP
from LinearRecognizer import LinearRecognizer

def RecognizerFactory(recognizer_data):
    """Initialize the recognizer.  Recognizer data may describe
        either an SLP or an MLP, so determine which and represent it
        appropriately."""
    if LinearRecognizer.validate_input(recognizer_data):
        return LinearRecognizer(recognizer_data)
    elif SLP.validate_input(recognizer_data): #Should never eval to true
        return SLP(recognizer_data)
    elif MLP.validate_input(recognizer_data):
        return MLP(recognizer_data)
    else:
        raise Exception("Couldn't parse recognizer data")

# def test():
# from matplotlib import pyplot as plt
# for i in range(1,8+1):
#     print "starting on:",i
#     with open("../correlation-%s.csv" % i) as f:
#  	paragraph = chomp_paragraph(f,False)
#  	new_paragraph = chomp_paragraph(f,False)
#  	while new_paragraph:
#  		paragraph = new_paragraph
# 		new_paragraph = chomp_paragraph(f,False)
#     exec("mlp%s = paragraph.recognizer" % i)
#     exec("sites%s = paragraph.motif()" % i)
# mlps = [eval("mlp%s" % i) for i in range(1,8+1)]
# siteses = [eval("sites%s" % i) for i in range(1,8+1)]
# for (mlp,sites) in zip(mlps,siteses):
#     for site in sites:
#         mlp.detect_correlations2(site)
    
