import matplotlib.pyplot as plt
import numpy as np

test_input = [{'query': 'YP_498987.1', 'qcov': 54.2, 'sstart': 349, 'send': 608, 'scov': 39.3, 'hit_length': 659}, {'query': 'YP_499313.1', 'qcov': 92.5, 'sstart': 366, 'send': 626, 'scov': 39.5, 'hit_length': 659}, {'query': 'YP_499364.1', 'qcov': 87.2, 'sstart': 417, 'send': 637, 'scov': 33.4, 'hit_length': 659}, {'query': 'YP_500562.1', 'qcov': 62.7, 'sstart': 448, 'send': 614, 'scov': 25.2, 'hit_length': 659}]

y_count = len(test_input)

plt.plot([0,test_input[0]['hit_length']],[y_count+1,y_count+1],label = "TCDB Reference" + "(0-" + str(test_input[0]['hit_length']) + ")")

for candidate in test_input:
    plt.plot([candidate['sstart'],candidate['send']], [y_count, y_count], label = candidate['query'] + "(" + str(candidate['sstart']) + "-" + str(candidate['send']) + ")")
    y_count -= 1

plt.legend()
plt.show()
