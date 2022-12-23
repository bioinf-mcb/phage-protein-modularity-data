# Phage pipeline repo
Contains pipeline going from protein data through mmseqs2 clustering and
profile construction, all vs all HMM comparison to pairwise coverage comparison.

### Required software:
- mmseqs2
- HHSuite (with databases)
- bgzip

### Required Python libraries (as in requirements.txt file):
- biopython==1.78
- csb==1.2.5
- notebook==6.3.0
- scikit-learn==0.24.2
- seaborn==0.11.1

### Environment
- Pipeline should be run in virtual environment created with virtualenv
- Pipeline was tested on Python3.9/MacOS and production run was made on Python3.6/Linux
