# Phage pipeline repo
Contains pipeline going from genome data (with Phanotate annotation) through mmseqs2 clustering and
profile construction to all vs all HMM comparison.

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
