

nucmer --prefix=ecoli_vs_ref ec_reference.fasta ./dastool_output/SemiBin_1.fa
dnadiff -d ecoli_vs_ref.delta
mummerplot --png ecoli_vs_ref.delta  # optional visualization


conda create -n quast
conda activate quest
conda install conda-forge::zlib
conda install pip 
pip install quast
pip install setuptools

quast.py ./dastool_output/SemiBin_1.fa -r ec_reference.fasta -o quast_output
