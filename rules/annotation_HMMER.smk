#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #

hmmscan --cpu 2 --domtblout TrinotatePFAM.out \
          ~/CourseData/RNA_data/trinity_trinotate_tutorial/trinotate_data/Pfam-A.hmm \
          Trinity.fasta.transdecoder.pep


# ----- rule ----- #