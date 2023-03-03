# Cancer_Immunotherapy
Challenge Overview:
This competition centers on the question of how to make T cells, the fighter cells of our immune system, better at killing cancer cells. The Eric and Wendy Schmidt Center at the Broad Institute of MIT and Harvard, in conjunction with Harvard’s Laboratory for Innovation Science and other collaborators, are hosting this competition.
While scientists have tested some genetic modifications, or “perturbations,” to T cells in the lab, there are too many possible perturbations — and combinations of perturbations — to solve this problem experimentally. Our goal in this competition is to bring machine learning to this problem to help identify which perturbations could make cancer immunotherapy more effective.
Challenge Requirements:
Predict the cell-state proportions of the held-out genes: For each of the 7 held-out knockouts (targeting genes 'Aqr', 'Bach2', 'Bhlhe40', 'Ets1', 'Fosb', 'Mafk', 'Stat3'), predict the resulting 5-dimensional vector of cell state proportions (a,b,c,d,e)(a,b,c,d,e), where
a=predicted proportion of progenitor cellsb=predicted proportion of effector cellsc=predicted proportion of terminal exhausted cellsd=predicted proportion of cycling cellse=predicted proportion of other cells
abcde=predicted proportion of progenitor cells=predicted proportion of effector cells=predicted proportion of terminal exhausted cells=predicted proportion of cycling cells=predicted proportion of other cells​
Cell state proportions should add to 1: the 5-dimensional vector of cell state proportions (a,b,c,d,e)(a,b,c,d,e) should be such that a+b+c+d+e=1a+b+c+d+e=1.
Validation-test split: The 7 held-out genes will be organised into 3 genes in the validation set (targeting genes 'Aqr', 'Bach2', 'Bhlhe40'), and 4 in the test set (targeting genes 'Ets1', 'Fosb', 'Mafk', 'Stat3').
Write-up: In addition, provide a short maximum 1-page write-up describing your approach.
More details on the submission structure and specific output requirements are shared in the relevant sections below.
Scoring:
The predicted cell state proportion vector for each of the 7 held-out knockouts will be evaluated based on the total variation distance
i.e., the l1l1​-loss, to the experimentally determined ground truth proportion vector.
For example, if the ground truth cell state proportion vector is (0.2,0.1,0.3,0,0.4)(0.2,0.1,0.3,0,0.4) and a submission's predicted cell state proportion vector is (0.19,0.11,0.28,0,0.42)(0.19,0.11,0.28,0,0.42), then the loss is ∣0.2−0.19∣+∣0.1−0.11∣+∣0.3−0.28∣+∣0−0∣+∣0.4−0.42∣=0.06.∣0.2−0.19∣+∣0.1−0.11∣+∣0.3−0.28∣+∣0−0∣+∣0.4−0.42∣=0.06.
