#GTA_Hunter
  It is a bioinformatics tool that classifies homologs from GTA structural gene as true GTAs or viruses. The tool is written in `Python v3.5.1` and requires following Python packages.

- NumPy v1.11
- CVXOPT v1.1.8

**GTA_Hunter.py is the main file, and is used through command line using the following commands:**
```
usage: GTA_Hunter.py [-h] -g GTA -v VIRUS [-q QUERIES] [-k [KMER]]
                     [-p [PSEAAC]] [-y] [-m] [-w WEIGHT WEIGHT]
                     [-t CLUSTER_TYPE] [-d DIST] [-c C] [-x [XVAL]]
                     [-e KERNEL KERNEL] [-s]

Gene Classification Using SVM.

optional arguments:
  -h, --help            show this help message and exit
  -g GTA, --GTA GTA     The .faa or .fna training file for GTA genes.
  -v VIRUS, --virus VIRUS
                        The .faa or .fna training file for viral genes.
  -q QUERIES, --queries QUERIES
                        The .faa or .fna query file to be classified.
  -k [KMER], --kmer [KMER]
                        The kmer size needed for feature generation
                        (default=4).
  -p [PSEAAC], --pseaac [PSEAAC]
                        Expand feature set to include pseudo amino acid
                        composition. Specify lamba (default=3). Weight = 0.05.
  -y, --physico         Expand feature set to include physicochemical
                        composition.
  -m, --min             Print bare minimum results.
  -w WEIGHT WEIGHT, --weight WEIGHT WEIGHT
                        Allows for weighting of training set. Will need to
                        specify the two pairwise distance files needed for
                        weighting (GTA first, then virus).
  -t CLUSTER_TYPE, --cluster_type CLUSTER_TYPE
                        Specify 'farthest' or 'nearest' neighbors clustering
                        (default='farthest').
  -d DIST, --dist DIST  Specify the cutoff distance for clustering in the
                        weighting scheme (default=0.01).
  -c C, --soft_margin C
                        The soft margin for the SVM (default=1.0).
  -x [XVAL], --xval [XVAL]
                        Performs cross validation of training set. Specify
                        folds over 10 repetitions (default=5).
  -e KERNEL KERNEL, --kernel KERNEL KERNEL
                        Specify kernel to be used and sigma if applicable
                        (i.e. gaussian) (default='linear', 0).
  -s, --svs             Show support vectors.
```
**Weight.py may also be used separately for cluster analysis:**
```
usage: Weight.py [-h] -p PROFILE_PATH -w PAIRWISE_PATH [-t CLUSTER_TYPE]
                 [-d CUTOFF] [-v] [-i] [-q] [-r R] [-c]

optional arguments:
  -h, --help            show this help message and exit
  -p PROFILE_PATH, --profile_path PROFILE_PATH
                        The .faa or .fna file used in calculating pairwise
                        distances
  -w PAIRWISE_PATH, --pairwise_path PAIRWISE_PATH
                        The .dist file of pairwise distances created by RAxML
  -t CLUSTER_TYPE, --cluter_type CLUSTER_TYPE
                        Specify 'farthest' or 'nearest' neighbors clustering.
  -d CUTOFF, --cutoff_distance CUTOFF
                        Specify the cutoff distance for clustering.
  -v, --visualize       Shows visualization of the clustered profiles.
  -i, --histogram       Shows histogram of pairwise distances.
  -q, --quartiles       Shows quartiles of pairwise distances.
  -r R, --threshold R   Show number of pairwise distances that fall below the
                        given threshold.
  -c, --count           Shows how many profiles clustered into how many
                        clusters.
```
