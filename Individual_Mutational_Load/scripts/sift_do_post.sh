OUTPUT=`echo $(basename $1) | sed 's/.glf.gz/.post/'`

python sift_post_estimator.py sift_priors_scores_db.tsv $1 ./sift/post/$OUTPUT
