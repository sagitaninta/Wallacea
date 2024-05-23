OUTPUT=`echo $(basename $1) | sed 's/.glf.gz/.post/'`

python posterior_estimator.py priors_scores_db.tsv $1 ./post/$OUTPUT
