#!/bin/sh



aws s3 sync --no-progress --only-show-errors s3://diversigen-virmap-databases/kraken2Custom /home/ec2-user/scratch/kraken2Custom &
aws s3 sync --no-progress --only-show-errors s3://diversigen-virmap-databases/VirmapDatabases.04.22.2019/ /home/ec2-user/scratch/VirmapDatabases.04.22.2019 &
aws s3 sync --no-progress --only-show-errors s3://diversigen-virmap-databases/hg38phix /home/ec2-user/scratch/hg38phix &
wait
