#!/bin/sh




sudo mkfs.ext4 /dev/md0
sudo mount /dev/md0 /home/ec2-user/scratch
sudo chmod -R 0775 /home/ec2-user/scratch
sudo chown ec2-user:ec2-user /home/ec2-user/scratch
