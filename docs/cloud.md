---
layout: page
title: Running on the cloud 
navigation: 5
---

## Amazon Web Service EC2

The simplest option is running an EC2 instance interactively where the pipeline can be installed as pointed in the previous documentation pages.

Last available Amazon Machine (AMI) we provide is:
* **ami-0034b83deda308802** (Ubuntu 18.04, CUDA compatible, Docker 19.x and Singularity 3.2.1 preinstalled)

### Share files in Amazon S3

Amazon Simple Storage Service (S3) is a convenient web service storage system for sharing raw input and final output files between your premises and your computing cloud instances.

Below we provide some instructions and advices to set up a S3 bucket.

You can include in /etc/fstab the following mounting point. Adapt according to your instance:

    frankfurt-nf    /mnt/frankfurt-nf fuse.s3fs _netdev,allow_other,passwd_file=/root/.passwd-s3fs,uid=1000,gid=1000   0 0
