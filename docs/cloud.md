---
layout: page
title: Running on the cloud 
navigation: 5
---

## Amazon Web Service EC2

The simplest option is running an EC2 instance interactively where the pipeline can be installed as pointed in the previous documentation pages.

Last available Amazon Machine (AMI) we provide is:
* **ami-06da745a6b33715f1** (Ubuntu 18.04, CUDA compatible, Docker 19.x and Singularity 3.2.1 preinstalled)

### Share files in Amazon S3

Amazon Simple Storage Service (S3) is a convenient web service storage system for sharing raw input and final output files between your premises and your computing cloud instances.

Below we provide some instructions and advices to set up a S3 bucket.

Some convenient instructions for S3 permisions in your EC2 instance can be found [here](https://cloudkul.com/blog/mounting-s3-bucket-linux-ec2-instance/). From the previous link you can learn how to retrieve the key and password to place in ```/root/.passwd-s3fs```

Ensure proper permission as well: ```chmod 600 /root/.passwd-s3fs```

You can include in ```/etc/fstab``` the following mounting point (adapt according to your case):

    frankfurt-nf    /mnt/frankfurt-nf fuse.s3fs _netdev,allow_other,passwd_file=/root/.passwd-s3fs,uid=1000,gid=1000   0 0


### Terraform

Place [terraform](https://www.terraform.io/downloads.html) binary in your path
    
    terraform init
    terraform validate
    terraform plan
    terraform apply
    

Connect to your EC2 instance:

    ssh -i "key-nf.pem" ubuntu@xxx.eu-central-1.compute.amazonaws.com