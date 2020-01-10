---
layout: page
title: Running on the cloud 
navigation: 6
---

## Amazon Web Service EC2

The simplest option is running an EC2 instance interactively where the pipeline can be installed as pointed in the previous documentation pages.

Last available Amazon Machine (AMI) we provide is:
* **ami-0bf3a9a6cb7a5ea9f** (Ubuntu 18.04, CUDA compatible, Docker 19.x, Singularity 3.2.1 and Nextflow 19.10 preinstalled)

When running an instance among the different [available types](https://aws.amazon.com/ec2/instance-types/), minimum CPU and memory requirements must be taken into account. These must fit with Nextflow executor process configuration values.

Keep in mind that not all [Amazon infrastructure regions](https://aws.amazon.com/about-aws/global-infrastructure/regions_az/) may have the same instance types. As a example, in January 2020 Frankfurt has GPU nodes, but not Paris. 

Launch an instane from the AMI image above (Go to EC2 > Images > AMI and copy-paste the ID provided above filtering in public images). Once you find that image, you can launch an instance from it.

You can connect to the launched instance by using this command below:

    ssh -i "key-nf.pem" ubuntu@xxx.eu-central-1.compute.amazonaws.com
    
where ```key-nf.pem``` is your private key ([reference](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html)) and host details can be obtained from *Connect popup* in EC2 instances dashboard.

### Terraform

For sake of commodity, you may prefer to automate deployment of EC2 instances and S3 buckets. Terraform is a convenient tool for this.

Place [terraform](https://www.terraform.io/downloads.html) binary in your local workstation path and move where your are keeping your tf files. Examples are provided in the terraform directory of this repository.

Adapt terraform configuration files to include your credentials, use your chosen instance types, which key pair they are associated with, or whether allow files in S3 bucket to be kept or not (```force_destroy``` parameter).

Initialize terraform directory:

    terraform init
    
Validate terraform files:

    terraform validate
    
Inspect what changes are going to be performed in your AWS account:

    terraform plan
    
Proceed:
    
    terraform apply
    
Once analyses are finished, infrastructure can be dismantled with:

    terraform destroy
    

### Share files in Amazon S3

Amazon Simple Storage Service (S3) is a convenient web service storage system for sharing raw input and final output files between your premises and your computing cloud instances.

Below we provide some instructions and advices to set up a S3 bucket.

Some convenient instructions for S3 permisions in your EC2 instance can be found [here](https://cloudkul.com/blog/mounting-s3-bucket-linux-ec2-instance/). From the previous link you can learn how to retrieve the key and password to place in ```/root/.passwd-s3fs```

Ensure proper permission as well: ```chmod 600 /root/.passwd-s3fs```

You can include in ```/etc/fstab``` the following mounting point (adapt according to your case):

    frankfurt-nf    /mnt/frankfurt-nf fuse.s3fs _netdev,allow_other,passwd_file=/root/.passwd-s3fs,uid=1000,gid=1000   0 0

If not mounted, you can mount it therefore straightforward by running:

    sudo mount /mnt/frankfurt-nf
    
Adapt your S3 bucket and mounting point names according to your choice.

Specially for huge amount of data, we suggest to use [AWS CLI](https://aws.amazon.com/cli/) to transfer files from your premises to a S3 Bucket ([Ref](https://docs.aws.amazon.com/en_us/cli/latest/userguide/cli-services-s3.html)). For instance, the commandline below uploads the data example file in a pre-existing bucket.

    aws s3 cp  multifast5_1.fast5 s3://frankfurt-nf


    
