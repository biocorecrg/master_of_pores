export TF_VAR_key_name=key-nf
# Modify instance type to fit more needs if desired: https://aws.amazon.com/ec2/instance-types/t2/
export TF_VAR_instance_type=t2.small
# Image used for entrypoint
export TF_VAR_ami_entrypoint=ami-0cf1f74891140b374
# Image used for setting in the cluster
export TF_VAR_ami_batch=ami-06b8c6e4fe388181d
export TF_VAR_ami_batch_gpu=ami-06b8c6e4fe388181d
# Region
export TF_VAR_region=eu-central-1
export TF_VAR_profile=default
# If lower bid percentage, it will take longer to run in AWS Batch, but it will be cheaper
export TF_VAR_bid_percentage=90
export TF_VAR_credentials=/home/myuser/.aws/credentials
export TF_VAR_ec2_password=sshpassword
export TF_VAR_instance_count=1
export TF_VAR_instance_volume_size=45
export TF_VAR_bucket_acl=public-read
export TF_VAR_bucket_prefix=mop2-bucket
export TF_VAR_compute_environment_name=nf-compute
export TF_VAR_compute_environment_name_gpu=nf-compute-gpu
export TF_VAR_queue_name=mop
export TF_VAR_queue_name_gpu=mop-gpu
export TF_VAR_compute_environment_type=SPOT
export TF_VAR_compute_environment_type_gpu=EC2
export TF_VAR_instance_batch='["optimal"]'
export TF_VAR_instance_batch_gpu='["p3"]'
export TF_VAR_subnets='["subnet-8a280df7", "subnet-c54d6588", "subnet-b85ab5d2"]'
export TF_VAR_repourl=https://github.com/biocorecrg/MOP2
export AWS_ACCOUNT_ID=$(aws sts get-caller-identity|jq .Account|tr -d \")
