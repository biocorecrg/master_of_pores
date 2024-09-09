//batch.tf

variable "ami_batch" {
  type = string
}

variable "ami_batch_gpu" {
  type = string
}

variable "bid_percentage" {
  type    = number
  default = 50
}

variable "bid_percentage_gpu" {
  type    = number
  default = 50
}

variable "max_vcpus" {
  type    = number
  default = 16
}

variable "min_vcpus" {
  type    = number
  default = 0
}

variable "desired_vcpus" {
  type    = number
  default = 0
}

variable "instance_batch" {
  type    = list(string)
  default = ["optimal"]
}

variable "instance_batch_gpu" {
  type    = list(string)
  default = ["p3"]
}


variable "compute_environment_name" {
  type    = string
  default = "nf-compute-spot"
}

variable "compute_environment_type" {
  type    = string
  default = "SPOT"
}

variable "compute_environment_name_gpu" {
  type    = string
  default = "nf-compute-spot-gpu"
}

variable "compute_environment_type_gpu" {
  type    = string
  default = "SPOT"
}


variable "subnets" {
  type    = list(string)
  default = ["subnet-8a280df7", "subnet-c54d6588", "subnet-b85ab5d2"]
}

variable "queue_name" {
  type    = string
  default = "spot"
}

variable "queue_name_gpu" {
  type    = string
  default = "spot-gpu"
}


resource "aws_batch_compute_environment" "nf-compute-spot" {

  compute_environment_name = format("%s-%s", var.compute_environment_name, random_string.rand.result)

  compute_resources {
    instance_role = aws_iam_instance_profile.ComputeInstanceProfile.arn

    image_id = var.ami_batch

    max_vcpus     = var.max_vcpus
    min_vcpus     = var.min_vcpus
    desired_vcpus = var.desired_vcpus

    instance_type = var.instance_batch

    subnets = var.subnets

    type = var.compute_environment_type

    spot_iam_fleet_role = (var.compute_environment_type == "SPOT" ? aws_iam_role.ClusterFleetRole.arn : null)

    bid_percentage = (var.compute_environment_type == "SPOT" ? var.bid_percentage : null)

    security_group_ids = [aws_security_group.allow_all.id]

  }

  service_role = aws_iam_role.ClusterRole.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_policy_attachment.AWSBatchServiceRole-policy-attachment]

  tags = {
    name = "nf-cluster"
  }
}

resource "aws_batch_compute_environment" "nf-compute-spot-gpu" {

  compute_environment_name = format("%s-%s", var.compute_environment_name_gpu, random_string.rand.result)

  compute_resources {
    instance_role = aws_iam_instance_profile.ComputeInstanceProfile.arn

    image_id = var.ami_batch_gpu

    max_vcpus     = var.max_vcpus
    min_vcpus     = var.min_vcpus
    desired_vcpus = var.desired_vcpus

    instance_type = var.instance_batch_gpu

    subnets = var.subnets

    type = var.compute_environment_type_gpu

    spot_iam_fleet_role = (var.compute_environment_type_gpu == "SPOT" ? aws_iam_role.ClusterFleetRole.arn : null)

    bid_percentage = (var.compute_environment_type_gpu == "SPOT" ? var.bid_percentage_gpu : null)


    security_group_ids = [aws_security_group.allow_all.id]

  }

  service_role = aws_iam_role.ClusterRole.arn
  type         = "MANAGED"
  depends_on   = [aws_iam_policy_attachment.AWSBatchServiceRole-policy-attachment]

  tags = {
    name = "nf-cluster-gpu"
  }
}

resource "aws_batch_job_queue" "spot" {

  name                 = var.queue_name
  state                = "ENABLED"
  priority             = 1
  compute_environments = [aws_batch_compute_environment.nf-compute-spot.arn]

  depends_on = [aws_batch_compute_environment.nf-compute-spot]

  tags = {
    name = "nf-queue"
  }
}

resource "aws_batch_job_queue" "spot-gpu" {

  name                 = var.queue_name_gpu
  state                = "ENABLED"
  priority             = 1
  compute_environments = [aws_batch_compute_environment.nf-compute-spot-gpu.arn]

  depends_on = [aws_batch_compute_environment.nf-compute-spot-gpu]
  tags = {
    name = "nf-queue-gpu"
  }
}

