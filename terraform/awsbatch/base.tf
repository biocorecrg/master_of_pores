// base.tf

variable "profile" {
  type    = string
  default = "default"
}

variable "credentials" {
  type = string
}

variable "region" {
  type = string
}

variable "key_name" {
  type = string
}

variable "ami_entrypoint" {
  type = string
}

variable "ec2_password" {
  type = string
}

variable "instance_type" {
  type    = string
  default = "t2.micro"
}

variable "instance_volume_size" {
  type    = number
  default = 10
}

variable "destroy_bucket" {
  type    = bool
  default = true
}

variable "bucket_acl" {
  type    = string
  default = "private"
}

variable "instance_count" {
  type    = number
  default = 2
}

variable "repourl" {
  type    = string
  default = "https://github.com/biocorecrg/MOP2"
}

variable "bucket_prefix" {
  type    = string
  default = "class-bucket"
}


provider "aws" {
  profile                 = var.profile
  shared_credentials_file = var.credentials
  region                  = var.region
}

// Random resource for naming
resource "random_string" "rand" {
  length  = 8
  special = false
}

// You may define an entry point for convenience

resource "aws_instance" "classroom" {

  ami                  = var.ami_entrypoint
  count                = var.instance_count
  instance_type        = var.instance_type
  iam_instance_profile = aws_iam_instance_profile.Multiprofile.name
  key_name             = var.key_name
  security_groups      = ["allow_ssh-${random_string.rand.result}", "allow_http-${random_string.rand.result}", "allow_shiny-${random_string.rand.result}"]
  user_data            = templatefile("ec2init.sh.tpl", { region = var.region, ec2_password = var.ec2_password, bucket_acl = var.bucket_acl, bucket_prefix = var.bucket_prefix, repourl = var.repourl, rand = random_string.rand.result, count = count.index + 1 })
  root_block_device {
    volume_size = var.instance_volume_size
  }

  // We add additional sleep time for allowing creation and proper set up of image
  provisioner "local-exec" {
    command = "sleep 5"
  }

  // Let's wait all buckets to be created first. It could be even tried one by one
  depends_on = [aws_s3_bucket.class-bucket, aws_iam_instance_profile.Multiprofile]

  tags = {
    name = "classroom-${count.index + 1}"
  }

}

resource "aws_s3_bucket" "class-bucket" {
  count         = var.instance_count
  bucket        = format("%s-%s", var.bucket_prefix, count.index + 1)
  acl           = var.bucket_acl
  force_destroy = var.destroy_bucket

  tags = {
    name = format("%s-%s", var.bucket_prefix, count.index + 1)
  }
}
