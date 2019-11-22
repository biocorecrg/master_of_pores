provider "aws" {
  access_key = "yourkeyhere"
  secret_key = "yoursecrethere"
  region     = "eu-central-1"
}


// You may define an entry point for convenience

resource "aws_instance" "entrypoint" {
 
  ami         = "ami-06a3c664bd6bb3fce"
  instance_type = "t2.micro"
  iam_instance_profile = "S3access"
  key_name = "key-nf"
  security_groups = [ "allow_ssh" ]
  tags = {
	name = "entrypoint"
  }

}

resource "aws_s3_bucket" "frankfurt-nf" {
  bucket = "frankfurt-nf"
  acl    = "private"
  force_destroy = true
  
  tags = {
    name = "S3 Frankfurt"
  }
}

