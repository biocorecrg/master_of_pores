
resource "aws_instance" "computing" {
 
  ami         = "ami-0bf3a9a6cb7a5ea9f"
  instance_type = "m4.2xlarge" // We suggest p3.2xlarge for GPU processing
  iam_instance_profile = "S3access"
  key_name = "key-nf"
  security_groups = [ "allow_ssh" ]
  tags = {
      name = "computing"
       
  }

}

resource "aws_s3_bucket" "frankfurt-nf" {
  bucket = "frankfurt-nf"
  acl    = "private"
  force_destroy = true // if we don't want files to be removed once tf destroy, change to false
  
  tags = {
    name = "S3 Frankfurt"
  }
}

