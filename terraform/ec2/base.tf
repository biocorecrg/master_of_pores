
resource "aws_instance" "computing" {
 
  ami         = "ami-06da745a6b33715f1"
  instance_type = "t2.micro"
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
  force_destroy = true
  
  tags = {
    name = "S3 Frankfurt"
  }
}

