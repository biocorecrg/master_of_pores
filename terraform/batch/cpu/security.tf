//security.tf
resource "aws_security_group" "allow_ssh" {

	name = "allow_ssh"
	//vpc_id = "${aws_vpc.nf-env.id}"
	ingress {
		cidr_blocks = [
			  "0.0.0.0/0"
		]
		from_port = 22
		to_port = 22
		protocol = "tcp"
	}
	
	// Terraform removes the default rule
	egress {
		from_port = 0
		to_port = 0
		protocol = "-1"
		cidr_blocks = ["0.0.0.0/0"]
	}
}


resource "aws_iam_policy_attachment" "AWSBatchServiceRole-policy-attachment" {

	name       = "AWSBatchServiceRole-policy-attachment"
	policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
	groups     = []
	users      = []
	roles      = ["S3access", "AWSBatchServiceRole"]

}



